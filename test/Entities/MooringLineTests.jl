import Mooring.ParameterHandler as PH
import Mooring.MooringLines as ML
import Mooring.MooringSegments as Seg
using Gridap.MultiField
using Gridap.Geometry: get_node_coordinates
using Gridap.TensorValues
using Test

@testset "MooringLine Tests" begin

ph = PH.ParameterHandler()

# Minimal topology
ph.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
ph.points[2] = PH.PointParameters(id=2, coords=[1.0, 0.0], mesh_size=0.2)
ph.segments[1] = PH.SegmentParameters(id=1, tag="line1", start_point=1, stop_point=2, length=1.5)
ph.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])

# Setup lines
mlines = ML.setup_lines(ph)

@test length(mlines) == 1
line = mlines[1]
@test line isa ML.MooringLine

segments = ML.get_segments(line)
@test length(segments) == 1
segment = segments[1]
@test segment isa Seg.MooringSegment
@test segment.tag == "line1"

points = Seg.get_points(segments[1])
start_point = points[ph.segments[1].start_point]
stop_point = points[ph.segments[1].stop_point]

@testset "Linear Map Tests" begin
    f_linear = ML.get_physical_linear_map(ph.segments[1], ph, start_point, stop_point)
    
    # Get reference coordinates
    p1_ref_node = start_point.btrian.glue.face_to_bgface[1]
    p2_ref_node = stop_point.btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(start_point.btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(stop_point.btrian)[p2_ref_node][1]
    
    # Test at endpoints
    start_pos = f_linear(VectorValue(x_ref_p1))
    @test isapprox(start_pos[1], 0.0; atol=1e-8)
    @test isapprox(start_pos[2], 0.0; atol=1e-8)
    
    end_pos = f_linear(VectorValue(x_ref_p2))
    @test isapprox(end_pos[1], 1.0; atol=1e-8)
    @test isapprox(end_pos[2], 0.0; atol=1e-8)
    
    # Test at midpoint (halfway along segment)
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    mid_pos = f_linear(VectorValue(mid_ref))
    @test isapprox(mid_pos[1], 0.5; atol=1e-8)
    @test isapprox(mid_pos[2], 0.0; atol=1e-8)
end

@testset "Quadratic Map - 2D Horizontal Cable (No Sag)" begin
    # Horizontal cable with exact length - should have no sag
    ph_horiz = PH.ParameterHandler()
    ph_horiz.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
    ph_horiz.points[2] = PH.PointParameters(id=2, coords=[10.0, 0.0], mesh_size=0.1)
    ph_horiz.segments[1] = PH.SegmentParameters(id=1, tag="horizontal", start_point=1, stop_point=2, length=10.0)
    ph_horiz.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_horiz = ML.setup_lines(ph_horiz)
    segments_horiz = ML.get_segments(mlines_horiz[1])
    points_horiz = Seg.get_points(segments_horiz[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_horiz.segments[1], ph_horiz, points_horiz[1], points_horiz[2])
    f_linear = ML.get_physical_linear_map(ph_horiz.segments[1], ph_horiz, points_horiz[1], points_horiz[2])
    
    # Get reference coordinates
    p1_ref_node = points_horiz[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_horiz[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_horiz[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_horiz[2].btrian)[p2_ref_node][1]
    
    # Test at several points - should match linear map (no excess length)
    for t in [0.0, 0.25, 0.5, 0.75, 1.0]
        r = VectorValue(x_ref_p1 + t * (x_ref_p2 - x_ref_p1))
        pos_quad = f_quad(r)
        pos_linear = f_linear(r)
        @test isapprox(pos_quad[1], pos_linear[1]; atol=1e-8)
        @test isapprox(pos_quad[2], pos_linear[2]; atol=1e-8)
    end
end

@testset "Quadratic Map - 2D Vertical Hanging Cable" begin
    # Vertical cable with excess length - should have sag
    ph_vert = PH.ParameterHandler()
    ph_vert.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
    ph_vert.points[2] = PH.PointParameters(id=2, coords=[0.0, -10.0], mesh_size=0.1)
    ph_vert.segments[1] = PH.SegmentParameters(id=1, tag="vertical", start_point=1, stop_point=2, length=12.0)
    ph_vert.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_vert = ML.setup_lines(ph_vert)
    segments_vert = ML.get_segments(mlines_vert[1])
    points_vert = Seg.get_points(segments_vert[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_vert.segments[1], ph_vert, points_vert[1], points_vert[2])
    
    # Get reference coordinates
    p1_ref_node = points_vert[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_vert[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_vert[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_vert[2].btrian)[p2_ref_node][1]
    
    # Endpoints should match exactly
    start_pos = f_quad(VectorValue(x_ref_p1))
    @test isapprox(start_pos[1], 0.0; atol=1e-8)
    @test isapprox(start_pos[2], 0.0; atol=1e-8)
    
    end_pos = f_quad(VectorValue(x_ref_p2))
    @test isapprox(end_pos[1], 0.0; atol=1e-8)
    @test isapprox(end_pos[2], -10.0; atol=1e-8)
    
    # Midpoint should have sag (but x-coordinate stays at 0 for vertical cable)
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    mid_pos = f_quad(VectorValue(mid_ref))
    @test isapprox(mid_pos[1], 0.0; atol=1e-8)
    # Vertical position should be lower than -5.0 due to sag
    @test mid_pos[2] < -5.0
end

@testset "Quadratic Map - 2D Inclined Cable with Sag" begin
    # Inclined cable with excess length
    ph_incl = PH.ParameterHandler()
    ph_incl.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
    ph_incl.points[2] = PH.PointParameters(id=2, coords=[10.0, -5.0], mesh_size=0.1)
    ph_incl.segments[1] = PH.SegmentParameters(id=1, tag="inclined", start_point=1, stop_point=2, length=15.0)
    ph_incl.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_incl = ML.setup_lines(ph_incl)
    segments_incl = ML.get_segments(mlines_incl[1])
    points_incl = Seg.get_points(segments_incl[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_incl.segments[1], ph_incl, points_incl[1], points_incl[2])
    f_linear = ML.get_physical_linear_map(ph_incl.segments[1], ph_incl, points_incl[1], points_incl[2])
    
    # Get reference coordinates
    p1_ref_node = points_incl[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_incl[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_incl[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_incl[2].btrian)[p2_ref_node][1]
    
    # Test endpoints match
    start_pos = f_quad(VectorValue(x_ref_p1))
    @test isapprox(start_pos[1], 0.0; atol=1e-8)
    @test isapprox(start_pos[2], 0.0; atol=1e-8)
    
    end_pos = f_quad(VectorValue(x_ref_p2))
    @test isapprox(end_pos[1], 10.0; atol=1e-8)
    @test isapprox(end_pos[2], -5.0; atol=1e-8)
    
    # Test midpoint
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    mid_pos_quad = f_quad(VectorValue(mid_ref))
    mid_pos_linear = f_linear(VectorValue(mid_ref))
    
    # Horizontal position should match
    @test isapprox(mid_pos_quad[1], 5.0; atol=1e-8)
    @test isapprox(mid_pos_quad[1], mid_pos_linear[1]; atol=1e-8)
    
    # Vertical position should be lower (more sag)
    @test mid_pos_quad[2] < mid_pos_linear[2]
    @test mid_pos_quad[2] < -2.5  # Linear midpoint is at -2.5
    
    # Test parabolic symmetry
    quarter_ref = x_ref_p1 + 0.25 * (x_ref_p2 - x_ref_p1)
    three_quarter_ref = x_ref_p1 + 0.75 * (x_ref_p2 - x_ref_p1)
    
    quarter_pos = f_quad(VectorValue(quarter_ref))
    three_quarter_pos = f_quad(VectorValue(three_quarter_ref))
    
    # Calculate sag relative to linear interpolation
    linear_quarter = f_linear(VectorValue(quarter_ref))
    linear_three_quarter = f_linear(VectorValue(three_quarter_ref))
    
    sag_quarter = quarter_pos[2] - linear_quarter[2]
    sag_three_quarter = three_quarter_pos[2] - linear_three_quarter[2]
    
    # Sag should be symmetric
    @test isapprox(sag_quarter, sag_three_quarter; atol=1e-8)
    
    # Sag at midpoint should be maximum
    sag_mid = mid_pos_quad[2] - mid_pos_linear[2]
    @test sag_mid < sag_quarter  # More negative = more sag
end

@testset "Quadratic Map - 3D Cable" begin
    # 3D cable with excess length
    ph_3d = PH.ParameterHandler()
    ph_3d.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0, 0.0], mesh_size=0.1)
    ph_3d.points[2] = PH.PointParameters(id=2, coords=[10.0, 5.0, -8.0], mesh_size=0.1)
    ph_3d.segments[1] = PH.SegmentParameters(id=1, tag="cable_3d", start_point=1, stop_point=2, length=18.0)
    ph_3d.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_3d = ML.setup_lines(ph_3d)
    segments_3d = ML.get_segments(mlines_3d[1])
    points_3d = Seg.get_points(segments_3d[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_3d.segments[1], ph_3d, points_3d[1], points_3d[2])
    f_linear = ML.get_physical_linear_map(ph_3d.segments[1], ph_3d, points_3d[1], points_3d[2])
    
    # Get reference coordinates
    p1_ref_node = points_3d[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_3d[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_3d[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_3d[2].btrian)[p2_ref_node][1]
    
    # Test endpoints
    start_pos = f_quad(VectorValue(x_ref_p1))
    @test isapprox(start_pos[1], 0.0; atol=1e-8)
    @test isapprox(start_pos[2], 0.0; atol=1e-8)
    @test isapprox(start_pos[3], 0.0; atol=1e-8)
    
    end_pos = f_quad(VectorValue(x_ref_p2))
    @test isapprox(end_pos[1], 10.0; atol=1e-8)
    @test isapprox(end_pos[2], 5.0; atol=1e-8)
    @test isapprox(end_pos[3], -8.0; atol=1e-8)
    
    # Test midpoint
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    mid_pos_quad = f_quad(VectorValue(mid_ref))
    mid_pos_linear = f_linear(VectorValue(mid_ref))
    
    # Horizontal positions (x, y) should match linear
    @test isapprox(mid_pos_quad[1], 5.0; atol=1e-8)
    @test isapprox(mid_pos_quad[2], 2.5; atol=1e-8)
    @test isapprox(mid_pos_quad[1], mid_pos_linear[1]; atol=1e-8)
    @test isapprox(mid_pos_quad[2], mid_pos_linear[2]; atol=1e-8)
    
    # Vertical position should have sag
    @test mid_pos_quad[3] < mid_pos_linear[3]
    @test mid_pos_quad[3] < -4.0  # Linear midpoint is at -4.0
    
    # Test horizontal span calculation
    horizontal_span = sqrt(10.0^2 + 5.0^2)
    @test isapprox(horizontal_span, sqrt(125.0); atol=1e-8)
    
    # Test symmetry in 3D
    quarter_ref = x_ref_p1 + 0.25 * (x_ref_p2 - x_ref_p1)
    three_quarter_ref = x_ref_p1 + 0.75 * (x_ref_p2 - x_ref_p1)
    
    quarter_pos = f_quad(VectorValue(quarter_ref))
    three_quarter_pos = f_quad(VectorValue(three_quarter_ref))
    linear_quarter = f_linear(VectorValue(quarter_ref))
    linear_three_quarter = f_linear(VectorValue(three_quarter_ref))
    
    sag_quarter = quarter_pos[3] - linear_quarter[3]
    sag_three_quarter = three_quarter_pos[3] - linear_three_quarter[3]
    
    @test isapprox(sag_quarter, sag_three_quarter; atol=1e-8)
end

@testset "Quadratic Map - Sag Calculation Verification" begin
    # Test case where we can calculate expected sag analytically
    ph_test = PH.ParameterHandler()
    ph_test.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
    ph_test.points[2] = PH.PointParameters(id=2, coords=[100.0, 0.0], mesh_size=0.1)
    ph_test.segments[1] = PH.SegmentParameters(id=1, tag="test", start_point=1, stop_point=2, length=110.0)
    ph_test.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_test = ML.setup_lines(ph_test)
    segments_test = ML.get_segments(mlines_test[1])
    points_test = Seg.get_points(segments_test[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_test.segments[1], ph_test, points_test[1], points_test[2])
    f_linear = ML.get_physical_linear_map(ph_test.segments[1], ph_test, points_test[1], points_test[2])
    
    # Get reference coordinates
    p1_ref_node = points_test[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_test[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_test[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_test[2].btrian)[p2_ref_node][1]
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    
    mid_pos_quad = f_quad(VectorValue(mid_ref))
    mid_pos_linear = f_linear(VectorValue(mid_ref))
    
    # Calculate expected sag
    horizontal_span = 100.0
    direct_distance = 100.0
    excess_length = 110.0 - 100.0
    expected_sag_depth = excess_length^2 / (8 * horizontal_span)  # = 100 / 800 = 0.125
    
    # At midpoint (s=0.5), parabola gives: -4*h*0.5*0.5 = -h
    expected_vertical_offset = -expected_sag_depth
    
    actual_sag = mid_pos_quad[2] - mid_pos_linear[2]
    @test isapprox(actual_sag, expected_vertical_offset; atol=1e-10)
end

@testset "Quadratic Map - Edge Cases" begin
    # Test with very small excess length (should behave like linear)
    ph_small = PH.ParameterHandler()
    ph_small.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
    ph_small.points[2] = PH.PointParameters(id=2, coords=[10.0, 0.0], mesh_size=0.1)
    ph_small.segments[1] = PH.SegmentParameters(id=1, tag="small", start_point=1, stop_point=2, length=10.001)
    ph_small.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_small = ML.setup_lines(ph_small)
    segments_small = ML.get_segments(mlines_small[1])
    points_small = Seg.get_points(segments_small[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_small.segments[1], ph_small, points_small[1], points_small[2])
    f_linear = ML.get_physical_linear_map(ph_small.segments[1], ph_small, points_small[1], points_small[2])
    
    p1_ref_node = points_small[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_small[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_small[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_small[2].btrian)[p2_ref_node][1]
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    
    mid_pos_quad = f_quad(VectorValue(mid_ref))
    mid_pos_linear = f_linear(VectorValue(mid_ref))
    
    # Should be very close to linear
    @test isapprox(mid_pos_quad[1], mid_pos_linear[1]; atol=1e-8)
    @test isapprox(mid_pos_quad[2], mid_pos_linear[2]; atol=1e-6)
end

# Transient FE spaces
X, Y = ML.get_transient_FE_spaces(line)
@test X isa TransientMultiFieldFESpace
@test Y isa MultiFieldFESpace
@test length(X.spaces) == 1
@test length(Y.spaces) == 1
U = X.spaces[1]
V = Y.spaces[1]
@test U isa TransientTrialFESpace
@test V isa SingleFieldFESpace

# Reference configuration
Xₕ = ML.get_reference_configuration(line, X(0.0))
@test Xₕ isa MultiFieldFEFunction
@test length(Xₕ) == 1
X1ₕ = Xₕ[1]
@test X1ₕ isa FEFunction
@test Seg.get_triangulation(segment) == get_triangulation(X1ₕ)

# Quasi-static residual
res = ML.get_quasi_static_residual(line, Xₕ, 0.0)
@test isa(res, Function)
op = FEOperator(res, X(0.0), Y)
uₕ, = solve(op)
dΩ = Seg.get_measures(segment)[1]
unorm = √(∑(∫(uₕ⋅uₕ)dΩ))
@test unorm == 0.0

# Solve quasi-static problem
xₕ = ML.solve_quasistatic(ph)
@test length(xₕ) == 2

end # testset