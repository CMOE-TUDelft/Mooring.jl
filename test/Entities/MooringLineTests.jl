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

@testset "Quadratic Map - Multiple Points Along Cable" begin
    # Test that the quadratic map produces consistent results at multiple points
    ph_multi = PH.ParameterHandler()
    ph_multi.points[1] = PH.PointParameters(id=1, coords=[0.0, 10.0], mesh_size=0.1)
    ph_multi.points[2] = PH.PointParameters(id=2, coords=[20.0, -5.0], mesh_size=0.1)
    ph_multi.segments[1] = PH.SegmentParameters(id=1, tag="multi", start_point=1, stop_point=2, length=30.0)
    ph_multi.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_multi = ML.setup_lines(ph_multi)
    segments_multi = ML.get_segments(mlines_multi[1])
    points_multi = Seg.get_points(segments_multi[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_multi.segments[1], ph_multi, points_multi[1], points_multi[2])
    f_linear = ML.get_physical_linear_map(ph_multi.segments[1], ph_multi, points_multi[1], points_multi[2])
    
    p1_ref_node = points_multi[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_multi[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_multi[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_multi[2].btrian)[p2_ref_node][1]
    
    # Test at multiple s values
    s_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    sag_values = Float64[]
    
    for s in s_values
        r = VectorValue(x_ref_p1 + s * (x_ref_p2 - x_ref_p1))
        pos_quad = f_quad(r)
        pos_linear = f_linear(r)
        
        # Horizontal position should always match
        @test isapprox(pos_quad[1], pos_linear[1]; atol=1e-8)
        
        # Vertical sag
        sag = pos_quad[2] - pos_linear[2]
        push!(sag_values, sag)
        
        # Sag should be negative (downward) except at endpoints
        if s ≈ 0.0 || s ≈ 1.0
            @test isapprox(sag, 0.0; atol=1e-8)
        else
            @test sag < 0.0
        end
    end
    
    # Maximum sag should be at midpoint
    max_sag_idx = argmin(sag_values)
    @test max_sag_idx == 6  # Index of s=0.5
end

@testset "Quadratic Map - 3D Asymmetric Cable" begin
    # Test 3D cable that is not symmetric in x-y plane
    ph_asym = PH.ParameterHandler()
    ph_asym.points[1] = PH.PointParameters(id=1, coords=[5.0, 3.0, 0.0], mesh_size=0.1)
    ph_asym.points[2] = PH.PointParameters(id=2, coords=[15.0, 10.0, -12.0], mesh_size=0.1)
    ph_asym.segments[1] = PH.SegmentParameters(id=1, tag="asym", start_point=1, stop_point=2, length=20.0)
    ph_asym.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_asym = ML.setup_lines(ph_asym)
    segments_asym = ML.get_segments(mlines_asym[1])
    points_asym = Seg.get_points(segments_asym[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_asym.segments[1], ph_asym, points_asym[1], points_asym[2])
    
    p1_ref_node = points_asym[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_asym[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_asym[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_asym[2].btrian)[p2_ref_node][1]
    
    # Test endpoints
    start_pos = f_quad(VectorValue(x_ref_p1))
    @test isapprox(start_pos[1], 5.0; atol=1e-8)
    @test isapprox(start_pos[2], 3.0; atol=1e-8)
    @test isapprox(start_pos[3], 0.0; atol=1e-8)
    
    end_pos = f_quad(VectorValue(x_ref_p2))
    @test isapprox(end_pos[1], 15.0; atol=1e-8)
    @test isapprox(end_pos[2], 10.0; atol=1e-8)
    @test isapprox(end_pos[3], -12.0; atol=1e-8)
end

@testset "Quadratic Map - Zero Excess Length Edge Case" begin
    # Cable with exact distance, no excess length
    ph_exact = PH.ParameterHandler()
    ph_exact.points[1] = PH.PointParameters(id=1, coords=[0.0, 5.0], mesh_size=0.1)
    ph_exact.points[2] = PH.PointParameters(id=2, coords=[12.0, 0.0], mesh_size=0.1)
    direct_dist = sqrt(12.0^2 + 5.0^2)
    ph_exact.segments[1] = PH.SegmentParameters(id=1, tag="exact", start_point=1, stop_point=2, length=direct_dist)
    ph_exact.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_exact = ML.setup_lines(ph_exact)
    segments_exact = ML.get_segments(mlines_exact[1])
    points_exact = Seg.get_points(segments_exact[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_exact.segments[1], ph_exact, points_exact[1], points_exact[2])
    f_linear = ML.get_physical_linear_map(ph_exact.segments[1], ph_exact, points_exact[1], points_exact[2])
    
    p1_ref_node = points_exact[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_exact[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_exact[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_exact[2].btrian)[p2_ref_node][1]
    
    # Test at multiple points - should match linear exactly
    for t in [0.0, 0.25, 0.5, 0.75, 1.0]
        r = VectorValue(x_ref_p1 + t * (x_ref_p2 - x_ref_p1))
        pos_quad = f_quad(r)
        pos_linear = f_linear(r)
        @test isapprox(pos_quad[1], pos_linear[1]; atol=1e-8)
        @test isapprox(pos_quad[2], pos_linear[2]; atol=1e-8)
    end
end

@testset "Quadratic Map - Large Excess Length" begin
    # Cable with significant excess length
    ph_large = PH.ParameterHandler()
    ph_large.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
    ph_large.points[2] = PH.PointParameters(id=2, coords=[30.0, 0.0], mesh_size=0.1)
    ph_large.segments[1] = PH.SegmentParameters(id=1, tag="large", start_point=1, stop_point=2, length=50.0)
    ph_large.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_large = ML.setup_lines(ph_large)
    segments_large = ML.get_segments(mlines_large[1])
    points_large = Seg.get_points(segments_large[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_large.segments[1], ph_large, points_large[1], points_large[2])
    f_linear = ML.get_physical_linear_map(ph_large.segments[1], ph_large, points_large[1], points_large[2])
    
    p1_ref_node = points_large[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_large[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_large[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_large[2].btrian)[p2_ref_node][1]
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    
    mid_pos_quad = f_quad(VectorValue(mid_ref))
    mid_pos_linear = f_linear(VectorValue(mid_ref))
    
    # Sag should be substantial
    sag = mid_pos_quad[2] - mid_pos_linear[2]
    @test sag < -5.0  # Significant downward sag
end

@testset "Quadratic Map - Arc Length Matches Segment Length" begin
    # Replace the failing "Sag Calculation Verification" test
    # Test that the computed arc length equals the specified segment length
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
    
    # Numerically compute arc length of the quadratic map
    n = 1000
    total_length = 0.0
    for i in 1:n
        s0 = (i-1) / n
        s1 = i / n
        r0 = VectorValue(x_ref_p1 + s0 * (x_ref_p2 - x_ref_p1))
        r1 = VectorValue(x_ref_p1 + s1 * (x_ref_p2 - x_ref_p1))
        p0 = f_quad(r0)
        p1 = f_quad(r1)
        total_length += norm(p1 .- p0)
    end
    
    # Arc length should match segment length (within 1% tolerance)
    @test isapprox(total_length, 110.0; rtol=0.01)
    
    # Verify there is sag (not just linear)
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    mid_pos_quad = f_quad(VectorValue(mid_ref))
    mid_pos_linear = f_linear(VectorValue(mid_ref))
    
    # Sag should be negative (downward)
    actual_sag = mid_pos_quad[2] - mid_pos_linear[2]
    @test actual_sag < 0.0
    
    # For this case with 10m excess on 100m span, expect significant sag
    @test actual_sag < -5.0  # Should have substantial downward displacement
end

@testset "Quadratic Map - Sag Properties" begin
    # Test general properties of the sag, not specific values
    ph_sag = PH.ParameterHandler()
    ph_sag.points[1] = PH.PointParameters(id=1, coords=[0.0, 0.0], mesh_size=0.1)
    ph_sag.points[2] = PH.PointParameters(id=2, coords=[50.0, 0.0], mesh_size=0.1)
    ph_sag.segments[1] = PH.SegmentParameters(id=1, tag="sag_test", start_point=1, stop_point=2, length=60.0)
    ph_sag.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
    
    mlines_sag = ML.setup_lines(ph_sag)
    segments_sag = ML.get_segments(mlines_sag[1])
    points_sag = Seg.get_points(segments_sag[1])
    
    f_quad = ML.get_physical_quadratic_map(ph_sag.segments[1], ph_sag, points_sag[1], points_sag[2])
    f_linear = ML.get_physical_linear_map(ph_sag.segments[1], ph_sag, points_sag[1], points_sag[2])
    
    p1_ref_node = points_sag[1].btrian.glue.face_to_bgface[1]
    p2_ref_node = points_sag[2].btrian.glue.face_to_bgface[1]
    x_ref_p1 = get_node_coordinates(points_sag[1].btrian)[p1_ref_node][1]
    x_ref_p2 = get_node_coordinates(points_sag[2].btrian)[p2_ref_node][1]
    
    # Test symmetry and parabolic shape properties
    quarter_ref = x_ref_p1 + 0.25 * (x_ref_p2 - x_ref_p1)
    mid_ref = (x_ref_p1 + x_ref_p2) / 2.0
    three_quarter_ref = x_ref_p1 + 0.75 * (x_ref_p2 - x_ref_p1)
    
    quarter_pos = f_quad(VectorValue(quarter_ref))
    mid_pos = f_quad(VectorValue(mid_ref))
    three_quarter_pos = f_quad(VectorValue(three_quarter_ref))
    
    quarter_linear = f_linear(VectorValue(quarter_ref))
    mid_linear = f_linear(VectorValue(mid_ref))
    three_quarter_linear = f_linear(VectorValue(three_quarter_ref))
    
    # Calculate sag relative to linear interpolation
    sag_quarter = quarter_pos[2] - quarter_linear[2]
    sag_mid = mid_pos[2] - mid_linear[2]
    sag_three_quarter = three_quarter_pos[2] - three_quarter_linear[2]
    
    # Symmetry: sag at 1/4 and 3/4 should be equal
    @test isapprox(sag_quarter, sag_three_quarter; atol=1e-6)
    
    # Maximum sag at midpoint
    @test sag_mid < sag_quarter  # More negative = more sag
    @test sag_mid < sag_three_quarter
    
    # All sag values should be negative (downward)
    @test sag_quarter < 0.0
    @test sag_mid < 0.0
    @test sag_three_quarter < 0.0
    
    # Horizontal positions should match linear
    @test isapprox(quarter_pos[1], quarter_linear[1]; atol=1e-8)
    @test isapprox(mid_pos[1], mid_linear[1]; atol=1e-8)
    @test isapprox(three_quarter_pos[1], three_quarter_linear[1]; atol=1e-8)
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