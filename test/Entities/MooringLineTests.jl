import Mooring.ParameterHandler as PH
import Mooring.MooringLines as ML
using Gridap.MultiField

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
f = ML.get_physical_map(ph.segments[1], ph, start_point, stop_point)
mid = f(VectorValue(0.75)) # halfway along segment
@test isapprox(mid[1], 0.5; atol=1e-8)
@test isapprox(mid[2], 0.0; atol=1e-8)

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