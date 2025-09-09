import Mooring.ParameterHandler as PH
import Mooring.MooringLines as ML

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

f = ML.get_physical_map(ph.segments[1], ph)
mid = f(VectorValue(5.0)) # halfway along segment
@test isapprox(mid[1], 5.0; atol=1e-8)
@test isapprox(mid[2], 0.0; atol=1e-8)