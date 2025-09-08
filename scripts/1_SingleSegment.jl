module SingleSegment
import Mooring.ParameterHandlers as PH
import Mooring.MooringLines as ML

# Define ParameterHandlers
ph = PH.ParameterHandler()
ph.points[1] = PH.PointParameters(id=1, coords=[0.0,0.0], motion_tag="fixed", mesh_size=1.0)
ph.points[2] = PH.PointParameters(id=2, coords=[10.0,5.0], motion_tag="custom_motion", mesh_size=1.0)
ph.segments[1] = PH.SegmentParameters(id=1, start_point=1, stop_point=2, material_tag="steel", area=0.01,
                                      density=7850.0, length=10.0, drag_tag="no_drag", seabed_tag="default_seabed")
ph.lines[1] = PH.LineParameters(id=1, points=[1,2], segments=[1])
ph.materials["steel"] = PH.MaterialParameters(tag="steel", E=2.0e11, μ=8.1e10)
ph.motions["custom_motion"] = PH.MotionParameters(tag="custom_motion", type="CustomMotion", wave_tag="(t,(x,y))->VectorValue(sin(2π * t), 0.0)")
ph.motions["fixed"] = PH.MotionParameters(tag="fixed")

# Setup lines
mlines = ML.setup_lines(ph)
end # module