import Mooring.PointMotion as PM
import WaveSpec.WaveTimeSeries as WTS

# Testing PointMotionType
custom_type = PM.CustomMotionType()
@test PM.get_point_motion(custom_type, nothing, 1.0, VectorValue(1.0,1.0)) == VectorValue(0.0,0.0)

# Testing WaveMotionType
wave_type = PM.WaveMotionType()
input_ramp = WTS.TimeRampType(0.0,1.0)
wave_motion = PM.get_point_motion(wave_type, input_ramp, 0.:0.5:1.0, VectorValue(1.0,1.0))
@test wave_motion(1.0) == VectorValue(0.0,0.0)