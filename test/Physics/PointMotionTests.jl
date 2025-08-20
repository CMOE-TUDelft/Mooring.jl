import Mooring.PointMotion as PM
import Mooring.EnvironmentalConditions as EC
import WaveSpec.WaveTimeSeries as WTS

# Testing PointMotionType
custom_type = PM.CustomMotionType()
@test PM.get_point_motion_function(custom_type)(1.0)(VectorValue(1.0,1.0)) == VectorValue(0.0,0.0)

# Testing WaveMotionType
wave_params = EC.WaveParameters()
input_ramp = WTS.TimeRampType(0.0,1.0)
wave_type = PM.WaveMotionType(wave_params, input_ramp, 0.:0.5:1.0, VectorValue(1.0,1.0))
wave_motion = PM.get_point_motion_function(wave_type)
@test wave_motion(1.0) == VectorValue(0.0,0.0)