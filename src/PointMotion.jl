module PointMotion

using Parameters
using Gridap.TensorValues
using Interpolations
import Mooring.EnvironmentalConditions as EC
import WaveSpec.Constants as WC
import WaveSpec.WaveTimeSeries as WTS

export MotionType, WaveMotionType, CustomMotionType
export get_point_motion

"""
MotionType Struct

This abstract type is used to define the type of motion. Possible implemented options are:
- `WaveMotion`: Motion given by a waves (regular or irregular)
- `CustomMotion`: Custom function
"""
abstract type MotionType end

"""
WaveMotion Struct 

This struct is used to define the properties of the wave motion. This includes:
    - wave_params::WaveParameters: wave parameters
"""
@with_kw struct WaveMotionType <: MotionType
  wave_params::EC.WaveParameters = EC.WaveParameters()
end

"""
CustomMotion Struct 

This struct is used to define the properties of the custom motion. This requires a function of space and time.
    - f::Function: custom functions
"""
@with_kw struct CustomMotionType <: MotionType
  f::Function = (t,x) -> VectorValue(0.0,0.0)
end

"""
get_point_motion

This function returns the motion of a point given the motion type, input ramp, time and position. For a `WaveMotionType`
the motion is calculated using the wave spectrum and the Airy wave theory. It calculates the position of a point in time
at certain given time instances (`t_range`) and constructs an interpolable function that returns the interpolated value at 
any given time.

This pre-computation is done to avoid evaluating the wave spectrum at each time step.

Input:
- `motion_type::MotionType`: Type of motion
- `input_ramp::Function`: Input ramp function
- `t_range::AbstractRange`: Vector of time instances from which the position in time is calculated
- `x::VectorValue`: Position

Output:
- `point_position::Function`: Interpolated position vector. It can be evaluated at any time instance
"""
function get_point_motion(motion_type::WaveMotionType, input_ramp, t_range::AbstractRange, x::VectorValue{2,Float64})
    sp = EC.get_input_spectrum(motion_type.wave_params)
    η, px, py = WTS.waveAiry1D_pPos(sp, t_range, x[1], x[2])
    t_ramp = WTS.timeRamp(t_range,input_ramp)
  
    itp = Interpolations.interpolate(
      px.*t_ramp, 
      BSpline(Cubic(Line(OnGrid()))) )
    sitp_px = scale(itp, t_range)
  
    itp = Interpolations.interpolate(
      py.*t_ramp, 
      BSpline(Cubic(Line(OnGrid()))) )
    sitp_py = scale(itp, t_range)

    point_position(t) = VectorValue(sitp_px(t), sitp_py(t))
  
    return point_position

end

function get_point_motion(motion_type::CustomMotionType, input_ramp, t::Real, x::VectorValue{2,Float64})
    return motion_type.f(t,x)
end

function get_point_motion(motion_type::MotionType, input_ramp, t::Real, x::VectorValue{3,Float64})
    error("`get_point_motion` not implemented for the 3D case")
end

end