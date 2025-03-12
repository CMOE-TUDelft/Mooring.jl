module PointMotion

using Parameters
import WaveSpec.Constants as WC
import WaveSpec.WaveTimeSeries as WTS

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
@with_kw struct IrregularWaveMotion <: MotionType
  wave_params::WaveParameters = WaveParameters()
end

"""
get_point_motion

This function returns the motion of a point given the motion type, input ramp, time and position.

Input:
- `motion_type::MotionType`: Type of motion
- `input_ramp::Function`: Input ramp function
- `t::Real`: Time
- `x::VectorValue`: Position

Output:
- `sitp_η::Function`: Interpolated function for the wave amplitude including ramping. It can be evaluated at any time instance
- `sitp_px::Function`: Interpolated function for the horizontal displacement. It can be evaluated at any time instance
- `sitp_py::Function`: Interpolated function for the vertical displacement. It can be evaluated at any time instance
"""
function get_point_motion(motion_type::RegularWaveMotion, input_ramp, t::Real, x::VectorValue{2,Float64})
    sp = get_input_spectrum(motion_type.wave_params)
    η, px, py = WTS.waveAiry1D_pPos(sp, t, x[1], x[2])
    t_ramp = timeRamp(t,input_ramp0)

    itp = Interpolations.interpolate(
      η.*t_ramp, 
      BSpline(Cubic(Line(OnGrid()))) )
    sitp_η = scale(itp, t)
  
    itp = Interpolations.interpolate(
      px.*t_ramp, 
      BSpline(Cubic(Line(OnGrid()))) )
    sitp_px = scale(itp, t)
  
    itp = Interpolations.interpolate(
      py.*t_ramp, 
      BSpline(Cubic(Line(OnGrid()))) )
    sitp_py = scale(itp, t)
  
    return sitp_η, sitp_px, sitp_py

end

function get_point_motion(motion_type::MotionType, input_ramp, t::Real, x::VectorValue{3,Float64})
    error("`get_point_motion` not implemented for the 3D case")
end

end