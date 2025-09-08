module PointMotion

using Parameters
using Gridap.TensorValues
using Interpolations
using RuntimeGeneratedFunctions

import Mooring.ParameterHandlers as PH
import Mooring.EnvironmentalConditions as EC
import WaveSpec.Constants as WC
import WaveSpec.WaveTimeSeries as WTS

export MotionType, WaveMotionType, CustomMotionType
export get_point_motion

RuntimeGeneratedFunctions.init(@__MODULE__)

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
    - wave_params::WaveParameters: wave Parameters
    - point_motion_function::Function: function that returns the position of a point in time

This struct precomputes the motion of a point given the motion type, input ramp, time and position. For a `WaveMotionType`
the motion is calculated using the wave spectrum and the Airy wave theory. It calculates the position of a point in time
at certain given time instances (`t_range`) and constructs an interpolable function that returns the interpolated value at 
any given time.

This pre-computation is done to avoid evaluating the wave spectrum at each time step.

"""
@with_kw struct WaveMotionType <: MotionType
  wave_params::EC.WaveParameters = EC.WaveParameters()
  point_motion_function::Function = t -> VectorValue(0.0, 0.0)
  function WaveMotionType(wave_params::EC.WaveParameters, input_ramp, t_range::AbstractRange, x::VectorValue{2,Float64})
    sp = EC.get_input_spectrum(wave_params)
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

    point_motion_function(t) = VectorValue(sitp_px(t), sitp_py(t))
    return new(wave_params, point_motion_function)
  end
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
This function returns the motion function of a point given its motion type.
Input:
- `motion_type::MotionType`: Type of motion (WaveMotionType or CustomMotionType)

Output:
- `point_motion_function::Function`: Interpolated position vector. It can be evaluated at any time instance
"""
function get_point_motion_function(motion_type::WaveMotionType) 
    return motion_type.point_motion_function

end

function get_point_motion_function(motion_type::CustomMotionType)
    f(t::Real) = x -> motion_type.f(t,x)
    return f
end

get_point_motion_function(::Nothing) = nothing

function get_point_motion_function(motion_type::MotionType)
    error("`get_point_motion` not implemented for this motion case")
end

"""
get_point_motion

This function returns the motion of a point given its parameters.
Input:
- `p_params::PointParameters`: Parameters of the point
Output:
- `motion::MotionType`: Motion of the point
"""
function get_point_motion(p_id::Int, ph::PH.ParameterHandler)
  p_params = ph.points[p_id]
  motion_params = ph.motions[p_params.motion_tag]
  motion_type = motion_params.type
  if motion_type === nothing
    return nothing
  elseif motion_type == "CustomMotion"
    f = @RuntimeGeneratedFunction(Meta.parse(motion_params.f))
    return CustomMotionType(f)
  elseif motion_type == "WaveMotion"
    error("Not implemented yet")
  end
end

end