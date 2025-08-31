module EnvironmentalConditions

using Parameters
using Gridap.TensorValues
using WaveSpec.Currents
using WaveSpec.WaveTimeSeries
using WaveSpec.Jonswap
using WaveSpec.Constants

"""
get_input_spectrum

This function returns the wave spectrum parameters for a given set of wave parameters. The wave spectrum
is calculated using the package [WaveSpec](https://github.com/shagun751/WaveSpec.jl).
"""
function get_input_spectrum(params::WaveParameters)
  @unpack Hs, Tp, h0, nω, seed, ωc, enableWaveSpec = params

  if(enableWaveSpec)
    if(ωc < 0)
      ω, S, A = jonswap(Hs, Tp,
        plotflag=false, nω = nω)
    else
      ω, S, A = jonswap(Hs, Tp,
        plotflag=false, nω = nω, ωc = ωc)
    end
  
  else
    ω, S, A = jonswap(0.0, 10.0,
        plotflag=false, nω = 64)
  end

  k = dispersionRelAng.(h0, ω; msg=false)
  α = randomPhase(ω, seed = seed)

  sp = SpecStruct( h0, ω, S, A, k, α; Hs = Hs, Tp = Tp )
  return sp
end

"""
get_current_field

This function returns the current field at a given point
of the undeformed configuration \$ r \\rightarrow X_h(r) \$.

Input:
- `r::Real`: Point in the undeformed configuration
- `Xh::Function`: Function that maps the undeformed configuration to the deformed configuration (VectorValue)
- `curObj::CurrentStat`: Current object

Output:
- `VectorValue`: Current field at the point `r`
"""
function get_current_field(r::Real, Xh, curObj::CurrentStat)
  
    X_qp = Xh(r)
    pz = X_qp ⋅ VectorValue(0.0,1.0) - curObj.h0
    pz = min(pz, 0.0)
    pz = max(pz, -curObj.h0)
  
    return VectorValue( curObj.itp( pz ), 0.0 )   
  
end

"""
get_wave_velocity

This function returns the wave velocity at a given point
in the space `x`. It uses the Airy wave theory, defined in the [WaveSpec](https://github.com/shagun751/WaveSpec.jl)
package.

Input:
- `t::Real`: Time
- `sp::SpecStruct`: Wave spectrum parameters
- `x::VectorValue`: Point in the space

Output:
- `VectorValue`: Wave velocity at the point `x`
"""
function get_wave_velocity(t, sp, x)
    w_u, w_w = waveAiry1D_vel(sp, t, x[1], x[2]-sp.h0 )
    return VectorValue(w_u, w_w)
end



end