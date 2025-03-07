module EnvironmentalConditions

using Parameters
using Gridap.TensorValues
using WaveSpec.Currents

export WaveParameters

"""
WaveParameters

This struct contains the parameters for the wave conditions.
The following parameters are included, with default values:
- `Hs::Real = 0.0`: Significant wave height
- `Tp::Real = 0.0`: Peak wave period
- `h0::Real = 0.0`: Water depth
- `nω::Int = 64`: Number of frequency components
- `seed::Int = 0`: Seed for random phase
- `ωc::Real = -1.0`: Cut-off frequency
- `enableWaveSpec::Bool = false`: Enable wave spectrum
"""
@with_kw struct WaveParameters

  Hs::Real = 0.0
  Tp::Real = 0.0
  h0::Real = 0.0
  nω::Int = 64
  seed::Int = 0
  ωc::Real = -1.0
  enableWaveSpec::Bool = false

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

function getWaveVel(t, sp, x)
    
    w_u, w_w = waveAiry1D_vel(sp, t, 
      x[1], x[2]-sp.h0 )
  
    return VectorValue(w_u, w_w)
  end



end