module EnvironmentalConditions

using Parameters
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




end