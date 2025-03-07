import Mooring.EnvironmentalConditions as EC
using Parameters
using WaveSpec.Currents
using WaveSpec.Jonswap
using WaveSpec.Constants
using Gridap.TensorValues

waveparams = WaveParameters()
function read_params(params::WaveParameters)
    @unpack Hs, Tp, h0, nω, seed, ωc, enableWaveSpec = waveparams
    return true
end
@test read_params(waveparams)

# Testing get_current_field
h0 = 23 # Water depth
h_vals = [-23.0, -11.0, 0.0] # Depths of the current field definition
u_vals = [0.0, 10.0, 20.0] # Current velocities at the prescribed depths
currentparams = CurrentStat(h0, h_vals, u_vals)
Xh(r) = VectorValue(0.0, r)
@test EC.get_current_field(0.0, Xh, currentparams) == VectorValue(0.0, 0.0)
@test EC.get_current_field(12.0, Xh, currentparams) == VectorValue(10.0, 0.0)
@test EC.get_current_field(23.0, Xh, currentparams) == VectorValue(20.0, 0.0)

# Testing get_wave_velocity
Hs = 3 # Significant wave height
Tp = 12 # Peak wave period
nω = 65 # Number of frequency components
seed = -1 # Seed for random phase
ωc = 1 # Cut-off frequency
enableWaveSpec = true # Enable wave spectrum
ω, S, A = jonswap(Hs, Tp, plotflag=false, nω = nω, ωc = ωc)
k = dispersionRelAng.(h0, ω; msg=false)
α = randomPhase(ω, seed = seed)
sp = SpecStruct( h0, ω, S, A, k, α; Hs = Hs, Tp = Tp )
@test typeof(EC.get_wave_velocity(0.0, sp, VectorValue(0.0,0.0))) == VectorValue{2,Float64}