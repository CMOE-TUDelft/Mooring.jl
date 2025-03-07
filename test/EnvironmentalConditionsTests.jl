import Mooring.EnvironmentalConditions as EC
using Parameters
using WaveSpec.Currents
using Gridap.TensorValues

waveparams = WaveParameters()
function read_params(params::WaveParameters)
    @unpack Hs, Tp, h0, nω, seed, ωc, enableWaveSpec = waveparams
    return true
end
@test read_params(waveparams)

# Testing getCurrentField
h0 = 23 # Water depth
h_vals = [-23.0, -11.0, 0.0] # Depths of the current field definition
u_vals = [0.0, 10.0, 20.0] # Current velocities at the prescribed depths
currentparams = CurrentStat(h0, h_vals, u_vals)
Xh(r) = VectorValue(0.0, r)
@test EC.getCurrentField(0.0, Xh, currentparams) == VectorValue(0.0, 0.0)
@test EC.getCurrentField(12.0, Xh, currentparams) == VectorValue(10.0, 0.0)
@test EC.getCurrentField(23.0, Xh, currentparams) == VectorValue(20.0, 0.0)
