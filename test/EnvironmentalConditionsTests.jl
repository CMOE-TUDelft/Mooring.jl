using Mooring.EnvironmentalConditions
using Parameters

waveparams = WaveParameters()
function read_params(params::WaveParameters)
    @unpack Hs, Tp, h0, nω, seed, ωc, enableWaveSpec = waveparams
    return true
end
@test read_params(waveparams)