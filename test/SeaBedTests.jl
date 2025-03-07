import Mooring.SeaBed as SB
using Parameters

sea_bed_params = SB.SeaBedParams()
function read_params(params::SB.SeaBedParams)
    @unpack kn, linDampRatio, quadDampRatio, od, A, tanh_ramp, penDepth_ramp, stillWei, cnstz = sea_bed_params
    return true
end
@test read_params(sea_bed_params)
