import Mooring.SeaBed as SB
using Parameters

# Testing SeaBedParams
sea_bed_params = SB.SeaBedParams()
function read_params(params::SB.SeaBedParams)
    @unpack kn, linear_damping_ratio, quadratic_damping_ratio, od, A, tanh_ramp, penetration_depth_ramp, still_weight, cnstz = sea_bed_params
    return true
end
@test read_params(sea_bed_params)
@test read_params(SB.SeaBedParams(od=1.0, A=1.0))

# Testing set_still_weight
@test sea_bed_params.still_weight == 0.0
modified_params = SB.set_still_weight(sea_bed_params, 1.0)
@test modified_params.still_weight == 1.0

# Testing ramp_tanh
@test SB.ramp_tanh(sea_bed_params, 1.0) == 2.0
@test SB.ramp_tanh(sea_bed_params, -1.0) == 0.0

# Testing ramp_linear
@test SB.ramp_linear(sea_bed_params, 1.0) == 1000.0
@test SB.ramp_linear(sea_bed_params, -1.0) == 0.0
