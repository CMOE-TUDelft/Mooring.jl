import Mooring.SeaBed as SB
using Parameters

# Testing SeaBedParameters
sea_bed_params = SB.SeaBedParameters()
function read_params(params::SB.SeaBedParameters)
    @unpack kn, linear_damping_factor, quadratic_damping_factor, od, A, tanh_ramp, penetration_depth_ramp, still_weight, cnstz = sea_bed_params
    return true
end
@test read_params(sea_bed_params)
@test read_params(SB.SeaBedParameters(od=1.0, A=1.0))

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

# Testing sea_bed_force
X = VectorValue(0.0, 0.0)
QTr = TensorValue(1.0, 0.0, 0.0, 1.0)
T1s = VectorValue(1.0, 0.0)
T1m = 2.0
u = VectorValue(0.0, 0.0)
∇u = TensorValue(1.0, 0.0, 0.0, 1.0)
v = VectorValue(0.0, 0.0)
@test SB.sea_bed_force(sea_bed_params, X, QTr, T1s, T1m, u, ∇u, v) == 0.0
u = VectorValue(0.0, -1.0)
@test SB.sea_bed_force(sea_bed_params, X, QTr, T1s, T1m, u, ∇u, v) == 2.0 * sea_bed_params.cnstz
