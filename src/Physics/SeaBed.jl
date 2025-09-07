module SeaBed

using Parameters
using Gridap.TensorValues
using Mooring.ParameterHandlers: SeaBedParameters

"""
set_still_weight

This function sets the still weight of the sea bed. It assumes the sea bed object
is already created and modifies the still weight parameter.

Input: 
- `params::SeaBedParameters`: Sea bed parameters
- `still_weight::Real`: Still weight [N]

Output:
- `SeaBedParameters`: Sea bed parameters with the still weight modified
"""
function set_still_weight(params::SeaBedParameters, still_weight::Real)
    return reconstruct(params, still_weight=still_weight)
end

"""
ramp_tanh

This function returns the tanh ramp function for a given excursion into the sea bed.
It assumes the excursion is positive when there is penetration to the sea bed and negative when 
the line is lifted from the soil. The function is defined as:

```math
x_{\\text{new}} = \\max(0, 2 \\tanh( \\text{tanh_ramp} \\cdot x ) ).
```

Input:
- `params::SeaBedParameters`: Sea bed parameters
- `exc::Real`: excursion value at a given time

Output:
- `Real`: Tanh ramp function value
"""
function ramp_tanh(params::SeaBedParameters, excursion::Real)
    return max(0.0, 2*tanh( params.tanh_ramp * excursion ) )
end

"""
ramp_linear

This function returns the linear ramp function for a given excursion into the sea bed. It assumes
the excursion is positive when there is penetration into the soil and computes the value based on 
the `penetration_depth_ramp` parameter: 

```math
x_{\\text{new}} = \\frac{x}{\\text{penetration_depth_ramp}}.
```

Input:
- `params::SeaBedParameters`: Sea bed parameters
- `exc::Real`: excursion value at a given time

Output:
- `Real`: Linear ramp function value
"""
function ramp_linear(params::SeaBedParameters, excursion::Real)
    if excursion > 0
        return excursion / params.penetration_depth_ramp
    end
    return 0.0
end


"""
sea_bed_force

This function computes the force exerted by the sea bed on the mooring line. It considers the
excursion of the line into the sea bed, the velocity of the line, and the properties of the sea bed.
The function is defined as:

```math
\\alpha_1 * w + \\alpha_2 * k * s_\\Lambda * ( x_z - \\beta_1 * v_z - \\beta_2 * v_z * \\abs(v_z) )
```

Input:
- `params::SeaBedParameters`: Sea bed parameters
- `X::VectorValue`: Position of the line
- `QTr::TensorValue`: Transformation matrix \$ Q^T \$
- `T1s::TensorValue`: Stress tensor
- `T1m::Real`: Maximum stress
- `u::VectorValue`: Displacement of the line
- `∇u::TensorValue`: Gradient of the displacement
- `v::VectorValue`: Velocity of the line

Output:
- `VectorValue`: Force exerted by the sea bed at a given point in the line
"""
function sea_bed_force(params::SeaBedParameters, X::VectorValue, 
    QTr::TensorValue, T1s::VectorValue, T1m::Real, 
    u::VectorValue, ∇u::TensorValue, v::VectorValue)
  
    # Define local variables
    local excursion, ramp_factor
    local FΓ, t1s, t1m2, sΛ        
  
    excursion = VectorValue(0.0,-1.0) ⋅ (X + u) # assumes flat bed in the x-y plane
    ramp_factor = ramp_tanh(params, excursion)
    half_ramp_factor = min(1.0, ramp_factor)
  
    vz = VectorValue(0.0, 1.0) ⋅ v # velocity in the z direction
  
    FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
    t1s = FΓ ⋅ T1s
    t1m2 = t1s ⋅ t1s    
  
    sΛ = (t1m2.^0.5) / T1m # Stretch along the line
    
    # params.cnstz = params.kn * params.od / params.A
    @unpack still_weight, cnstz, linear_damping_factor, quadratic_damping_factor = params

    return half_ramp_factor * still_weight  + 
      ramp_factor * cnstz * sΛ * ( 
        excursion +
        -linear_damping_factor * vz +
        -quadratic_damping_factor * vz * abs(vz) 
      )
  
  end

end