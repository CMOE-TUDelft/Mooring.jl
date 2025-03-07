module SeaBed

using Parameters
using Gridap.TensorValues

export SeaBedParams

""" 
SeaBedParams Struct

This struct contains the properties of the seabed.
The following parameters are included, with default values:
- `kn::Real = 30e3`: Normal stiffness [N/m2]
- `linear_damping_ratio::Real = 0.05`: Linear damping ratio [s]
- `quadratic_damping_ratio::Real = 0.0`: Quadratic damping ratio [s^2/m]
- `od::Real = 0.1`: Outer diameter of the line [m]
- `A::Real = 0.008`: Area of the line [m^2]
- `tanh_ramp::Real = 1e2`: Tanh ramp function parameter 
- `penetration_depth_ramp::Real = 1e-3`: Penetration depth ramp function parameter [m]
- `still_weight::Real = 0.0`: Still weight [N]
- `cnstz::Real = 0.0`: Constant spring stiffness of the sea bed [N/m]

Relevant references:
- Quadratic law impact damping: https://doi.org/10.1080/0020739X.2021.1954253
- Critical damping of Moordyn: https://moordyn.readthedocs.io/en/latest/troubleshooting.html#model-stability-and-segment-damping
"""
@with_kw struct SeaBedParams
    kn::Real = 30e3
    linear_damping_ratio::Real = 0.05
    quadratic_damping_ratio::Real = 0.0  
    od::Real = 0.1
    A::Real = 0.008 
    tanh_ramp::Real = 1.0e2
    penetration_depth_ramp::Real = 1.0e-3
    still_weight::Real = 0.0
    cnstz::Real = kn * od / A
end

"""
set_still_weight

This function sets the still weight of the sea bed. It assumes the sea bed object
is already created and modifies the still weight parameter.

Input: 
- `params::SeaBedParams`: Sea bed parameters
- `still_weight::Real`: Still weight [N]

Output:
- `SeaBedParams`: Sea bed parameters with the still weight modified
"""
function set_still_weight(params::SeaBedParams, still_weight::Real)
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
- `params::SeaBedParams`: Sea bed parameters
- `exc::Real`: excursion value at a given time

Output:
- `Real`: Tanh ramp function value
"""
function ramp_tanh(params::SeaBedParams, excursion::Real)
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
- `params::SeaBedParams`: Sea bed parameters
- `exc::Real`: excursion value at a given time

Output:
- `Real`: Linear ramp function value
"""
function ramp_linear(params::SeaBedParams, excursion::Real)
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

[TO DO]
```math
````

Input:
- `params::SeaBedParams`: Sea bed parameters
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
function sea_bed_force(params::SeaBedParams, X::VectorValue, 
    QTr::TensorValue, T1s::VectorValue, T1m::Real, 
    u::VectorValue, ∇u::TensorValue, v::VectorValue)
  
    # Define local variables
    local excursion, lSpng
    local FΓ, t1s, t1m2, sΛ        
  
    excursion = VectorValue(0.0,-1.0) ⋅ (X + u) # assumes flat bed in the x-y plane
    lSpng = ramp_tanh(params, excursion)
    lstill_weight = min(1.0, lSpng)
  
    vz = VectorValue(0.0, 1.0) ⋅ v # velocity in the z direction
  
    FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
    t1s = FΓ ⋅ T1s
    t1m2 = t1s ⋅ t1s    
  
    sΛ = (t1m2.^0.5) / T1m
    
    # params.cnstz = params.kn * params.od / params.A
  
    return lstill_weight * params.still_weight  + 
      lSpng * params.cnstz * sΛ * ( 
        excursion +
        -params.linear_damping_ratio * vz +
        -params.quadratic_damping_ratio * vz * abs(vz) 
      )
  
  end

end