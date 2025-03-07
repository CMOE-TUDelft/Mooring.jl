module SeaBed

using Parameters

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

This function returns the tanh ramp function for a given excitation.

Input:
- `bedObj::SeaBedParams`: Sea bed parameters
- `exc::Real`: Excitation value at a given time

Output:
- `Real`: Tanh ramp function value
"""
function ramp_tanh(params::SeaBedParams, excitation::Real)
    return max(0.0, 2*tanh( params.tanh_ramp * excitation ) )
end

end