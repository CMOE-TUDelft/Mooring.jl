module SeaBed

using Parameters

export SeaBedParams

""" 
SeaBedParams Struct

This struct contains the properties of the seabed.
The following parameters are included, with default values:
- `kn::Real = 30e3`: Normal stiffness [N/m2]
- `linDampRatio::Real = 0.05`: Linear damping ratio [s]
- `quadDampRatio::Real = 0.0`: Quadratic damping ratio [s^2/m]
- `od::Real` = 0.1: Outer diameter of the line [m]
- `A::Real` = 0.008: Area of the line [m^2]
- `tanh_ramp::Real = 1e2`: Tanh ramp function parameter 
- `penDepth_ramp::Real = 1e-3`: Penetration depth ramp function parameter [m]
- `stillWei::Real = 0.0`: Still weight [N]
- `cnstz::Real` = 0.0 : Constant spring stiffness of the sea bed [N/m]

Relevant references:
- Quadratic law impact damping: https://doi.org/10.1080/0020739X.2021.1954253
- Critical damping of Moordyn: https://moordyn.readthedocs.io/en/latest/troubleshooting.html#model-stability-and-segment-damping
"""
@with_kw struct SeaBedParams
    kn::Real = 30e3
    linDampRatio::Real = 0.05
    quadDampRatio::Real = 0.0  
    od::Real = 0.1
    A::Real = 0.008 
    tanh_ramp::Real = 1.0e2
    penDepth_ramp::Real = 1.0e-3
    stillWei::Real = 0.0
    cnstz::Real = 0.0
end

"""
SeaBedParams Constructor

This function creates a new instance of the SeaBedParams struct when only given
the outer diameter `od` and the area `A` of the line.
"""
function SeaBedParams( od::Real, A::Real)    
    kn = 30e3
    cnstz = kn * od / A
    return SeaBedParams(od=od, A=A, cnstz=cnstz)
end

end