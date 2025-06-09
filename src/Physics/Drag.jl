module Drag

using Parameters
using Gridap.TensorValues
import WaveSpec.WaveTimeSeries as WTS

export DragType, NoDrag, Custom, ChainStudless, ChainStudlink
export DragProperties

ρWater = 1025 #Kg/m3

"""
DragType Struct

This abstract type is used to define the type of drag. Possible implemented options are:
- `NoDrag`: No drag
- `Custom`: Custom drag
- `ChainStudless`: Studless chain drag
- `ChainStudlink`: Studlink chain drag
"""
abstract type DragType end

"""
NoDrag Struct

Using this type all drag properties are set to zero.
"""
struct NoDrag <:DragType end

"""
Custom Struct

This type is used to define custom drag properties.
"""
struct Custom <:DragType end

"""
ChainStudless Struct

This type is used to define default drag properties for a studless chain.
"""
struct ChainStudless <: DragType end

"""
ChainStudlink Struct

This type is used to define default drag properties for a studlink chain.
"""
struct ChainStudlink <: DragType end


"""
DragProperties Struct

This struct contains the properties of the drag. The following parameters 
are included, no default values are assumed.
- `dragType::DragType`: Type of drag
- `ρw::Real`: Water density [Kg/m3]
- `nd::Real`: Nominal diameter [m]
- `od::Real`: Outer diameter [m]
- `id::Real`: Inner diameter [m]
- `AStr::Real`: Area (structural) [m2]
- `Cd_n::Real`: Drag coefficient normal
- `Cd_t::Real`: Drag coefficient tangent
- `dd_n::Real`: Drag diameter normal
- `dd_t::Real`: Drag diameter tangent
- `Ca_n::Real`: Added mass coefficient normal
- `Ca_t::Real`: Added mass coefficient tangent
- `Cfd_n::Real`: Friction coefficient normal
- `Cfd_t::Real`: Friction coefficient tangent
"""
@with_kw struct DragProperties
	
	dragType::DragType
	ρw::Real

	nd::Real  # Nominal diameter
	od::Real	# Outer diameter
	id::Real	# Inner diameter

	AStr::Real		# Area (structural)

	Cd_n::Real # Drag coeff normal	
	Cd_t::Real # Drag coeff tangent

	dd_n::Real 	# Drag dia normal
	dd_t::Real 	# Drag dia tangent

	Ca_n::Real  # Added mass normal
	Ca_t::Real  # Added mass tangent

	Cfd_n::Real
	Cfd_t::Real

  function DragProperties( dragType::DragType, ρw, 
    nd, od, id, AStr,
    Cd_n, Cd_t, dd_n, dd_t,
    Ca_n, Ca_t )
  
    Cfd_n = 0.5 * ρw * dd_n * Cd_n / AStr #kg/m4
    Cfd_t = 0.5 * ρw * π * dd_t * Cd_t / AStr #kg/m4
  
    new( dragType, ρw, 
      nd, od, id, AStr,
      Cd_n, Cd_t, dd_n, dd_t,
      Ca_n, Ca_t, Cfd_n, Cfd_t )
  
  end
end

function DragProperties(dragType::DragType) 
  error("Invalid DragType for default constructor.")
end

"""
DragProperties Constructor

This function constructs the drag properties for the `NoDrag` type.

Input:
- `dragType::NoDrag`: Type of drag

Output:
- `DragProperties`: Drag properties
"""
function DragProperties( dragType::NoDrag )

  DragProperties( NoDrag(), 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0 )
end

"""
DragProperties Constructor

This function constructs the drag properties for the `ChainStudless` type.

Input:
- `dragType::ChainStudless`: Type of drag
- `nd::Real`: Nominal diameter [m]
- `AStr::Real`: Area (structural) [m2]
Optional parameters (with default values):
- `Cd_n::Real = 2.4`: Drag coefficient normal
- `Cd_t::Real = 1.15`: Drag coefficient tangent
- `Ca_n::Real = 1.0`: Added mass coefficient
- `Ca_t::Real = 0.5`: Added mass coefficient
- `ρw::Real = 1025`: Water density [Kg/m3]
- `od::Real = 1.80 * nd`: Outer diameter [m]
- `id::Real = 0.0`: Inner diameter [m]
- `dd_n::Real = nd`: Drag diameter normal
- `dd_t::Real = nd / π`: Drag diameter tangent

Output:
- `DragProperties`: Drag properties
"""
function DragProperties( dragType::ChainStudless, nd, AStr;
  Cd_n = 2.4, Cd_t = 1.15, Ca_n = 1.0, Ca_t = 0.5, 
  ρw = ρWater,
  od = 1.80 * nd, id = 0.0, dd_n = nd, dd_t = nd / π )

  DragProperties( dragType, ρw, nd, od, id, AStr,
    Cd_n, Cd_t, dd_n, dd_t,
    Ca_n, Ca_t )
end

"""
DragProperties Constructor

This function constructs the drag properties for the `ChainStudlink` type.

Input:
- `dragType::ChainStudlink`: Type of drag
- `nd::Real`: Nominal diameter [m]
- `AStr::Real`: Area (structural) [m2]
Optional parameters (with default values):
- `Cd_n::Real = 2.6`: Drag coefficient normal
- `Cd_t::Real = 1.4`: Drag coefficient tangent
- `Ca_n::Real = 1.0`: Added mass coefficient
- `Ca_t::Real = 0.5`: Added mass coefficient
- `ρw::Real = 1025`: Water density [Kg/m3]
- `od::Real = 1.89 * nd`: Outer diameter [m]
- `id::Real = 0.0`: Inner diameter [m]
- `dd_n::Real = nd`: Drag diameter normal
- `dd_t::Real = nd / π`: Drag diameter tangent
"""
function DragProperties( dragType::ChainStudlink, nd, AStr;
  Cd_n = 2.6, Cd_t = 1.4, Ca_n = 1.0, Ca_t = 0.5, 
  ρw = ρWater,
  od = 1.89 * nd, id = 0.0, dd_n = nd, dd_t = nd / π )  

  DragProperties( dragType, ρw, nd, od, id, AStr,
    Cd_n, Cd_t, dd_n, dd_t,
    Ca_n, Ca_t )
end

"""
DragProperties Constructor

This function constructs the drag properties for the `Custom` type.

Input:
- `dragType::Custom`: Type of drag
- `nd::Real`: Nominal diameter [m]
- `AStr::Real`: Area (structural) [m2]
Optional parameters (with default values):
- `Cd_n::Real`: Drag coefficient normal
- `Cd_t::Real`: Drag coefficient tangent
- `Ca_n::Real`: Added mass coefficient
- `Ca_t::Real`: Added mass coefficient
- `ρw::Real = 1025`: Water density [Kg/m3]
- `od::Real = nd`: Outer diameter [m]
- `id::Real = 0.0`: Inner diameter [m]
- `dd_n::Real = nd`: Drag diameter normal
- `dd_t::Real = nd / π`: Drag diameter tangent
"""
function DragProperties( dragType::Custom, nd, AStr;
  Cd_n = 1.0, Cd_t = 1.0, Ca_n = 1.0, Ca_t = 1.0, 
  ρw = ρWater,
  od = nd, id = 0.0, dd_n = nd, dd_t = nd / π )
  
  DragProperties( dragType, ρw, nd, od, id, AStr,
    Cd_n, Cd_t, dd_n, dd_t,
    Ca_n, Ca_t )
end

""" 
relative_velocity

This function calculates the relative velocity between the line and the fluid. The default implementation 
is the negative of the line velocity, \$ v_{\\text{rel}} = -v \$.

Input:
- `v::VectorValue`: Line velocity

Output:
- `VectorValue`: Relative velocity
"""
relative_velocity(v::VectorValue) = -v

"""
relative_velocity

This function calculates the relative velocity between the line and the fluid. This implementation
computes the difference between the line (\$ v \$) and the fluid (\$ u \$) velocities, \$ v_{\\text{rel}} = u - v \$.

Input:
- `u::VectorValue`: Fluid velocity
- `v::VectorValue`: Line velocity

Output:
- `VectorValue`: Relative velocity
"""
relative_velocity(u::VectorValue,v::VectorValue) = u - v

"""
relative_velocity

This function calculates the relative velocity between the line and the fluid. This implementation
computes the difference between the line (\$ v \$) and the fluid (\$ u \$) velocities, assuming a time ramp
function for the fluid velocity (\$ \\alpha(t) \$). The relative velocity is given by \$ v_{\\text{rel}} = \\alpha(t)*u - v \$.

Input:
- `t::Real`: Time
- `input_ramp::WTS.TimeRampTypeAll`: Time ramp function for the fluid velocity
- `u::VectorValue`: Fluid velocity
- `v::VectorValue`: Line velocity

Output:
- `VectorValue`: Relative velocity
"""
function relative_velocity(t::Real, input_ramp::WTS.TimeRampTypeAll,
  u::VectorValue,v::VectorValue)
  t_ramp = timeRamp(t, inputRamp)
  return t_ramp*u - v
end


"""
drag_ΓX

This function calculates the drag force on the line. The drag force is given by:

```math
F_{\\text{drag}} = (C_{\\text{fd}_n} \\cdot v_n \\cdot \\|v_n\\| + C_{\\text{fd}_t} \\cdot v_t \\cdot \\|v_t\\|) \\cdot s_\\Lambda
```

Input:
- `dragProp::DragProperties`: Drag properties
- `sΛ::Real`: line stretch
- `t::VectorValue`: Line tangent
- `vr::VectorValue`: Relative velocity

Output:
- `VectorValue`: Drag force
"""
function drag_ΓX(dragProp::DragProperties, sΛ::Real, t::VectorValue, vr::VectorValue)

  t_norm2 = t ⋅ t
  t_dir = t / √(t_norm2)    
  vt = (vr ⋅ t_dir) * t_dir
  vt_norm = √(vt ⋅ vt)
  vn = vr - vt
  vn_norm = √(vn ⋅ vn)    

  @unpack Cfd_n, Cfd_t = dragProp

  return (Cfd_n * vn * vn_norm + Cfd_t * vt * vt_norm) * sΛ 

end    


# function drag_ΓX(dragProp, QTr, T1s, T1m, ∇u, v) 

#   local FΓ, t1s, t1m2, vn, vnm, sΛ, vt, vtm

#   FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)

#   t1s = FΓ ⋅ T1s
#   t1m2 = t1s ⋅ t1s    
#   #t1 = t1s / ((t1s ⋅ t1s).^0.5)

#   sΛ = (t1m2.^0.5) / T1m
  
#   vt = -(v ⋅ t1s) * t1s / t1m2
#   vtm = (vt ⋅ vt).^0.5
#   vn = -v - vt
#   vnm = (vn ⋅ vn).^0.5    

#   return (dragProp.Cfd_n * vn * vnm + 
#     dragProp.Cfd_t * vt * vtm) * sΛ 

# end    


# # Self Drag + Current, No Wave
# function drag_ΓX(t, dragProp, inputRamp, 
#   UCur, waveVel, QTr, T1s, T1m, ∇u, v)

#   local FΓ, t1s, t1m2, vn, vnm, sΛ, vt, vtm, tRamp

#   tRamp = timeRamp(t, inputRamp)

#   FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)

#   t1s = FΓ ⋅ T1s
#   t1m2 = t1s ⋅ t1s    
#   #t1 = t1s / ((t1s ⋅ t1s).^0.5)

#   sΛ = (t1m2.^0.5) / T1m
  
#   vr = UCur*tRamp + waveVel*tRamp - v
#   vt = (vr ⋅ t1s) * t1s / t1m2
#   vtm = (vt ⋅ vt).^0.5
#   vn = vr - vt
#   vnm = (vn ⋅ vn).^0.5    

#   return (dragProp.Cfd_n * vn * vnm + 
#     dragProp.Cfd_t * vt * vtm) * sΛ 

# end    
# ----------------------End----------------------

end