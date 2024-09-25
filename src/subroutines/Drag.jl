module Drag

using Revise
using Parameters
using Gridap
using WaveSpec.Constants

ρWater = 1025 #Kg/m3

"""
Custom Structs
=============

"""
# ---------------------Start---------------------
abstract type DragType end

struct NoDrag <:DragType end
struct Custom <:DragType end
struct ChainStudless <: DragType end
struct ChainStudlink <: DragType end


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


function DragProperties( dragType::NoDrag )

  DragProperties( NoDrag(), 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0 )
end


function DragProperties( dragType::ChainStudless, nd, AStr;
  Cd_n = 2.4, Cd_t = 1.15, Ca_n = 1.0, Ca_t = 0.5, 
  ρw = ρWater,
  od = 1.80 * nd, id = 0.0, dd_n = nd, dd_t = nd / π )

  DragProperties( dragType, ρw, nd, od, id, AStr,
    Cd_n, Cd_t, dd_n, dd_t,
    Ca_n, Ca_t )
end


function DragProperties( dragType::ChainStudlink, nd, AStr;
  Cd_n = 2.6, Cd_t = 1.4, Ca_n = 1.0, Ca_t = 0.5, 
  ρw = ρWater,
  od = 1.89 * nd, id = 0.0, dd_n = nd, dd_t = nd / π )  

  DragProperties( dragType, ρw, nd, od, id, AStr,
    Cd_n, Cd_t, dd_n, dd_t,
    Ca_n, Ca_t )
end


function DragProperties( dragType::Custom, nd, AStr;
  Cd_n = 1.0, Cd_t = 1.0, Ca_n = 1.0, Ca_t = 1.0, 
  ρw = ρWater,
  od = nd, id = 0.0, dd_n = nd, dd_t = nd / π )
  
  DragProperties( dragType, ρw, nd, od, id, AStr,
    Cd_n, Cd_t, dd_n, dd_t,
    Ca_n, Ca_t )
end
# ----------------------End----------------------



"""
Functions
=============

"""
# ---------------------Start---------------------
function drag_ΓX(dragProp, QTr, T1s, T1m, ∇u, v) #No wave and current

  local FΓ, t1s, t1m2, vn, vnm, sΛ, vt, vtm

  FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)

  t1s = FΓ ⋅ T1s
  t1m2 = t1s ⋅ t1s    
  #t1 = t1s / ((t1s ⋅ t1s).^0.5)

  sΛ = (t1m2.^0.5) / T1m
  
  vt = -(v ⋅ t1s) * t1s / t1m2
  vtm = (vt ⋅ vt).^0.5
  vn = -v - vt
  vnm = (vn ⋅ vn).^0.5    

  return (dragProp.Cfd_n * vn * vnm + 
    dragProp.Cfd_t * vt * vtm) * sΛ 

end    
# ----------------------End----------------------

end