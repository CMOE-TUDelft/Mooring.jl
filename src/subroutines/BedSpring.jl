module BedSpring

using Revise
using Parameters
using Gridap



"""
Custom Structs
=============

"""
# ---------------------Start---------------------
@with_kw struct Bed
  kn::Real 
  # Based on Marco's suggestion of 30kN/m2/m in the OrcaFlex manual  
  linear_damping_ratio::Real  #sec
  quadratic_damping_ratio::Real  #sec2/m 
  # qudratic law impact damping
  # https://doi.org/10.1080/0020739X.2021.1954253
  # Critical damping of Moordyn
  # https://moordyn.readthedocs.io/en/latest/troubleshooting.html#model-stability-and-segment-damping

  od::Real
  A::Real      
  tanh_ramp::Real
  penetration_depth_ramp::Real #Penetration depth ramp #m

  still_weight::Real    
  cnstz::Real

  function Bed( od::Real, A::Real;
    kn = 30e3, 
    linear_damping_ratio = 0.05,
    quadratic_damping_ratio = 0.0,
    tanh_ramp = 1e2,
    penetration_depth_ramp = 1e-3,    
    still_weight::Real=0.0 )    

    cnstz = kn * od / A
    
    new(kn, linear_damping_ratio, quadratic_damping_ratio,
      od, A, 
      tanh_ramp, penetration_depth_ramp, 
      still_weight, cnstz)
  end
end
# ----------------------End----------------------



"""
Functions
=============

"""
# ---------------------Start---------------------
function set_still_weight(bedObj::Bed, still_weight)
  @unpack od, A, kn, tanh_ramp, 
    linear_damping_ratio, quadratic_damping_ratio, penetration_depth_ramp = bedObj  

  Bed(od, A; 
    kn, linear_damping_ratio, quadratic_damping_ratio,
    tanh_ramp, penetration_depth_ramp, 
    still_weight)
end


function ramp_tanh(bedObj::Bed, exc)
  # return 0.5 + 0.5*( tanh( bedObj.tanh_ramp * exc ) )
  return max(0.0, 2*tanh( bedObj.tanh_ramp * exc ) )
end


function ramp_linear(bedObj::Bed, exc)
  if(exc >0) 
    return exc / bedObj.penetration_depth_ramp
  end
  return 0.0
end


# Bed1
# function forceFnc(bedObj::Bed, X, QTr, T1s, T1m, u, ∇u, v)
  
#   local exc, lspng
#   local FΓ, t1s, t1m2, sΛ        

#   exc = VectorValue(0.0,-1.0) ⋅ (X + u)
#   lspng = 0.5 + 0.5*( tanh( bedObj.tanh_ramp * exc ) )

#   vz = VectorValue(0.0, 1.0) ⋅ v

#   FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
#   t1s = FΓ ⋅ T1s
#   t1m2 = t1s ⋅ t1s    

#   sΛ = (t1m2.^0.5) / T1m
  
#   # bedObj.cnstz = bedObj.kn * bedObj.od / bedObj.A

#   return lspng * bedObj.still_weight  + 
#     # lspng * bedObj.kn * bedObj.od / bedObj.A * exc * sΛ -
#     # lspng * bedObj.linear_damping_ratio* bedObj.kn * bedObj.od / bedObj.A * vz * sΛ
#     lspng * bedObj.cnstz * sΛ * ( exc - bedObj.linear_damping_ratio * vz )

# end


# # Bed2
# function forceFnc(bedObj::Bed, X, QTr, T1s, T1m, u, ∇u, v)
  
#   local exc, lSpng, lstill_weight
#   local FΓ, t1s, t1m2, sΛ        

#   exc = VectorValue(0.0,-1.0) ⋅ (X + u)
#   # lSpng = 0.5 + 0.5*( tanh( bedObj.tanh_ramp * exc ) )

#   lSpng = 0.0
#   if(exc >0) 
#     lSpng = exc / bedObj.penetration_depth_ramp
#   end
#   lstill_weight = min(1.0, lSpng)
  
#   vz = VectorValue(0.0, 1.0) ⋅ v  

#   FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
#   t1s = FΓ ⋅ T1s
#   t1m2 = t1s ⋅ t1s    

#   sΛ = (t1m2.^0.5) / T1m
  
#   return lstill_weight * bedObj.still_weight  + 
#     lSpng * bedObj.cnstz * sΛ * (exc - bedObj.linear_damping_ratio * vz) 

# end


# Bed1b
function forceFnc(bedObj::Bed, X, QTr, T1s, T1m, u, ∇u, v)
  
  local exc, lSpng
  local FΓ, t1s, t1m2, sΛ        

  exc = VectorValue(0.0,-1.0) ⋅ (X + u)
  lSpng = ramp_tanh(bedObj, exc)
  lstill_weight = min(1.0, lSpng)

  vz = VectorValue(0.0, 1.0) ⋅ v

  FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
  t1s = FΓ ⋅ T1s
  t1m2 = t1s ⋅ t1s    

  sΛ = (t1m2.^0.5) / T1m
  
  # bedObj.cnstz = bedObj.kn * bedObj.od / bedObj.A

  return lstill_weight * bedObj.still_weight  + 
    lSpng * bedObj.cnstz * sΛ * ( 
      exc +
      -bedObj.linear_damping_ratio * vz +
      -bedObj.quadratic_damping_ratio * vz * abs(vz) 
    )

end
# ----------------------End----------------------

end 