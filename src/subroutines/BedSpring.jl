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
  linDampRatio::Real  #sec
  quadDampRatio::Real  #sec2/m 
  # qudratic law impact damping
  # https://doi.org/10.1080/0020739X.2021.1954253

  od::Real
  A::Real      
  tanh_ramp::Real
  penDepth_ramp::Real #Penetration depth ramp #m

  stillWei::Real    
  cnstz::Real

  function Bed( od::Real, A::Real;
    kn = 30e3, 
    linDampRatio = 0.05,
    quadDampRatio = 0.0,
    tanh_ramp = 1e2,
    penDepth_ramp = 1e-3,    
    stillWei::Real=0.0 )    

    cnstz = kn * od / A
    
    new(kn, linDampRatio, quadDampRatio,
      od, A, 
      tanh_ramp, penDepth_ramp, 
      stillWei, cnstz)
  end
end
# ----------------------End----------------------



"""
Functions
=============

"""
# ---------------------Start---------------------
function setStillWei!(bedObj::Bed, stillWei)
  @unpack od, A, kn, tanh_ramp, 
    linDampRatio, quadDampRatio, penDepth_ramp = bedObj  

  Bed(od, A; 
    kn, linDampRatio, quadDampRatio,
    tanh_ramp, penDepth_ramp, 
    stillWei)
end


function rampTanh(bedObj::Bed, exc)
  # return 0.5 + 0.5*( tanh( bedObj.tanh_ramp * exc ) )
  return max(0.0, 2*tanh( bedObj.tanh_ramp * exc ) )
end


function rampLin(bedObj::Bed, exc)
  if(exc >0) 
    return exc / bedObj.penDepth_ramp
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

#   return lspng * bedObj.stillWei  + 
#     # lspng * bedObj.kn * bedObj.od / bedObj.A * exc * sΛ -
#     # lspng * bedObj.linDampRatio* bedObj.kn * bedObj.od / bedObj.A * vz * sΛ
#     lspng * bedObj.cnstz * sΛ * ( exc - bedObj.linDampRatio * vz )

# end


# # Bed2
# function forceFnc(bedObj::Bed, X, QTr, T1s, T1m, u, ∇u, v)
  
#   local exc, lSpng, lStillWei
#   local FΓ, t1s, t1m2, sΛ        

#   exc = VectorValue(0.0,-1.0) ⋅ (X + u)
#   # lSpng = 0.5 + 0.5*( tanh( bedObj.tanh_ramp * exc ) )

#   lSpng = 0.0
#   if(exc >0) 
#     lSpng = exc / bedObj.penDepth_ramp
#   end
#   lStillWei = min(1.0, lSpng)
  
#   vz = VectorValue(0.0, 1.0) ⋅ v  

#   FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
#   t1s = FΓ ⋅ T1s
#   t1m2 = t1s ⋅ t1s    

#   sΛ = (t1m2.^0.5) / T1m
  
#   return lStillWei * bedObj.stillWei  + 
#     lSpng * bedObj.cnstz * sΛ * (exc - bedObj.linDampRatio * vz) 

# end


# Bed1b
function forceFnc(bedObj::Bed, X, QTr, T1s, T1m, u, ∇u, v)
  
  local exc, lSpng
  local FΓ, t1s, t1m2, sΛ        

  exc = VectorValue(0.0,-1.0) ⋅ (X + u)
  lSpng = rampTanh(bedObj, exc)
  lStillWei = min(1.0, lSpng)

  vz = VectorValue(0.0, 1.0) ⋅ v

  FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
  t1s = FΓ ⋅ T1s
  t1m2 = t1s ⋅ t1s    

  sΛ = (t1m2.^0.5) / T1m
  
  # bedObj.cnstz = bedObj.kn * bedObj.od / bedObj.A

  return lStillWei * bedObj.stillWei  + 
    lSpng * bedObj.cnstz * sΛ * ( 
      exc +
      -bedObj.linDampRatio * vz +
      -bedObj.quadDampRatio * vz * abs(vz) 
    )

end
# ----------------------End----------------------

end 