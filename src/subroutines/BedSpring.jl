module BedSpring

using Revise
using Gridap



"""
Custom Structs
=============

"""
# ---------------------Start---------------------
mutable struct Bed
  kn::Real 
  # Based on Marco's suggestion of 30kN/m2/m in the OrcaFlex manual  
  od::Real
  A::Real  
  tanh_ramp::Real
  stillWei::Real
  dampRatio::Real

  cnst1::Real

  function Bed( od::Real, A::Real,     
    kn = 30e3, tanh_ramp = 1e3,
    dampRatio = 0.05,
    stillWei::Real=0.0 )    

    cnst1 = kn * od / A
    
    new(kn, od, A, tanh_ramp, stillWei, dampRatio, cnst1)
  end
end
# ----------------------End----------------------



"""
Functions
=============

"""
# ---------------------Start---------------------
function forceFnc(bedObj, X, QTr, T1s, T1m, u, ∇u, v)
  
  local exc, lspng
  local FΓ, t1s, t1m2, sΛ        

  exc = VectorValue(0.0,-1.0) ⋅ (X + u)
  lspng = 0.5 + 0.5*( tanh( bedObj.tanh_ramp * exc ) )

  vz = VectorValue(0.0, 1.0) ⋅ v

  FΓ = ( ∇u' ⋅ QTr ) + TensorValue(1.0,0.0,0.0,1.0)
  t1s = FΓ ⋅ T1s
  t1m2 = t1s ⋅ t1s    

  sΛ = (t1m2.^0.5) / T1m
  
  # bedObj.cnst1 = bedObj.kn * bedObj.od / bedObj.A

  return lspng * bedObj.stillWei  + 
    # lspng * bedObj.kn * bedObj.od / bedObj.A * exc * sΛ -
    # lspng * bedObj.dampRatio* bedObj.kn * bedObj.od / bedObj.A * vz * sΛ
    lspng * bedObj.cnst1 * sΛ * ( exc - bedObj.dampRatio * vz )

end
# ----------------------End----------------------

end 