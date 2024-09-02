"""
Custom Structs
=============

"""
# ---------------------Start---------------------
mutable struct bedSpringStruct
  kn::Real 
  # Based on Marco's suggestion of 30kN/m2/m in the OrcaFlex manual  
  od::Real
  A::Real  
  tanh_ramp::Real
  stillWei::Real
  dampRatio::Real

  function bedSpringStruct( od::Real, A::Real,     
    kn = 30e3, tanh_ramp = 1e3,
    dampRatio = 0.05,
    stillWei::Real=0.0 )    
    
    new(kn, od, A, tanh_ramp, stillWei, dampRatio)
  end

end
# ----------------------End----------------------

