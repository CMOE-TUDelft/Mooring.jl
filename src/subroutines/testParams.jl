"""
Warmup and Test params
===========

"""
@with_kw struct Test_params

  ρw = 1025 #Kg/m3 density of water    

  # initCSV::String = "models/catShape_xfl60_zfl20.csv" 
  resDir = "data/"

  # Material properties
  E = 64.2986e9 #N
  L = 75 #m
  A_str = π*0.048*0.048/4 #m2 Str cross-section area
  ρcDry = 7.8e3 #kg/m3 Density of steel    

  xz_fl = (60, 20)

  # Parameter Domain
  nx = 100
  order = 1  

  # bedSpring setup
  bedObj = bedSpring.Bed( 0.048, π*0.048*0.048/4 )
  
  outFreeSurface = false

  # Time Parameters
  t0 = 0.0
  simT = 0.04
  simΔt = 0.02
  outΔt = 0.2
  maxIter = 200

  # Drag coeff
  C_dn = 2.6 # Normal drag coeff
  d_dn = 0.048 #m Normal drag projection diameter
  C_dt = 1.4 # Tangent drag coff
  d_dt = 0.048/pi #m Tangent drag projection diameter

  # Added mass coeff
  C_an = 1.0 # Normal added-mass coeff
  C_at = 0.5 # Tangent added-mass coff

  # Time signal ramp up (t0 t1)
  startRamp = (0.0, 15)

  # Wave spectrum
  h0 = 23 #m
  Hs = 3 #m
  Tp = 12 #s
  nω = 257 #including 0
  seed = -1
  ωc = -1

  # Current
  strCur = CurrentStat(23, [-23.0, -11.0, 0.0], [0.0, 0.0, 0.0])

  # Forced fairlead motion
  ffm_η = 0.1 #m
  ffm_ω = 0.5 #Hz
  ϵ0 = 0.1

end