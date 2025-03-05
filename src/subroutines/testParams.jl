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
  AStr = π*0.048*0.048/4 #m2 Str cross-section area
  nd = 0.048  #m Nominal diameter
  ρcDry = 7.8e3 #kg/m3 Density of steel    
  materialDampCoeff = 0.0 # Material damping coeff (sec)
  dragProp = Drag.DragProperties(Drag.NoDrag()) 

  # Schapery characteristics
  sch = StressNLVE.Schapery(true,
    D0 = 1/64.2986e9,
    Dn = [0.0, 0.0, 0.0, 0.0],
    λn = [1.0, 1e-1, 1e-2, 1e-3],
    g0 = [1, 0.0],
    g1 = [0.0, 0.0],
    g2 = [0.0, 0.0]
  )

  xz_fl = (60, 20)

  # Parameter Domain
  nx = 100
  order = 1  

  # bedSpring setup
  bedObj = BedSpring.Bed( 0.048, π*0.048*0.048/4 )
  
  outFreeSurface = false

  # Time Parameters
  t0 = 0.0
  simT = 0.04
  simΔt = 0.02
  outΔt = 0.2
  maxIter = 200

  # Time signal ramp up (t0 t1)
  inputRamp = TimeRampType(0.0, 15)

  # Wave spectrum
  enableWaveSpec = true
  h0 = 23 #m
  Hs = 3 #m
  Tp = 12 #s
  nω = 257 #including 0
  seed = -1
  ωc = -1

  # Current
  curObj = CurrentStat(23, [-23.0, -11.0, 0.0], [0.0, 0.0, 0.0])

  FLMotion::FairLeadMotion.MotionType = FairLeadMotion.Regular()

  # Forced fairlead motion
  ffm_η = 0.1 #m
  ffm_ω = 0.5 #Hz
  ϵ0 = 0.1

end