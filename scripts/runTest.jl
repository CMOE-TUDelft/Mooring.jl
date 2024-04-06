using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants
using .Currents

include(srcdir("gnlPara2D_rieke.jl"))

# Warmup run
params = gnlPara2D.Test_params( 
  resDir = datadir("sims","run"),
  
  # Parameter Domain
  nx = 5,
  order = 1,

  # Time Parameters
  t0 = 0.0,
  simT = 0.04,
  simΔt = 0.02,
  outΔt = 0.2,
  maxIter = 2

)
gnlPara2D.main(params)


# Sea state
h0 = 22.81 #m
Hs = 3.8179 #m
Tp = 8.37 #s
nω = 257 #including 0

strCur = CurrentStat(22.81, 
  [-22.81, -21.81, -16.81, -11.81, -5.81, -0.81, 0.0], 
  [0.0, 0.3104, 0.4269, 0.487, 0.4906, 0.5781, 0.5781])
@show strCur.itp

# Production run
params = gnlPara2D.Test_params( 

  initCSV = "models/catShape_xfl60_zfl20.csv",
  # resDir = datadir("sims",
  #   "testCurrent",
  #   "run_cur1_dt0p10_k50p0_rr"),
  resDir = datadir("sims","run"),

  # Material properties
  E = 64.2986e9, #N
  L = 75, #m
  A_str = 2*(π*0.048*0.048/4), #m2 Str cross-section area
  ρcDry = 7.8e3, #kg/m3 Dry Density of steel  

  # Parameter Domain
  nx = 100,

  outFreeSurface = true,

  # Time Parameters
  t0 = 0.0,
  simT = 50.0,
  simΔt = 0.10,
  outΔt = 0.2,
  
  # Drag coeff
  C_dn = 2.6, # Normal drag coeff
  d_dn = 0.048, #m Normal drag projection diameter
  C_dt = 1.4, # Tangent drag coff
  d_dt = 0.048/pi, #m Trangent drag projection diameter

  # Added mass coeff
  C_an = 1.0, # Normal added-mass coeff
  C_at = 0.5, # Tangent added-mass coff

  # Wave spectrum
  h0 = h0,
  Hs = Hs,
  Tp = Tp,
  nω = nω,
  seed = 1234,

  # Current
  strCur = strCur
)
gnlPara2D.main(params)


