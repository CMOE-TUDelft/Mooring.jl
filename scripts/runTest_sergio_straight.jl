using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants
using .Currents

include(srcdir("gnlPara2D_sergio.jl"))

resDir = datadir("sims_sergio_202407",
  "run")

ffm_η = 0.1 #m
ffm_f = 0.80 #Hz
ϵ0 = 0.01

# Warmup run
params = gnlPara2D.Test_params( 
  initCSV = "models/initVertical.csv",
  resDir = resDir,

  # # Material properties
  # E = 100e9, #N
  # L = 10, #m
  # A_str = 0.01, #m2 Str cross-section area
  # ρcDry = 980, #kg/m3 Density of steel  

  # Material properties
  E = 1e9, #N
  L = 100, #m
  A_str = 0.01, #m2 Str cross-section area
  ρcDry = 1414.5, #kg/m3 Density of steel  
  
  # Parameter Domain
  nx = 200,
  order = 1,

  outFreeSurface = false,

  # Time Parameters
  t0 = 0.0,
  simT = 100/ffm_f,
  simΔt = 1/ffm_f/80.0,
  outΔt = 1/ffm_f,
  maxIter = 300,

  # Drag coeff
  C_dn = 0.1, # Normal drag coeff
  d_dn = sqrt(4*0.01/pi), #m Normal drag projection diameter
  C_dt = 0.0, # Tangent drag coff
  d_dt = sqrt(4*0.01/pi), #m Tangent drag projection diameter
  
  # Time signal ramp up (t0 t1)
  startRamp = (0.0, 10/ffm_f),

  # Forced fairlead motion
  ffm_η = ffm_η, #m
  ffm_ω = 2*pi*ffm_f, #rad/s
  ϵ0 = ϵ0
)
gnlPara2D.main(params)


