using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants
using .Currents

include(srcdir("gnlPara2D_simple.jl"))

resDir = datadir("sims_202406","run")

ffm_η = 0.1 #m
ffm_f = 0.2 #Hz

# Warmup run
params = gnlPara2D.Test_params( 
  initCSV = "models/initShape.csv",
  resDir = resDir,

  # Material properties
  E = 99.82198e9, #N
  L = 75, #m
  A_str = π*0.05*0.05/4, #m2 Str cross-section area
  ρcDry = 7.8e3, #kg/m3 Density of steel  
  
  # Parameter Domain
  nx = 80,
  order = 1,

  outFreeSurface = false,

  # Time Parameters
  t0 = 0.0,
  simT = 100/ffm_f,
  simΔt = 0.05,
  outΔt = 1.25,
  maxIter = 100,
  
  # Time signal ramp up (t0 t1)
  startRamp = (0.0, 10/ffm_f),

  # Forced fairlead motion
  ffm_η = ffm_η, #m
  ffm_ω = 2*pi*ffm_f #rad/s

)
gnlPara2D.main(params)


