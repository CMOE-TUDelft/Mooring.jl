using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants

include(srcdir("gnlPara2D_rieke.jl"))

# Warmup run
params = gnlPara2D.Warmup_params( 
  resDir = datadir("sims","run")
)
gnlPara2D.main(params)


# Production run
params = gnlPara2D.Warmup_params( 

  initCSV = "models/catShape_xfl60_zfl20.csv",
  resDir = datadir("sims","run"),

  # Material properties
  E = 64.2986e9, #N
  mDry = 52.8,  #kg/m Dry weight per unit len
  mSub = 45.936,  #kg/m Submerged weight per unit len
  L = 75, #m
  ρc = 7.8e3, #kg/m3 Density of steel

  # Parameter Domain
  nx = 100,

  # Time Parameters
  t0 = 0.0,
  simT = 20.0,
  simΔt = 0.2,
  outΔt = 0.2,

  # Fairlead Excitation
  fairLead_η = 0.1,
  fairLead_T = 4.0
)
gnlPara2D.main(params)


