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
  L = 75, #m
  A_str = 2*(π*0.048*0.048/4), #m2 Str cross-section area
  ρcDry = 7.8e3, #kg/m3 Dry Density of steel

  # Parameter Domain
  nx = 100,

  # Time Parameters
  t0 = 0.0,
  simT = 10.0,
  simΔt = 0.2,
  outΔt = 0.2,

  # Fairlead Excitation
  fairLead_η = 0.1,
  fairLead_T = 4.0,

  # Drag coeff
  Cdn = 2.6, # Normal drag coeff
  d_dn = 0.048, #m Normal drag projection diameter
  Cdt = 1.4, # Tangent drag coff
  d_dt = 0.048/pi #m Trangent drag projection diameter
)
gnlPara2D.main(params)


