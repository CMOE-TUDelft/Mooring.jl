using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants
using .Currents

include(srcdir("gnlPara2D_moordyn.jl"))
include(srcdir("aux","gnlStructs.jl"))

resDir = datadir("sims_202409",
  "run")

ffm_η = 0.1 #m
ffm_f = 1.0 #Hz
ϵ0 = 0.0

dia = 0.048
A_str = π*dia*dia/4

# Warmup run
params = gnlPara2D.Test_params( 
  # initCSV = "models/catShape_xfl60_zfl20.csv",
  resDir = resDir,
  
  # Material properties
  E = 64.2986e9, #N
  L = 75, #m
  A_str = A_str, #m2 Str cross-section area
  ρcDry = 7.8e3, #kg/m3 Dry Density of steel   

  # Fairlead position
  xz_fl = (60, 20),
  
  # Parameter Domain
  nx = 100,
  order = 1,  

  # bedSpring setup
  bedObj = bedSpringStruct( dia, A_str ),

  outFreeSurface = false,

  # Time Parameters
  t0 = 0.0,
  simT = 20/ffm_f,
  simΔt = 1/ffm_f/100.0,
  outΔt = 1/ffm_f/4.0,
  maxIter = 300,

  # Drag coeff
  C_dn = 0.01, # Normal drag coeff
  d_dn = sqrt(4*dia/pi), #m Normal drag projection diameter
  C_dt = 0.0, # Tangent drag coff
  d_dt = sqrt(4*dia/pi), #m Tangent drag projection diameter
  
  # Time signal ramp up (t0 t1)
  startRamp = (0.0, 2/ffm_f),

  # Forced fairlead motion
  ffm_η = ffm_η, #m
  ffm_ω = 2*pi*ffm_f, #rad/s
  ϵ0 = ϵ0
)
gnlPara2D.main(params)


