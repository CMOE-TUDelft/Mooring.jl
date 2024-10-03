using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants
using .Currents
using .WaveTimeSeries
using Mooring.BedSpring
using Mooring.Drag
using Mooring.FairLeadMotion

include(srcdir("gnlPara2D_moordyn.jl"))


resDir = datadir("sims_202409",
  "run")

ffm_η = 0.1 #m
ffm_f = 0.50 #Hz
ϵ0 = 0.0

dia = 0.048
AStr = π*dia*dia/4


# Warmup run
params = gnlPara2D.Test_params( 
  # initCSV = "models/catShape_xfl60_zfl20.csv",
  resDir = resDir,
  
  # Material properties
  E = 64.2986e9, #N
  L = 75, #m
  AStr = AStr, #m2 Str cross-section area
  nd = dia,  #m Nominal diameter
  ρcDry = 7.8e3, #kg/m3 Dry Density of steel   
  dragProp = Drag.DragProperties(
    Drag.Custom(), dia, AStr,
    Cd_n = 0.051503226936425275, Cd_t = 0.0 ),

  # Fairlead position
  xz_fl = (60, 20),
  
  # Parameter Domain
  nx = 400,
  order = 1,  

  # bedSpring setup
  bedObj = BedSpring.Bed( dia, AStr ),

  outFreeSurface = false,

  # Time Parameters
  t0 = 0.0,
  simT = 20/ffm_f,
  simΔt = 1/ffm_f/100.0,
  outΔt = 1/ffm_f/4.0,
  maxIter = 300,
  
  # Time signal ramp up (t0 t1)
  inputRamp = TimeRampType(0.0, 2/ffm_f),  

  Hs = 0.0,

  # Current
  curObj = CurrentStat(23, 
    [-23.0, -22.0, -15.0, 0.0], 
    # [0.0, 0.2, 0.3, 0.3]),
    [0.0, 0.0, 0.0, 0.0]),

  
  FLMotion = FairLeadMotion.Regular(),

  # Forced fairlead motion
  ffm_η = ffm_η, #m
  ffm_ω = 2*pi*ffm_f, #rad/s
  ϵ0 = ϵ0
)
gnlPara2D.main(params)

# # Production run
# params = gnlPara2D.Test_params(params, 
#   simT = 20/ffm_f)
# gnlPara2D.main(params)


