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
# using Mooring.StressNLVE

include(srcdir("gnlPara2D_moordyn.jl"))


ffm_η = 3.0 #m
ffm_f = 1/12.1 #Hz
ϵ0 = 0.0
# materialDampCoeff = 0.07
materialDampCoeff = 0.0
tStepsPerT = 100

caseName = savename(@dict ϵ0 ffm_f ffm_η; digits=3)
typeName = "LE"

resDir = datadir("sims_202501","catenary",
  # "run")  
  caseName, typeName)

dia = 0.052
AStr = π*dia*dia/4

depth = 186


# Warmup run
params = gnlPara2D.Test_params( 
  # initCSV = "models/catShape_xfl60_zfl20.csv",
  resDir = resDir,
  
  # Material properties
  E = 5.076e9, #N
  L = 835.35, #m
  AStr = AStr, #m2 Str cross-section area
  nd = dia,  #m Nominal diameter
  ρcDry = 1380, #kg/m3 Dry Density of steel   
  materialDampCoeff = materialDampCoeff, # Material damping coeff (sec)
  dragProp = Drag.DragProperties(
    Drag.Custom(), dia, AStr,
    Cd_n = 0.01, Cd_t = 0 ),

  # Fairlead position
  xz_fl = (796.732, depth),
  
  # Parameter Domain
  nx = 100,
  order = 2,
  # nx = 24,
  # order = 10,  

  # bedSpring setup
  bedObj = BedSpring.Bed( 
    dia, AStr, 
    kn = 3.0e6, 
    linDampRatio = 0.10,
    quadDampRatio = 0.10,
    tanh_ramp = 1e2),

  outFreeSurface = false,

  # Time Parameters
  t0 = 0.0,
  simT = 20/ffm_f,
  simΔt = 1/ffm_f/tStepsPerT,
  outΔt = 1/ffm_f/4.0,
  maxIter = 100,
  
  # Time signal ramp up (t0 t1)
  inputRamp = TimeRampType(0.0, 5/ffm_f),  
  
  h0 = depth,
  enableWaveSpec = false,

  # Current
  curObj = CurrentStat(depth, 
    [-depth, -22.0, -15.0, 0.0], 
    # [0.0, 0.2, 0.3, 0.3]),
    [0.0, 0.0, 0.0, 0.0]),

  
  FLMotion = FairLeadMotion.Regular(),

  # Forced fairlead motion
  ffm_η = ffm_η, #m
  ffm_ω = 2*pi*ffm_f, #rad/s
  ϵ0 = ϵ0
)
gnlPara2D.main(params)

