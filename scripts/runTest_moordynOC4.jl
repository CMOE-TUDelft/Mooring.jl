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


resDir = datadir("sims_202501",
  # "run")
  "res_T12.1")
  # "res_T07.5")
  # "res_T20.0")
  # "res_T00.2")

ffm_η = 3.0 #m
ffm_f = 1/12.1 #Hz
ϵ0 = 0.0
# materialDampCoeff = 0.07
materialDampCoeff = 0.0
tStepsPerT = 100

# ffm_η = 0.5 #m
# ffm_f = 1/7.5 #Hz
# ϵ0 = 0.0
# # materialDampCoeff = 0.07
# materialDampCoeff = 0.0
# tStepsPerT = 100

# ffm_η = 10 #m
# ffm_f = 1/20 #Hz
# ϵ0 = 0.0
# # materialDampCoeff = 1e-1
# materialDampCoeff = 0.07
# tStepsPerT = 200

# ffm_η = 0.5 #m
# ffm_f = 5 #Hz 1/0.2 sec
# ϵ0 = 0.0
# # materialDampCoeff = 2e-3
# # materialDampCoeff = 0.07
# materialDampCoeff = 0.0007
# tStepsPerT = 100

dia = 0.0766
AStr = π*dia*dia/4

depth = 186


# Warmup run
params = gnlPara2D.Test_params( 
  # initCSV = "models/catShape_xfl60_zfl20.csv",
  resDir = resDir,
  
  # Material properties
  E = 1.63529e11, #N
  L = 835.35, #m
  AStr = AStr, #m2 Str cross-section area
  nd = dia,  #m Nominal diameter
  ρcDry = 24596.5, #kg/m3 Dry Density of steel   
  materialDampCoeff = materialDampCoeff, # Material damping coeff (sec)
  dragProp = Drag.DragProperties(
    Drag.Custom(), dia, AStr,
    Cd_n = 2.0, Cd_t = 0.8 ),

  # Fairlead position
  xz_fl = (796.732, depth),
  
  # # Parameter Domain
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
  simT = 2/ffm_f,
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

# # Production run
# params = gnlPara2D.Test_params(params, 
#   simT = 20/ffm_f)
# gnlPara2D.main(params)

"""
Other Setting

"""

# ffm_η = 3.0 #m
# ffm_f = 1/12.1 #Hz
# ϵ0 = 0.0

# dia = 0.0766*1.717 
# # 0.0766 is normail dia, this is equivalent dia
# # See Table 5 in Hall (2015)
# AStr = π*dia*dia/4

# depth = 186


# # Warmup run
# params = gnlPara2D.Test_params( 
#   # initCSV = "models/catShape_xfl60_zfl20.csv",
#   resDir = resDir,
  
#   # Material properties
#   E = 5.546932471920144e10, #N
#   L = 835.35, #m
#   AStr = AStr, #m2 Str cross-section area
#   nd = dia,  #m Nominal diameter
#   ρcDry = 8343.2165, #kg/m3 Dry Density of steel   
#   dragProp = Drag.DragProperties(
#     Drag.Custom(), dia, AStr,
#     Cd_n = 2.0, Cd_t = 0.8 ),

#   # Fairlead position
#   xz_fl = (796.732, depth),
  
#   # Parameter Domain
#   nx = 200,
#   order = 1,  

#   # bedSpring setup
#   bedObj = BedSpring.Bed( 
#     dia, AStr, 
#     kn = 3.0e6, dampRatio = 0.10 ),

#   outFreeSurface = false,

#   # Time Parameters
#   t0 = 0.0,
#   simT = 20/ffm_f,
#   simΔt = 1/ffm_f/100.0,
#   outΔt = 1/ffm_f/4.0,
#   maxIter = 300,
  
#   # Time signal ramp up (t0 t1)
#   inputRamp = TimeRampType(0.0, 2/ffm_f),  

#   Hs = 0.0,
#   h0=depth,

#   # Current
#   curObj = CurrentStat(depth, 
#     [-depth, -22.0, -15.0, 0.0], 
#     # [0.0, 0.2, 0.3, 0.3]),
#     [0.0, 0.0, 0.0, 0.0]),

  
#   FLMotion = FairLeadMotion.Regular(),

#   # Forced fairlead motion
#   ffm_η = ffm_η, #m
#   ffm_ω = 2*pi*ffm_f, #rad/s
#   ϵ0 = ϵ0
# )
# gnlPara2D.main(params)


