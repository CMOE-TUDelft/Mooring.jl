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
using Mooring.StressNLVE

include(srcdir("gnlPara2D_moordyn.jl"))


resDir = datadir("sims_202501",
  "run")
  # "run_sch_veStiff_f3p05")

ffm_η = 1.0 #m
ffm_f = 3.05 #Hz
ϵ0 = 0.01
# materialDampCoeff = 0.07
materialDampCoeff = 0.0
tStepsPerT = 100


## StressNLVE
# ---------------------Start---------------------   
# sch = StressNLVE.Schapery(true,
#   D0 = 1.97e-10,
#   Dn = [  1e-10,  1.5e-10,  1e-10,  1.5e-10, 
#           20e-10, 1e-10,    1e-10,  100e-10],
#   λn = [  1e-1,   1e-2,     1e-3,   1e-4,
#           1e-5,   1e-6,     1e-7,   1e-8],
#   g0 = [0.55, -1e-10, -3e-18, 1e-26],
#   g1 = [0.10, 6.6e-11, 4.2e-19, -2.7e-27],
#   g2 = [1.7, -8.8e-9, 4.4e-17, 5.0e-28]
# )

sch = StressNLVE.Schapery(true,
  D0 = 1.97e-10,
  Dn = [0.0, 0.0, 0.0, 0.0],
  λn = [1e-1, 1e-2, 1e-3, 1e-4],
  g0 = [1, 0],
  g1 = [0, 0],
  g2 = [0, 0]
)
# ----------------------End----------------------  

dia = 0.052
AStr = π*dia*dia/4

depth = 110
L0 = 100


# Warmup run
params = gnlPara2D.Test_params( 
  # initCSV = "models/catShape_xfl60_zfl20.csv",
  resDir = resDir,
  
  # Material properties
  E = 5.076e9, #N
  L = L0, #m
  AStr = AStr, #m2 Str cross-section area
  nd = dia,  #m Nominal diameter
  ρcDry = 1380, #kg/m3 Dry Density of steel   
  materialDampCoeff = materialDampCoeff, # Material damping coeff (sec)
  dragProp = Drag.DragProperties(
    Drag.Custom(), dia, AStr,
    Cd_n = 0.01, Cd_t = 0 ),

  # Schapery characteristics
  sch = sch,

  # Fairlead position
  xz_fl = (0.0, L0),
  
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


