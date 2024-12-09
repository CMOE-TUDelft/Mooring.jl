using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants
using .Currents
using TickTock
using Printf

include(srcdir("gnlPara2D_conv.jl"))

allparams = Dict(
  "epsilon0" => 0.01,
  "ffm_amp" => [0.0],
  "ffm_f" => [0.5], #Hz  
  "nx" => [2,4,8,16,32,64],
  "TBydt" => [2000000],
  "order" => [1]
)
folderName = @sprintf("conv_h_o1")


# allparams = Dict(
#   "epsilon0" => 0.01,
#   "ffm_amp" => [0.0],
#   "ffm_f" => [0.5], #Hz  
#   "nx" => [2,4,8,16,32,64],
#   "TBydt" => [2000000],
#   "order" => [2]
# )
# folderName = @sprintf("conv_h_o2")


# allparams = Dict(
#   "epsilon0" => 0.01,
#   "ffm_amp" => [0.0],
#   "ffm_f" => [0.5], #Hz  
#   "nx" => [2,4,8,16,32,64],
#   "TBydt" => [2000000],
#   "order" => [3]
# )
# folderName = @sprintf("conv_h_o3")


# allparams = Dict(
#   "epsilon0" => 0.01,
#   "ffm_amp" => [0.0],
#   "ffm_f" => [0.5], #Hz  
#   "nx" => [2],
#   "TBydt" => [2000000],
#   "order" => collect(range(1,10))
# )
# folderName = @sprintf("conv_p")


dicts = dict_list(allparams)


function makesim(d::Dict)  

  # resDir = datadir("sims_sergio_202406","run")

  ffm_η = d["ffm_amp"] #m
  ffm_f = d["ffm_f"] #Hz
  ϵ0 = d["epsilon0"]
  nx = d["nx"]
  ndt = d["TBydt"]
  order = d["order"]
  @show ffm_f

  
  # folderName = @sprintf("runSet_epsilon0=%0.02f_ffm_amp=%0.02f",
  #   ϵ0, ffm_η)
  # folderName = @sprintf("conv_p")

  caseName = savename(d)
  # caseName = @sprintf("epsilon0=%0.02f_ffm_amp=%0.02f_ffm_f=%0.02f_TBydt=%04i_order=%02i",
  #   ϵ0, ffm_η, ffm_f,ndt,order)
  println(caseName)
  
  resDir = datadir(
    "sims_202412","run",
    folderName,
    caseName)
  mkdir(resDir)  

  # Warmup run
  params = gnlPara2D.Test_params( 
    initCSV = "models/initVertical.csv",
    resDir = resDir,

    # Material properties
    E = 1e9, #N
    L = 100, #m
    A_str = 0.01, #m2 Str cross-section area
    ρcDry = 1414.5, #kg/m3 Density of steel  
    
    xz_fl = (0, 100),
    
    # Parameter Domain
    nx = nx,
    order = order,

    outFreeSurface = false,

    # Time Parameters
    t0 = 0.0,
    simT = 2/ffm_f/ndt,
    simΔt = 1/ffm_f/ndt,
    outΔt = 1000/ffm_f/ndt,
    maxIter = 300,

    h0 = 1000,
    
    # Time signal ramp up (t0 t1)
    # startRamp = (0.0, 2/ffm_f),
    startRamp = (0.0, 1/ffm_f/4),

    uhAna_amp = 0.1,
    # uhAna_ampForCmp = 0.1*0.9963389469, #0.10 1e4 
    # uhAna_ampForCmp = 0.1*0.999963368, #0.10 1e3
    # uhAna_ampForCmp = 0.1*0.9996492768, #0.01 1e4
    uhAna_ampForCmp = 0.1*0.999996493, #0.01 1e3

    # Forced fairlead motion
    ffm_η = ffm_η, #m
    ffm_ω = 2*pi*ffm_f, #rad/s
    ϵ0 = ϵ0
  )
  gnlPara2D.main(params)

  # Production run
  params = gnlPara2D.Test_params(params, 
    # simT = 2/ffm_f
    simT = 1000/ffm_f/ndt
    )
  gnlPara2D.main(params)
end


for (i, d) in enumerate(dicts)
  tick()
  makesim(d)   
  tock()
end 