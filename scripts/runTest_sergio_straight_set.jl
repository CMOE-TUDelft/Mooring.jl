using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants
using .Currents
using TickTock
using Printf

include(srcdir("gnlPara2D_sergio.jl"))


allparams = Dict(
  "epsilon0" => 0.01,
  "ffm_amp" => [0.4],
  "ffm_f" => [0.2], #Hz  
)

dicts = dict_list(allparams)


function makesim(d::Dict)  

  # resDir = datadir("sims_sergio_202406","run")

  ffm_η = d["ffm_amp"] #m
  ffm_f = d["ffm_f"] #Hz
  ϵ0 = d["epsilon0"]
  @show ffm_f

  
  folderName = @sprintf("runSet_epsilon0=%0.02f_ffm_amp=%0.02f",
    ϵ0, ffm_η)

  #caseName = savename(d)
  caseName = @sprintf("epsilon0=%0.02f_ffm_amp=%0.02f_ffm_f=%0.02f",
    ϵ0, ffm_η, ffm_f)
  println(caseName)
  
  resDir = datadir(
    "sims_202410",
    folderName,
    caseName)
  mkdir(resDir)  

  A_str = 0.01
  nd = sqrt(4*A_str/π)

  # Warmup run
  params = gnlPara2D.Test_params( 
    initCSV = "models/initVertical.csv",
    resDir = resDir,

    # Material properties
    E = 1e9, #N
    L = 100, #m
    A_str = A_str, #m2 Str cross-section area
    ρcDry = 1414.5, #kg/m3 Density of steel  
    
    # Parameter Domain
    nx = 100,
    order = 2,

    outFreeSurface = false,

    # Time Parameters
    t0 = 0.0,
    simT = 200/ffm_f,
    simΔt = 1/ffm_f/100.0,
    outΔt = 10/ffm_f,
    maxIter = 300,

    h0 = 1000,
    
    # Time signal ramp up (t0 t1)
    startRamp = (0.0, 10/ffm_f),

    # Forced fairlead motion
    ffm_η = ffm_η, #m
    ffm_ω = 2*pi*ffm_f, #rad/s
    ϵ0 = ϵ0,

    # Drag coeff
    C_dn = 0.01, # Normal drag coeff
    d_dn = nd, #m Normal drag projection diameter
    C_dt = 0.0, # Tangent drag coff
    d_dt = nd/pi #m Tangent drag projection diameter

    #no Bed 
    #no wave
    #yes Self Drag
  )
  gnlPara2D.main(params)
end


for (i, d) in enumerate(dicts)
  tick()
  makesim(d)   
  tock()
end 