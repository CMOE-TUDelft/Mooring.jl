using DrWatson
@quickactivate "Mooring.jl"

using Parameters
using WaveSpec
using .Constants
using CSV
using DataFrames
using Printf

include(srcdir("gnlPara2D_rieke.jl"))

testID = parse(Int,ENV["CASE_ID"])
setList = CSV.read(scriptsdir("meshConvSetup.csv"), DataFrame)
resDir = datadir("sims", "runBatch", @sprintf("case_%04i",testID))
mkpath(resDir)

@show testID, setList[testID,:]
@show resDir

# Warmup run
params = gnlPara2D.Test_params( 
  resDir = resDir,
  
  # Parameter Domain
  nx = 20,
  order = 1,

  # Time Parameters
  t0 = 0.0,
  simT = 0.04,
  simΔt = 0.02,
  outΔt = 0.2,
  maxIter = 2

)
gnlPara2D.main(params)


# Production run
params = gnlPara2D.Test_params( 

  initCSV = "models/catShape_xfl60_zfl20.csv",
  # resDir = datadir("sims","testBedProperties","run_tanh1000_k0100_damp05_pPos"),
  resDir = resDir,

  # Material properties
  E = 64.2986e9, #N
  L = 75, #m
  A_str = 2*(π*0.048*0.048/4), #m2 Str cross-section area
  ρcDry = 7.8e3, #kg/m3 Dry Density of steel

  # Parameter Domain
  nx = setList[testID,2],
  order = setList[testID,3]

  # Time Parameters
  t0 = 0.0,
  simT = 15.0,
  simΔt = 0.01,
  outΔt = 0.25,
  
  # Drag coeff
  C_dn = 2.6, # Normal drag coeff
  d_dn = 0.048, #m Normal drag projection diameter
  C_dt = 1.4, # Tangent drag coff
  d_dt = 0.048/pi, #m Trangent drag projection diameter

  # Added mass coeff
  C_an = 1.0, # Normal added-mass coeff
  C_at = 0.5 # Tangent added-mass coff
)
gnlPara2D.main(params)


