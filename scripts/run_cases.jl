module RunCases
using Mooring
using DrWatson
@quickactivate "Mooring.jl"

# Warmup run
params = Mooring.gnl_simple_params(nx=3,sampling_points=[0.5])
Mooring.main_static(params)
Mooring.main_dynamic(params)

# Test cases
η = [ηi for ηi in 0.05:0.05:0.3]
ω = [ωi for ωi in 0.2:0.2:2.0]
all_params = @strdict η ω
cases = dict_list(all_params)

# Execute case function
function run_case(case)
  @unpack η, ω = case
  case_name = savename(case)
  params = Mooring.gnl_simple_params(η₀=η,ω=ω,testname=case_name,T=200.0,Δt=2.0,outΔt=2.0)
  Mooring.main_dynamic(params)
end

# Execute cases
params = Mooring.gnl_simple_params() # initial condition
Mooring.main_static(params)
for case in cases
  # path = datadir()
  # data,file = produce_or_load(path,case,run_case)
  run_case(case)
end

end
