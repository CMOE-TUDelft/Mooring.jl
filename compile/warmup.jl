module Warmup
using Mooring
using DrWatson
@quickactivate "Mooring.jl"

params = Mooring.gnl_simple_params(nx=3,sampling_points=[0.5,0.75])
Mooring.main_static(params)
Mooring.main_dynamic(params)
end
