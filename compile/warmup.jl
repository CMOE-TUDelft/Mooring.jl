module Warmup
using Mooring
using DrWatson
@quickactivate "Mooring.jl"

params = Mooring.gnl_simple_params(nx=3)
Mooring.main_static(params)
end
