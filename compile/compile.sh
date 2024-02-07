source ./modules.sh

julia --project=.. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
julia --project=.. -O3 --check-bounds=no --color=yes -e "include("./warmup.jl")"