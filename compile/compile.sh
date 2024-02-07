source ./modules.sh

julia --project=.. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"