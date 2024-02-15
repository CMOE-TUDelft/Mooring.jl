#!/bin/bash

#SBATCH --job-name="compile"
#SBATCH -p compute
#SBATCH -t 00:50:00
#SBATCH -n 1
#SBATCH -o stdout
#SBATCH -e stderr

source ./modules.sh

srun julia --project=.. -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
srun julia --project=.. -O3 --check-bounds=no --color=yes -e 'include("warmup.jl")'