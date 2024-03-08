#!/bin/sh
#
#SBATCH --job-name="gnlMoor"
#SBATCH --partition=compute
#SBATCH --time=20:00:00
#SBATCH -n 8
#SBATCH -o stdout-batch/slurm-%j-%4t-%n.out
#SBATCH -e stdout-batch/slurm-%j-%4t-%n.err

source ./compile/modules.sh

INITIAL_CASE=1
FINAL_CASE=5
for i in $(seq $INITIAL_CASE $FINAL_CASE)
do
    echo "case: $i"
    export CASE_ID=$i
    # mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_serial.jl")' &
    srun -N1 -n1 -c1 --exact julia --project=. -O3 --check-bounds=no ./scripts/runConv.jl &
done
wait