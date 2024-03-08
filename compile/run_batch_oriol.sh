#!/bin/sh
#
#SBATCH --job-name="gnlMoor"
#SBATCH --partition=genoa
#SBATCH --time=2-00:00:00
#SBATCH -n 576
#SBATCH -o stdout-batch/slurm-%j-%4t-%n.out
#SBATCH -e stdout-batch/slurm-%j-%4t-%n.err

source ../compile/modules_snellius.sh
INITIAL_CASE=1
FINAL_CASE=588
for i in $(seq $INITIAL_CASE $FINAL_CASE)
do
    echo "case: $i"
    export CASE_ID=$i
    # mpiexecjl --project=../ -n 1 julia -J ../PerforatedCylinder_parallel.so -O3 --check-bounds=no -e 'include("run_case_serial.jl")' &
    srun -N1 -n1 -c1 --exact julia --project=../ -J ../PerforatedCylinder_parallel_genoa.so -O3 --check-bounds=no -e 'include("run_case_serial.jl")' &
done
wait