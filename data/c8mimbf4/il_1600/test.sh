#!/bin/sh
#SBATCH --partition=pre         # default "univ" if not specified
#SBATCH --time=0-05:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1          # require 2 nodes
#SBATCH --ntasks-per-node=20            # default 16 if this line not specified
#SBATCH --mem-per-cpu=4000      # RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


module load compile/intel-2016
#module load intel-2017.2.0 

export OMP_NUM_THREADS=20
export OMP_STACKSIZE=200M
#cd /home/cson2/program/utility/trj_analysis/ionicLiquid/structure_factor
#bash compile_omp.sh
#env
#ldd calc_struct_omp

#cd data/c8mimbf4/il_1600
export FOR_PRINT='test.out'
export FORT0='test.err'
../../../calc_struct_omp -f param.dat
