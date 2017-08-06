#!/bin/bash
#SBATCH --partition=pre         # default "univ" if not specified
#SBATCH --time=0-03:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1          # require 2 nodes
#SBATCH --ntasks-per-node=16            # default 16 if this line not specified
#SBATCH --mem-per-cpu=4000      # RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

module load compile/intel-2016
#module load intel-2017.2.0 

export OMP_NUM_THREADS=16
export OMP_STACKSIZE=500M
sleep 10800

#../../../calc_struct_omp -f param.dat
