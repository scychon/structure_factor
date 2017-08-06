#!/bin/sh
#SBATCH --partition=pre         # default "univ" if not specified
#SBATCH --time=0-05:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1          # require 2 nodes
#SBATCH --ntasks-per-node=20            # default 16 if this line not specified
#SBATCH --mem-per-cpu=4000      # RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#module load compile/intel
#module load mpi/intel/mvapich2-1.9

# Set some environment variables 
#export PATH=$PROGRAM/gmx-4.6.5/bin:$PATH

module load intel-2017.2.0
#module load intel-2017.2.0 

export OMP_NUM_THREADS=20
export OMP_STACKSIZE=2000M

#export FOR_PRINT='test.out'
#export FORT0='test.err'
#ldd ../../../calc_struct_omp
../../../calc_struct_omp -f param.dat
