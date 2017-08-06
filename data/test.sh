module load compile/intel-2016
#module load intel-2017.2.0 

export OMP_NUM_THREADS=10
export OMP_STACKSIZE=500M
../calc_struct_omp -f param.dat
