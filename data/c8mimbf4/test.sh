#module load compile/intel-2016
module load intel-2017.2.0 

export OMP_NUM_THREADS=16
export OMP_STACKSIZE=200M
../../calc_struct_omp -f param.dat
