export GMXXTC=$HOME/program/utility/trajmod/gmx_xtc
export OPENMMDCD=$HOME/program/utility/trajmod/openmm_dcd

#gfortran -openmp $GMXXTC/xtc-interface-wrap.o variables.f90 topol.f90 fvector.f90 filehandle.f90 conductivity.f90 -o calc_cond_omp -lxdrfile -L $GMXXTC/lib -I $GMXXTC -CB
#gfortran -fopenmp $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 filehandle.f90 structure_factor.f90 -o calc_struct_omp -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -ffree-line-length-256
#gfortran -fopenmp $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 filehandle.f90 sq_from_gr.f90 -o sq_from_gr  -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -ffree-line-length-256
#gfortran -fopenmp $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 filehandle.f90 qq_sq_from_gr.f90 -o qq_sq_from_gr  -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -ffree-line-length-256
gfortran -fopenmp $GMXXTC/xtc-interface-wrap.o $OPENMMDCD/dcdmod.o $OPENMMDCD/trajmod.o variables.f90 topol.f90 fvector.f90 filehandle.f90 qq_read_frame.f90 qq_structure_factor.f90 -o calc_qq_struct_omp -lxdrfile -L $GMXXTC/lib -I $GMXXTC -I $OPENMMDCD -ffree-line-length-256
