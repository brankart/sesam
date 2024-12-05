#!/bin/bash
#

module purge
module load nvidia-compilers/24.3
module load cuda/12.2.0  # not needed before, worked with cuda/12.4.1
module load openmpi/4.0.5-cuda
module load netcdf-fortran/4.5.3-mpi-cuda

ln -sf ../macro/make.jean-zay_gpu Makefile.macro

target="$HOME/bin/sesam_gpu"

./mkmf -t Makefile.macro -p $target ../src/*.[Ffh]90

make

#make install

