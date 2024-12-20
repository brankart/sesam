#!/bin/bash
#

module purge
module load intel-all/19.0.4
module load hdf5/1.8.21-mpi
module load netcdf/4.7.2-mpi
module load netcdf-fortran/4.5.2-mpi
module load netcdf-cxx4/4.3.1-mpi

ln -sf ../macro/make.jean-zay_mpiifort Makefile.macro

target="$HOME/bin/sesam_mpi"

./mkmf -t Makefile.macro -p $target ../src/*.[Ffh]90

make

#make install

