#!/bin/bash
#

cluster="jean-zay"
cluster="adastra"

if [ $cluster = "jean-zay" ] ; then

  module purge
  module load nvidia-compilers/24.3
  module load cuda/12.2.0  # not needed before, worked with cuda/12.4.1
  module load openmpi/4.0.5-cuda
  module load netcdf-fortran/4.5.3-mpi-cuda
  #module load hdf5/1.14.5-mpi-cuda # not needed before (+hdf5_fortran in Makefile.macro)
  module load hdf5/1.12.0-mpi-cuda # not needed before (+hdf5_fortran in Makefile.macro)

  ln -sf ../macro/make.jean-zay_gpu Makefile.macro

elif [ $cluster = "adastra" ] ; then

  module purge
  #module load cpe/24.07
  module load cpe/25.09
  module load craype-x86-trento
  module load craype-accel-amd-gfx90a
  module load PrgEnv-cray
  module load rocm
  module load cray-mpich
  module load cray-hdf5-parallel
  module load cray-netcdf-hdf5parallel

  ln -sf ../macro/make.adastra_gpu Makefile.macro

else

  echo "Bad cluster name"
  exit

fi

target="$HOME/bin/sesam_gpu"

./mkmf -t Makefile.macro -p $target ../src/*.[Ffh]90

make

#make install

