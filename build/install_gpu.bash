#!/bin/bash
#

cluster="jean-zay"
cluster="adastra"

if [ $cluster = "jean-zay" ] ; do

  module purge
  module load nvidia-compilers/24.3
  module load cuda/12.2.0  # not needed before, worked with cuda/12.4.1
  module load openmpi/4.0.5-cuda
  module load netcdf-fortran/4.5.3-mpi-cuda
  #module load hdf5/1.14.5-mpi-cuda # not needed before (+hdf5_fortran in Makefile.macro)
  module load hdf5/1.12.0-mpi-cuda # not needed before (+hdf5_fortran in Makefile.macro)

  ln -sf macro/make.jean-zay_gpu Makefile.macro

elif [ $cluster = "adastra" ] ; do

  module purge
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed

  ln -sf macro/make.adastra_gpu Makefile.macro

else

  echo "Bad cluster name"
  exit

fi

target="$HOME/bin/sesam_gpu"

./mkmf -t Makefile.macro -p $target ../src/*.[Ffh]90

make

#make install

