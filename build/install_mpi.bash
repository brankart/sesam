#!/bin/bash
#

cluster='jean-zay'
cluster='adastra'

if [ $cluster = 'jean-zay' ] ; then

  module purge
  module load intel-all/19.0.4
  module load hdf5/1.8.21-mpi
  module load netcdf/4.7.2-mpi
  module load netcdf-fortran/4.5.2-mpi
  module load netcdf-cxx4/4.3.1-mpi

  ln -sf ../macro/make.jean-zay_mpiifort Makefile.macro

elif [ $cluster = 'adastra' ] ; then

  module purge
  module load cpe/24.07
  module load craype-x86-genoa
  module load PrgEnv-cray
  module load cray-mpich
  module load cray-hdf5-parallel
  module load cray-netcdf-hdf5parallel

  ln -sf ../macro/make.adastra_genoa Makefile.macro

else

  echo "Bad cluster name"
  exit

fi


target="$HOME/bin/sesam_mpi"

./mkmf -t Makefile.macro -p $target ../src/*.[Ffh]90

make

#make install

