#!/bin/bash

# Test case to compile
example="test5_mpi"

# Type of compilation
compile="ifort"        # Good
compile="pgi"          # Good
compile="nvidia"       # Good
compile="gpu_pgi"      # Fail
compile="gpu_nvi"      # Good
compile="mpi_ifort"    # Good
compile="mpi_pgi"      # Slow
compile="mpi_nvi"      # Slow
compile="mpigpu_pgi"   # Fail
compile="mpigpu_nvi"   # Good

# purge all modules
module purge

# --------------
# No MPI, no GPU
# --------------
if [ $compile = "ifort" ] ; then        # Good
  module load intel-compilers/19.0.4
  ifort $example.f90 -o $example
elif [ $compile = "pgi" ] ; then        # Good
  module load pgi/20.4
  pgfortran $example.f90 -o $example
elif [ $compile = "nvidia" ] ; then     # Good
  module load nvidia-compilers/24.3
  nvfortran $example.f90 -o $example
# ----------------
# MPI only, no GPU
# ----------------
elif [ $compile = "mpi_ifort" ] ; then  # Good
  module load intel-all/19.0.4
  mpiifort $example.f90 -o $example
elif [ $compile = "mpi_pgi" ] ; then    # Slow
  module load pgi/20.4
  module load openmpi/4.0.5-cuda
  mpifort $example.f90 -o $example
elif [ $compile = "mpi_nvi" ] ; then    # Slow
  module load nvidia-compilers/24.3
  module load openmpi/4.0.5-cuda
  mpifort $INC $LIB $example.f90 -o $example
# ----------------
# GPU only, no MPI
# ----------------
elif [ $compile = "gpu_pgi" ] ; then    # Fail
  module load pgi/20.4
  #pgfortran -ta=tesla:cc70 $example.f90 -o $example
  #pgfortran -acc=autopar,sync $example.f90 -o $example
  pgfortran -acc $example.f90 -o $example
elif [ $compile = "gpu_nvi" ] ; then    # Good
  module load nvidia-compilers/24.3
  nvfortran -acc -Minfo=accel $example.f90 -o $example
# --------------------
# MPI and GPU together
# --------------------
elif [ $compile = "mpigpu_pgi" ] ; then # Fail
  module load pgi/20.4
  module load openmpi/4.0.5-cuda
  mpifort -ta=tesla:cc70 -Minfo=accel $example.f90 -o $example
elif [ $compile = "mpigpu_nvi" ] ; then # Good
  module load nvidia-compilers/24.3
  module load openmpi/4.0.5-cuda
  module load netcdf-fortran/4.5.3-mpi-cuda
  mpifort -acc -Minfo=accel $example.f90 -o $example
else
  echo "Bad compilation option"
fi

