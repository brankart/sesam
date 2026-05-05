#!/bin/bash

# Test case to compile
example="test2"
example="test5_mpi"

# Type of compilation
compile="ifort"          # Good
compile="gfortran"       # Good
compile="crayftn"        # Good
compile="flang"          # Good
compile="acc_crayftn"    # Good
compile="omp_crayftn"    # Good
compile="omp_flang"      # Good
compile="mpi_crayftn"    # Good
compile="mpi_crayftn"    # Good
compile="mpiomp_crayftn" # Good

# purge all modules
module purge

# --------------
# No MPI, no GPU
# --------------
if [ $compile = "ifort" ] ; then        # Good
  module load PrgEnv-intel
  ifort $example.f90 -o $example
elif [ $compile = "gfortran" ] ; then   # Good
  module load PrgEnv-gnu
  gfortran $example.f90 -o $example
elif [ $compile = "crayftn" ] ; then    # Good
  module load PrgEnv-cray
  ftn $example.f90 -o $example
elif [ $compile = "flang" ] ; then      # Good
  module load PrgEnv-aocc
  flang $example.f90 -o $example
# ----------------
# MPI only, no GPU
# ----------------
elif [ $compile = "mpi_crayftn" ] ; then  # Good
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed
  ftn -hnoacc $example.f90 -o $example
# ----------------
# GPU only, no MPI
# ----------------
elif [ $compile = "acc_crayftn" ] ; then    # Good
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed
  ftn -hacc -hnoomp $example.f90 -o $example
elif [ $compile = "omp_crayftn" ] ; then    # Good
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed
  ftn -homp -hnoacc $example.f90 -o $example
elif [ $compile = "omp_flang" ] ; then    # Good
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-amd
  module load amd-mixed
  ftn -fopenmp $example.f90 -o $example
# --------------------
# MPI and GPU together
# --------------------
elif [ $compile = "mpiomp_crayftn" ] ; then # Fail
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed
  ftn -homp -hnoacc $example.f90 -o $example
else
  echo "Bad compilation option"
fi

