#!/bin/bash

# Test case to compile
example="test2"
example="test5_mpi"

# Type of compilation
compile="crayftn"        # Good
compile="ifort"          # Slower (factor 5 with cray)
compile="gfortran"       # Very slow (factor 20 with cray!)
compile="flang"          # Very slow (factor 20 with cray!)
compile="omp_crayftn"    # Good
compile="acc_crayftn"    # Slow (factor 10 with omp!)
compile="omp_flang"      # Very slow (factor 20 with cray!)
compile="mpi_crayftn"    # Good (factor 1/8 with jean-zay?)
compile="mpiomp_crayftn" # Good

module purge
# --------------
# No MPI, no GPU
# --------------
if [ $compile = "ifort" ] ; then
  module load PrgEnv-intel
  ./$example
elif [ $compile = "gfortran" ] ; then
  module load PrgEnv-gnu
  ./$example
elif [ $compile = "crayftn" ] ; then
  module load PrgEnv-cray
  ./$example
elif [ $compile = "flang" ] ; then
  module load PrgEnv-aocc
  ./$example
# ----------------
# MPI only, no GPU
# ----------------
elif [ $compile = "mpi_crayftn" ] ; then
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed
  srun -A hmg2840 --constraint=GENOA --nodes=1 --time=00:10:00 --ntasks-per-node=2 --hint=nomultithread ./$example
# ----------------
# GPU only, no MPI
# ----------------
elif [ $compile = "acc_crayftn" ] ; then
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed
  srun -A hmg2840 --constraint=MI250 --nodes=1 --time=00:10:00 --ntasks-per-node=1 --cpus-per-task=1 --gres=gpu:1 --hint=nomultithread ./$example
elif [ $compile = "omp_crayftn" ] ; then
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed
  srun -A hmg2840 --constraint=MI250 --nodes=1 --time=00:10:00 --ntasks-per-node=1 --cpus-per-task=1 --gres=gpu:1 --hint=nomultithread ./$example
elif [ $compile = "omp_flang" ] ; then
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-amd
  module load amd-mixed
  srun -A hmg2840 --constraint=MI250 --nodes=1 --time=00:10:00 --ntasks-per-node=1 --cpus-per-task=1 --gres=gpu:1 --hint=nomultithread ./$example
# --------------------
# MPI and GPU together
# --------------------
elif [ $compile = "mpiomp_crayftn" ] ; then
  module load cpe/24.07
  module load craype-accel-amd-gfx90a craype-x86-trento
  module load PrgEnv-cray
  module load amd-mixed
  #srun -A hmg2840 --constraint=MI250 --nodes=1 --time=00:10:00 --ntasks-per-node=2 --cpus-per-task=1 --gres=gpu:2 --hint=nomultithread ./$example
  srun -A hmg2840 --constraint=MI250 --nodes=1 --time=00:10:00 --ntasks-per-node=1 --cpus-per-task=1 --gres=gpu:1 --hint=nomultithread ./$example
else
  echo "Bad compilation option"
fi

