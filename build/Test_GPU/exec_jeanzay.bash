#!/bin/bash

# Test case to compile
example="test3_mpi"

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

module purge
# --------------
# No MPI, no GPU
# --------------
if [ $compile = "ifort" ] ; then
  module load intel-compilers/19.0.4
  ./$example
elif [ $compile = "pgi" ] ; then
  module load pgi/20.4
  ./$example
elif [ $compile = "nvidia" ] ; then
  module load nvidia-compilers/24.3
  ./$example
# ----------------
# MPI only, no GPU
# ----------------
elif [ $compile = "mpi_ifort" ] ; then
  module load intel-all/19.0.4
  srun -A egi@cpu --ntasks=2 --hint=nomultithread ./$example
elif [ $compile = "mpi_pgi" ] ; then
  module load pgi/20.4
  module load openmpi/4.0.5-cuda
  #export MPI_HOME=/gpfslocalsys/pgi/20.4/linux86-64-llvm/20.4/mpi/openmpi-3.1.3
  #export PATH=$MPI_HOME/bin:$PATH
  #export LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH
  srun -A egi@cpu --ntasks=2 --hint=nomultithread ./$example --mca mpi_cuda_support 0
elif [ $compile = "mpi_nvi" ] ; then
  module load nvidia-compilers/24.3
  export LD_LIBRARY_PATH=/gpfslocalsys/nvhpc/24.3/Linux_x86_64/24.3/comm_libs/12.3/hpcx/hpcx-2.17.1/ompi/lib:$LD_LIBRARY_PATH
  srun -A egi@cpu --ntasks=2 --hint=nomultithread ./$example
# ----------------
# GPU only, no MPI
# ----------------
elif [ $compile = "gpu_pgi" ] ; then
  module load pgi/20.4
  srun -A egi@v100 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --gres=gpu:1 --hint=nomultithread ./$example
elif [ $compile = "gpu_nvi" ] ; then
  module load nvidia-compilers/24.3
  srun -A egi@v100 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --gres=gpu:1 --hint=nomultithread ./$example
# --------------------
# MPI and GPU together
# --------------------
elif [ $compile = "mpigpu_pgi" ] ; then
  module load pgi/20.4
  srun -A egi@v100 --nodes=1 --ntasks-per-node=2 --cpus-per-task=1 --gres=gpu:2 --hint=nomultithread ./$example
elif [ $compile = "mpigpu_nvi" ] ; then
  module load nvidia-compilers/24.3
  module load openmpi/4.0.5-cuda
  module load netcdf-fortran/4.5.3-mpi-cuda
  srun -A egi@v100 --nodes=1 --ntasks-per-node=2 --cpus-per-task=1 --gres=gpu:2 --hint=nomultithread ./$example
else
  echo "Bad compilation option"
fi

