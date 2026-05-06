program gpu_acc
    use mpi
    implicit none
    !include "mpif.h"
    ! Problem variables
    logical, parameter :: do_cpu_control=.false.
    integer, parameter :: N = 2000000, M = 6, L = 1, Niter = 10000
    real, allocatable :: ens(:,:,:)
    real, allocatable :: new(:)
    integer, allocatable :: sample(:)
    integer :: i, j, k
    real(kind=8) :: val
    ! Timing variables
    real :: start_cpu, end_cpu, start_gpu, end_gpu
    ! MPI variables
    integer :: jproc, jpproc, mpi_code
    integer :: jproc0, valsize

    ! Initialize MPI
    CALL mpi_init(mpi_code)
    CALL mpi_comm_size(mpi_comm_world,jpproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_world,jproc,mpi_code)
    print *, 'MPI initialized, proc=', jproc

    ! Broadcast test
    val=0.
    if (jproc==0) val=1.
    valsize = 1 ; jproc0 = 0
    CALL mpi_bcast(val,valsize,mpi_double_precision,jproc0,mpi_comm_world,mpi_code)
    print *, 'Value, proc=', jproc,val

    allocate(ens(N, M, L))
    allocate(new(N))
    allocate(sample(M))

    ! Initialize ens and sample arrays
    ! Fill ens with some values and set sample
    do i = 1, N
    do j = 1, M
       ens(i, j, :) = real(i+j,4)/real(N+M,4)
    end do
    end do

    do j = 1, M  ! M assumed even
       if (mod(j,2) == 0) then
         sample(j) = j/2
       else
         sample(j) = (M+j+1)/2
       endif
    enddo

    ! CPU timing
     call cpu_time(start_cpu)
     do k = 1, Niter
       call product_cpu(ens,sample,new,N,M,L)
     end do
     call cpu_time(end_cpu)
     print *, 'Proc: ,',jproc,' CPU Time: ', end_cpu - start_cpu
     print *, 'Value: ', minval(new(:)), maxval(new(:)),sum(new(:)),new(1)

    ! GPU timing
    call cpu_time(start_gpu)
    !$acc data copyin(ens, sample) create(new)
    !$omp target data map(to:ens,sample) map(alloc:new)
    do k = 1, Niter
       call product_gpu(ens,sample,new,N,M,L)
    end do
    !$omp end target data
    !$acc end data
    call cpu_time(end_gpu)
    print *, 'Proc: ,',jproc,' GPU Time: ', end_gpu - start_gpu
    print *, 'Value: ', minval(new(:)), maxval(new(:)),sum(new(:)),new(1)

    ! CPU timing
    if (do_cpu_control) then
      call cpu_time(start_cpu)
      do k = 1, Niter
        call product_cpu(ens,sample,new,N,M,L)
      end do
      call cpu_time(end_cpu)
      print *, 'Proc: ,',jproc,' CPU Time: ', end_cpu - start_cpu
      print *, 'Value: ', minval(new(:)), maxval(new(:)),sum(new(:)),new(1)
    endif

    deallocate(ens)
    deallocate(new)

    ! Finalize MPI
    CALL mpi_finalize(mpi_code)

end program gpu_acc

subroutine product_cpu(ens,sample,new,N,M,L)
    use mpi
    implicit none
    real, intent(in) :: ens(N,M,L)
    integer, intent(in) :: sample(M)
    real, intent(out) :: new(N)
    integer, intent(in) :: N, M, L

    integer :: i, j
    real(kind=8) :: a, b, cost_jobs, cost_one

    integer :: mpi_code
    real(kind=8) :: val

    a = 0.3 ; b = 0.7

    val=1.
    CALL mpi_bcast(val,1,mpi_double_precision,0,mpi_comm_world,mpi_code)

    ! Compute Schur products
    ! ----------------------
    do i = 1, N
      new(i) = ens(i, sample(1), 1)
      do j = 2, M
        new(i) = new(i) * ens(i, sample(j), 1)
      end do
    end do

    ! Apply Schur product
    ! -------------------
    do i = 1, N
      new(i) = a * ens(i,sample(1),1) + b * new(i)
    end do

    ! Evaluation of cost function
    ! ---------------------------
    cost_jobs = 0.
    do i = 1, N
      cost_one  = ( new(i) - ens(i,sample(M),1) ) / b
      cost_jobs = cost_jobs + 0.5 * cost_one * cost_one
    end do

    new(1) = cost_jobs

end subroutine product_cpu

subroutine product_gpu(ens,sample,new,N,M,L)
    use mpi
    implicit none
    real, intent(in) :: ens(N,M,L)
    integer, intent(in) :: sample(M)
    real, intent(out) :: new(N)
    integer, intent(in) :: N, M, L

    integer :: i, j
    real :: a, b, cost_jobs, cost_one

    integer :: mpi_code
    real(kind=8) :: val

    a = 0.3 ; b = 0.7

    val=1.
    CALL mpi_bcast(val,1,mpi_double_precision,0,mpi_comm_world,mpi_code)

    ! Compute Schur products
    ! ----------------------
    !$acc data copyin(sample) present(ens, new)
    !$acc parallel loop
    !$omp target update to(sample)
    !$omp target teams distribute parallel do
    do i = 1, N
      new(i) = ens(i, sample(1), 1)
      do j = 2, M
        new(i) = new(i) * ens(i, sample(j), 1)
      end do
    end do
    !$omp end target teams distribute parallel do
    !$acc end parallel loop
    !$acc end data

    ! Apply Schur product
    ! -------------------
    !$acc data present(ens, sample, new)
    !$acc parallel loop
    !$omp target teams distribute parallel do
    do i = 1, N
      new(i) = a * ens(i,sample(1),1) + b * new(i)
    end do
    !$omp end target teams distribute parallel do
    !$acc end parallel loop
    !$acc end data

    ! Evaluation of cost function
    ! ---------------------------
    cost_jobs = 0.
    !$acc data present(ens, sample, new)
    !$acc parallel loop private(cost_one) reduction(+:cost_jobs)
    !$omp target teams distribute parallel do reduction(+:cost_jobs) private(cost_one)
    do i = 1, N
      cost_one  = ( new(i) - ens(i,sample(M),1) ) / b
      cost_jobs = cost_jobs + 0.5 * cost_one * cost_one
    end do
    !$omp end target teams distribute parallel do
    !$acc end parallel loop
    !$acc end data

    new(1) = cost_jobs

end subroutine product_gpu

