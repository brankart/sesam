program gpu_acc
    !use openacc
    implicit none
    integer, parameter :: N = 1000000, M = 6, L = 1, Niter = 1000
    real, allocatable :: ens(:,:,:)
    real, allocatable :: new(:)
    integer, allocatable :: sample(:)
    integer :: i, j, k
    real :: start_cpu, end_cpu, start_gpu, end_gpu

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
    print *, 'CPU Time: ', end_cpu - start_cpu
    print *, 'Value: ', minval(new(:)), maxval(new(:)),sum(new(:))

    ! GPU timing
    call cpu_time(start_gpu)
    !$acc data copyin(ens, sample) create(new)
    !$omp target data map(to:ens,sample) map(alloc:new)
    do i = 1, Niter
       call product_gpu(ens,sample,new,N,M,L)
    end do
    !$omp end target data
    !$acc end data
    call cpu_time(end_gpu)
    print *, 'GPU Time: ', end_gpu - start_gpu
    print *, 'Value: ', minval(new(:)), maxval(new(:)),sum(new(:))

    ! CPU timing
    call cpu_time(start_cpu)
    do i = 1, Niter
      call product_cpu(ens,sample,new,N,M,L)
    end do
    call cpu_time(end_cpu)
    print *, 'CPU Time: ', end_cpu - start_cpu
    print *, 'Value: ', minval(new(:)), maxval(new(:)),sum(new(:))

    deallocate(ens)
    deallocate(new)
end program gpu_acc

subroutine product_cpu(ens,sample,new,N,M,L)
    implicit none
    real, intent(in) :: ens(N,M,L)
    integer, intent(in) :: sample(M)
    real, intent(out) :: new(N)
    integer, intent(in) :: N, M, L

    integer :: i, j

    do i = 1, N
      new(i) = 1.
      do j = 1, M
        new(i) = new(i) * ens(i, sample(j), 1)
      end do
    end do

end subroutine product_cpu

subroutine product_gpu(ens,sample,new,N,M,L)
    implicit none
    real, intent(in) :: ens(N,M,L)
    integer, intent(in) :: sample(M)
    real, intent(out) :: new(N)
    integer, intent(in) :: N, M, L

    integer :: i, j

    !$acc data present(ens, sample, new)
    !$acc parallel loop
    !! OMP5.x : $omp target teams distribute parallel do map(present:ens,sample,new)
    !$omp target teams distribute parallel do
    do i = 1, N
      new(i) = 1.
      do j = 1, M
        new(i) = new(i) * ens(i, sample(j), 1)
      end do
    end do
    !$omp end target teams distribute parallel do
    !$acc end parallel loop
    !$acc end data

end subroutine product_gpu

