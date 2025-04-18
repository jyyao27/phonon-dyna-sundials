program crash
   implicit none
   integer, parameter :: dp = selected_real_kind(14,200)
   real(dp), allocatable :: arr(:)

!$omp parallel default(shared) private(arr)
   allocate( arr(50) )
   deallocate( arr )
!$omp end parallel

end program crash
