module boltz_scatter_sizest
   use kinds, only: i8b
   use iso_c_binding, only: c_size_t
   use pert_const, only: dp

   use qe_mpi_mod, only: stdout, ionode, npool, my_pool_id

   use yaml_utils, only: ymlout
   use sys_utils, only: bytes_to_gib

   implicit none
   public

   ! Total bytes allocated for Boltzmann scattering calculations code
   integer(kind=i8b), save :: boltz_scatter_totalbytes_cpu = 0

   ! Total bytes allocated on GPU
   integer(kind=i8b), save :: boltz_scatter_totalbytes_gpu = 0

contains


subroutine boltzscat_accum_cpu_bytes(count)
   implicit none
   integer(kind=i8b), intent(in) :: count

   boltz_scatter_totalbytes_cpu = boltz_scatter_totalbytes_cpu + count
end subroutine


subroutine boltzscat_accum_gpu_bytes(count)
   implicit none
   integer(kind=i8b), intent(in) :: count

   boltz_scatter_totalbytes_gpu = boltz_scatter_totalbytes_gpu + count
end subroutine


subroutine boltzscat_accum_cpugpu_bytes(count)
   implicit none
   integer(kind=i8b), intent(in) :: count

   boltz_scatter_totalbytes_cpu = boltz_scatter_totalbytes_cpu + count
   boltz_scatter_totalbytes_gpu = boltz_scatter_totalbytes_gpu + count
end subroutine


subroutine boltz_scatter_estimate_pools(max_cpumem_gb, max_gpumem_gb)
   implicit none
   !
   real(dp), intent(in) :: max_cpumem_gb, max_gpumem_gb
   !
   real(dp) :: cpumem_gb_pool, cpumem_gb_total
   real(dp) :: gpumem_gb_pool, gpumem_gb_total
   integer :: pools_cpu, pools_gpu, pools_est
   character(len=40) :: pool_limit

   ! Approximate the amount of memory needed for both the CPU data structures
   ! and the GPU data structures.  We compute the total amounts so we can
   ! guess how many pools we might need.

   cpumem_gb_pool = bytes_to_gib(boltz_scatter_totalbytes_cpu)
   cpumem_gb_total = cpumem_gb_pool * real(npool)

   gpumem_gb_pool = bytes_to_gib(boltz_scatter_totalbytes_gpu)
   gpumem_gb_total = gpumem_gb_pool * real(npool)

   pools_cpu = ceiling(cpumem_gb_total / max_cpumem_gb)
   pools_gpu = ceiling(gpumem_gb_total / max_gpumem_gb)

   if (npool < pools_cpu .or. npool < pools_gpu) then
      write (stdout, *) 'WARNING:  number of pools ', npool, ' may be too low; OOM is likely'
   endif

   ! Figure out an estimate for the number of pools (aka. MPI tasks) required.
   ! Increasing the number of pools decreases the memory required per pool, so
   ! we need the MAXIMUM pool-count predicted from CPU or GPU memory usage.
   if (pools_cpu > pools_gpu) then
      pools_est = pools_cpu
      pool_limit = '(CPU memory limited)'
   elseif (pools_gpu > pools_cpu) then
      pools_est = pools_gpu
      pool_limit = '(GPU memory limited)'
   else
      ! CPU and GPU requirements seem to be the same
      pools_est = pools_cpu
      pool_limit = '(CPU+GPU memory limited)'
   endif

   if (pools_est == 0) pools_est = 1

   write (stdout, '(5X,A,I5,2X,A)') 'Estimated minimum pools required:  ', pools_est, trim(pool_limit)

end subroutine boltz_scatter_estimate_pools


end module boltz_scatter_sizest

