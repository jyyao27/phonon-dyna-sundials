module phonon_anharm_prop
   use hdf5_utils
   use boltz_grid, only: grid
   use pert_const, only: dp, pi, czero, cone, ryd2mev
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
      mp_split_pools, distribute_points
   use pert_param, only: delta_smear_ph, phfreq_cutoff_ph, prefix, boltz_qdim
   use pert_utils, only: abs2, bose, fermi, gauss, match_table_all, find_free_unit
   use boltz_utils, only: kpt2num, num2kpt, kpts_minus, kpts_plus
   use pert_output,only: progressbar_init, progressbar, set_progress_step, stopwatch
   use phonon_dispersion, only: lattice_ifc, lattice_ifct, solve_phonon_modes, &
                solve_phi3mat, solve_gruneisen
   implicit none

   public :: calc_gruneisen
contains

subroutine calc_gruneisen (ph, pht, kpts, wq, grun)
   implicit none
   type(lattice_ifc), intent(in) :: ph
   type(lattice_ifct), intent(in) :: pht
   !kpts: k-points in the current pool. kpts(3,:)
   real(dp), intent(in) :: kpts(:,:)
   !grun(nmodes, nkpt), wq(nmodes, nkpt)
   real(dp), intent(out) :: wq(:,:), grun(:,:) 
   ! local variables
   integer :: mem_max, step_max, nk_pool, nstep_k, ikc
   integer :: istep_k, kst, kend, step_k, n, icount, iter 
   integer :: numk, nmodes, ik, nph3, ndim(2)
   real(dp) :: xk(3), wcut, e_thr
   real(dp) :: ek(ph%nm), g(ph%nm)
   complex(dp), allocatable :: uk(:,:,:)
   integer, pointer :: ik_loc(:)

   call start_clock('calc_gruneisen')
   wq = 0.0_dp;  grun = 0.0_dp; ndim = shape(grun)
   !#. of phonon modes, #. of kpts
   nmodes = ndim(1);  numk = ndim(2);
   !array sanity check
   if(size(kpts,2).ne.numk .or. size(wq,1).ne.nmodes .or. size(wq,2) .ne. numk) &
      call errore('calc_gruneisen','arguments dimension mismatch.',1)
   if(pht%nm .ne. ph%nm .or. pht%na .ne. ph%na) &
      call errore('calc_gruneisen','fc array dimension mismatch',1)

   e_thr = delta_smear_ph*3.0_dp  !exp(-3*3) ~ 10^-4
   wcut  = phfreq_cutoff_ph
   nph3 =  (ph%na * 3)**3
   !max memory per process for g_kerp, hard-coded here, not sure what is the optimal value
   mem_max  = 1024*1024*1024/2 !8GB
   step_max = mem_max / ( nph3 * pht%max_nrt )
   
   ik_loc => null()
!   if(numk > 2*step_max*npool) then
      !interleave distribution of kpoints over pools, for better load balance
   call distribute_points(numk, nk_pool, ik_loc, .true.)
   !
   call set_progress_step(nk_pool, step_k, nstep_k, step_max)
   !output debug info
   !write(stdout,'(5x, 3(a,i8,1x))') &
   !   'step_k:', step_k, 'nstep_k:', nstep_k, 'step_max:', step_max
   
   !allocate work sapce
   allocate( uk(ph%nm, ph%nm, step_k))

   if(ionode) call progressbar_init('Gruneisen:')
   icount = 0
   do istep_k = 1, nstep_k
      kst = (istep_k-1)*step_k + 1
      kend = min(kst+step_k-1, nk_pool)
      
      uk = czero;
!$omp parallel do schedule(guided) default(shared) private(n, ikc, ik, xk, ek, g)
      do n = kst, kend
         ikc = n - kst + 1
         !electronic wavefunction at ik
         ik = ik_loc(n)
         xk = kpts(:, ik)
         call solve_phonon_modes(ph, xk, ek, uk(:,:,ikc))
         call solve_gruneisen(pht, xk, ek, uk(:,:,ikc), g)
         wq(:,ik) = ek
         grun(:,ik) = g
      enddo
!$omp end parallel do
      !track progress

         icount = icount + 1
         iter = icount
         if( mod(iter, nstep_k).eq.0 .or. iter.eq.(nstep_k) ) then
            write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/(nstep_k), '%'
         endif

   enddo

   deallocate(uk)
   if(numk .ne. nk_pool) then
          call mp_sum(grun, inter_pool_comm)
          call mp_sum(wq, inter_pool_comm)
   endif
   !release space
   if( associated(ik_loc) ) deallocate( ik_loc )

   call stop_clock('calc_gruneisen')

end subroutine calc_gruneisen

end module phonon_anharm_prop
