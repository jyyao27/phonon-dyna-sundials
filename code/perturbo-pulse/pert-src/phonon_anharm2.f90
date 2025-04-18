module phonon_anharm2
   use hdf5_utils
   use boltz_grid, only: grid
   use pert_const, only: dp, pi, czero, cone, ryd2mev
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
      mp_split_pools, distribute_points
   use pert_param, only: delta_smear_ph, phfreq_cutoff, prefix, boltz_qdim
   use pert_utils, only: abs2, bose, fermi, gauss, match_table_all, find_free_unit, &
      get_exp_ikr
   use boltz_utils, only: kpt2num, num2kpt, kpts_minus, kpts_plus, kpts_plus_invert
   use pert_output,only: progressbar_init, progressbar, set_progress_step, stopwatch
   use phonon_dispersion, only: lattice_ifc, lattice_ifct, solve_phonon_modes, &
                solve_phi3mat, solve_gruneisen, solve_phi3mat_fast
   implicit none

   public :: calc_selfph, calc_gruneisen
contains

!subroutine calc_selfph (ph, pht, kpts, qset, qwt, tmpr, enk, imsgm)
subroutine calc_selfph(ph, pht, kpts, qg, tmpr, enk, imsgm)
   implicit none
   type(lattice_ifc), intent(in) :: ph
   type(lattice_ifct), intent(in) :: pht
   !tmpr: temperature array; 
   !qwt: q-points (qset) weight. kpts: k-points in the current pool. kpts(3,:)   
   !real(dp), intent(in) :: kpts(:,:), qset(:,:), qwt(:), tmpr(:)
   real(dp), intent(in) :: kpts(:,:), tmpr(:)
   type(grid), intent(in) :: qg

   !imsgm(nband, nmod, ntmpr, nkpt), enk(nband, nkpt)
   real(dp), intent(out) :: enk(:,:), imsgm(:,:,:,:) 
   ! local variables
   integer :: mem_max, step_max, nk_pool, qst_pool, qend_pool, nq_pool, nstep_k, ikc
   integer :: step_q, nstep_q, istep_k, kst, kend, step_k, n, icount, iter
   integer :: numk, numq, ntmpr, numm, ik, iq, it, im, ib, jb, i, ndim(4), nph3
   real(dp) :: xk(3), xq(3), xkq(3), wcut, e_thr,  ph_n2, ph_n3, dt1, dt2, scat
   
   real(dp) :: ek(ph%nm), wq(ph%nm), ekq(ph%nm), wq1, wq2, wq3
   complex(dp) :: ukq(ph%nm, ph%nm), uq(ph%nm, ph%nm)
   logical,  allocatable :: ltable(:,:,:)
   !complex(dp),  allocatable :: exp_ikr1(:), expiqr(:)
   integer :: num_kq_pair
   real(dp) :: g2
   complex(dp), allocatable :: uk(:,:,:)
   complex(dp) :: ph3
   integer, pointer :: ik_loc(:)
   integer :: iprocess

   call start_clock('calc_selfph')
   num_kq_pair = 0
   enk = 0.0_dp;  imsgm = 0.0_dp; ndim = shape(imsgm)
   !#. of selected bands, #.temperature, #. of kpts, number of q-points
   numm = ndim(1);  ntmpr = ndim(3);  numk = ndim(4);  !numq = size(qset,2)
   numq = qg%nk

   !array sanity check
   if(size(tmpr).ne.ntmpr .or. size(kpts,2).ne.numk &
      .or. size(enk,1).ne.numm .or. size(enk,2) .ne. numk) &
      call errore('calc_selfph','arguments dimension mismatch.',1)
    
   if(pht%nm .ne. ph%nm .or. pht%na .ne. ph%na) &
      call errore('calc_selfph','fc array dimension mismatch',1)

   e_thr = delta_smear_ph*3.0_dp  !exp(-3*3) ~ 10^-4
   wcut  = phfreq_cutoff
   nph3 =  (ph%na * 3)**3
   !max memory per process for g_kerp, hard-coded here, not sure what is the optimal value
   !mem_max  = 1024*1024*1024/4 !4GB
   mem_max  = 1024*1024*1024/2 !8GB
   step_max = mem_max / ( nph3 * pht%max_nrt )
   ! kelly - not importing pht
   !step_max = mem_max / ( nph3 * ph%max_nr )
   
   ik_loc => null()
   !if(numk > 2*step_max*npool) then
      !interleave distribution of kpoints over pools, for better load balance
      call distribute_points(numk, nk_pool, ik_loc, .true.)
     !call mp_split_pools(numk, kst_pool, kend_pool, nk_pool)
      qst_pool = 1;     qend_pool = numq;    nq_pool = numq
      
   !else
      !if (ionode) write(*,*) 'now spliting q pools'
   !   call distribute_points(numk, nk_pool, ik_loc, .false.)
      !kst_pool = 1;     kend_pool = numk;    nk_pool = numk
      !distrbite q among pools
   !   call mp_split_pools(numq, qst_pool, qend_pool, nq_pool)
      !if (ionode) write(*,*) qst_pool, qend_pool, nq_pool
   !endif
   
   !
   call set_progress_step(nk_pool, step_k, nstep_k, step_max)

   !if(nstep_k >= 8) then
   !   !tracking progess on k only
   !   step_q = nq_pool;    nstep_q = 1
   !else
      !track step on k and q combined
   !   call set_progress_step(nq_pool, step_q, nstep_q)
   !endif
   !output debug info
   write(stdout,'(5x, 5(a,i8,1x))') &
      'nq_pool',nq_pool,'nk_pool', nk_pool, 'step_k:', step_k, 'nstep_k:', nstep_k, 'step_max:', step_max
   
   !allocate work sapce
   allocate( uk(ph%nm, ph%nm, step_k))

   if(ionode) call progressbar_init('Im[Sigma]:')
   do istep_k = 1, nstep_k
      kst = (istep_k-1)*step_k + 1
      kend = min(kst+step_k-1, nk_pool)
      
      uk = czero;
!$omp parallel do schedule(guided) default(shared) private(n, ikc, ik, xk, ek)
      do n = kst, kend
         ikc = n - kst + 1
         !electronic wavefunction at ik
         ik = ik_loc(n)
         xk = kpts(:, ik)
         call solve_phonon_modes(ph, xk, ek, uk(:,:,ikc))
         
         enk(:,ik) = ek(:)
         
      enddo
!$omp end parallel do

   iter = 0
      ! 
! add iprocess to private
!$omp parallel default(shared) private(iq, xq, wq, wq1, wq2, wq3, uq, i, n, ik, ikc, xkq, &
!$omp& ekq, ukq, ltable, ph3, g2, im, ib, jb, dt1, dt2, it, ph_n2, ph_n3, &
!$omp& scat, iprocess) 
      allocate( ltable(ph%nm, ph%nm, ph%nm) )
      !allocate(exp_ikr1(pht%nrvect), expiqr(pht%nrvect))
!$omp do schedule(guided) 
      do iq = qst_pool, qend_pool
         ! set g2 as constant
         ! convert from mev**2 to ryd**2
         !g2 = 0.05/ryd2mev/ryd2mev
         g2 = 0.0_dp
         xq = num2kpt(qg%kpt(iq), qg%ndim)
         !xq = qset(1:3, iq)
         !get phonon frequcies and eigen-displacement at iq
         call solve_phonon_modes(ph, xq, wq, uq)
         !call get_exp_ikr(xq, pht%rvect_set(:,1,:), exp_ikr1)   
         !get e-ph matrix elements in wannier gauge and cartesian coord.
         do n = kst, kend
            ikc = n - kst + 1
            !phonon eigenmodes at ik
            ik = ik_loc(n)
            
            ! -xtong 
            ! q2      q0          q1
            xkq = - kpts(:,ik) - xq
            !write(*,*) 'ik, iq, ikq', ik, iq, kpt2num(xkq, qg%ndim)
            ! fold xkq in BZ
            ! I think xkq is in crystal coord
            ! if xkq(i) >= 0.5, xkq(i) = xkq(i) - 1
            
            ! recover full BZ kpt
            do i = 1, 3
               xkq(i) = merge(xkq(i)-1.0_dp, xkq(i), xkq(i) > 0.5_dp - 1.0E-10_dp)
               xkq(i) = merge(xkq(i)+1.0_dp, xkq(i), xkq(i) < -0.5_dp - 1.0E-10_dp)
            enddo

               !get electronic wavefunction at ikq
            call solve_phonon_modes(ph, xkq, ekq, ukq)
            !call get_exp_ikr(xkq, pht%rvect_set(:,2,:), expiqr)   
            !expiqr = expiqr * exp_ikr1
            
            ! -kelly using 3 different momentum conserving scheme
            do iprocess = 1, 3
               ! select the corresponding process
               !select case iprocess


               !call start_clock('match tab')
               !check if there is any scattering channel perserve energy conservation
               if (iprocess .eq. 1) then
                  call match_table_all(enk(:, ik), ekq, wq, wcut, e_thr, ltable)
               elseif (iprocess .eq. 2) then
                  call match_table_all(wq, enk(:, ik), ekq, wcut, e_thr, ltable)
               elseif (iprocess .eq. 3) then
                  call match_table_all(ekq, enk(:, ik), wq, wcut, e_thr, ltable)
               endif
               !call stop_clock('match tab')
               ! if everyone is zero
               if(.not. any(ltable)) then
                  cycle
               endif

               do im = 1, ph%nm
                  wq2 = wq(im)
                  if (wq2 .le. wcut) cycle
               do ib = 1, ph%nm
                  wq1 = enk(ib, ik)
                  if (wq1 .le. wcut) cycle
               do jb = 1, ph%nm
                  wq3 = ekq(jb)
                  if (wq3 .le. wcut) cycle

                  if (iprocess .eq. 1) then
                     if(.not. ltable(ib,jb,im)) cycle
                  elseif (iprocess .eq. 2) then
                     if(.not. ltable(im,ib,jb)) cycle
                  elseif (iprocess .eq. 3) then
                     if(.not. ltable(jb,ib,im)) cycle
                  endif

                  ph3 = czero;
                  !call solve_phi3mat(pht, kpts(:,ik), xq, xkq, enk(:,ik), wq, ekq, &
                  !                 uk(:,:,ik), uq, ukq, ib, im, jb, ph3)
               
                  ! - kelly: i think uk(:, :, ik) is wrong, here is the fix
                  call solve_phi3mat(pht, kpts(:,ik), xq, xkq, enk(:,ik), wq, ekq, &
                                   uk(:,:,ikc), uq, ukq, ib, im, jb, ph3)
                  !call solve_phi3mat_fast(pht, expiqr, wq1, wq2, wq3, &
                  !                 uk(:,ib,ikc), uq(:,im), ukq(:,jb), ph3)
               
                  !transfor to phonon mode and bloch gauge.
                  !compute |g|^2
               
                  g2 = abs2( ph3 )
                  
                  ! it=1
                  ! delete do it
                  do it = 1, ntmpr
                     ph_n2 = bose(tmpr(it), wq2) ! N_q\mu
                     ph_n3 = bose(tmpr(it), wq3)
                     !call start_clock('gauss')
                     !scat = 2.0_dp * (ph_n2 - ph_n3)*dt1 + (1.0_dp + ph_n2 + ph_n3)*dt2
                     if (iprocess .eq. 1) then
                        dt2 = gauss( delta_smear_ph, wq2 + wq3 - wq1) ! + freq_shift)
                        scat = (1.0_dp + ph_n2 + ph_n3)* dt2
                     elseif (iprocess .eq. 2) then
                        dt2 = gauss( delta_smear_ph, wq1 + wq3 - wq2) ! + freq_shift)
                        scat = (ph_n3 - ph_n2)* dt2
                     elseif (iprocess .eq. 3) then
                        dt2 = gauss( delta_smear_ph, wq1 + wq2 - wq3) ! + freq_shift)
                        scat = (ph_n2 - ph_n3)* dt2
                     endif
                     !call stop_clock('gauss')
                     !NOTE: the order of jb and ib are important! Tr[g^\dagger g]
                     scat = 18.0_dp * scat * pi * qg%kweight * g2

!!NB: only one thread can update this shared variable at a time. negligible overhead.
!$omp atomic update
                     imsgm(ib, im, it, ik) = imsgm(ib, im, it, ik) + scat
                  enddo
               !if (ionode) write(*,*) dt1, d2, scat, ph_n2, phn3
               enddo; enddo; enddo
            enddo
         enddo

         !track progress
!$omp atomic update
         iter = iter + 1
         
         if( (mod(iter, int(0.1*nstep_k*nq_pool)).eq.0) .or. (iter.eq. int(0.1*nstep_k*nq_pool)) ) then
!$omp critical (calc_selfph_progress)
            write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/(1.0_dp*nq_pool*nstep_k), '%'
!            write(stdout,'(8x, f7.2, a1)') (iter*1.0_dp/nstep_k), '%'
!$omp end critical (calc_selfph_progress)
         endif

      enddo
!$omp end do
      deallocate(ltable)
!$omp end parallel
 enddo

!   write(stdout,*) 'num_kq_pair', num_kq_pair

   deallocate(uk)
   call mp_sum(imsgm, inter_pool_comm)
   if(numk .ne. nk_pool) call mp_sum(enk, inter_pool_comm)
   !release space
   if( associated(ik_loc) ) deallocate( ik_loc )

   call stop_clock('calc_selfph')

end subroutine calc_selfph


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
   wcut  = phfreq_cutoff
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

end module phonon_anharm2
