!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
! Copyright (C) 2021-2023 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Ivan Maliyov
!                         Dhruv Desai, Sergio Pineda Flores, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  calculate the imaginary part of the electronic self-energy
!
!  NOTE: We use openmp parallel for qset.
!  NOTE: Segment fault may occur if openmp private space runs out.
!    we allocated memory for large arrays in each threads. however, the 
!    private stack space each thread has is limited (defined by OMP_STACKSIZE)
!    default value is implememention dependent and is usually quite samll. 
!    Segmentation fault may appear if run of of space. We can specify the size 
!    of the private space by 'export OMP_STACKSIZE=xxK (or xxM, xxG)'
!    or use openmp runtime subroutines: kmp_{set, get}_stacksize_s().
!
! Maintenance:
!===============================================================================

module transition_dipole
   use pert_const, only: dp, pi, czero, cone
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
      mp_split_pools, distribute_points
   use pert_param, only: delta_smear, phfreq_cutoff, prefix
   use pert_utils, only: abs2, bose, fermi, gauss, match_table_pulse
   use pert_output,only: progressbar_init, progressbar, set_progress_step, stopwatch
   use polar_correction, only: eph_wan_longrange
   use band_structure, only: electron_wann, solve_eigenvalue_vector, solve_band_velocity_mat
   use phonon_dispersion, only: lattice_ifc, solve_phonon_modes
   use elphon_coupling_matrix, only: elph_mat_wann, eph_fourier_el, eph_fourier_elph, &
      eph_transform_fast, eph_transform_polar
   implicit none

   public :: calc_trans_dipole
contains


subroutine calc_trans_dipole &
      (elph, el, ph, kpts, qset, qwt, tmpr, ef, cond_bmin, bmin, pulse, enk, dipole)
   implicit none
   type(electron_wann), intent(in) :: el
   type(lattice_ifc),   intent(in) :: ph
   type(elph_mat_wann), intent(in) :: elph
   !tmpr: temperature array; ef: chemical potential array; their size is the same.
   !qwt: q-points (qset) weight. kpts: k-points in the current pool. kpts(3,:)
   real(dp), intent(in) :: kpts(:,:), qset(:,:), qwt(:), tmpr, ef(:)
   !compute imsigma for bands start from bmin.
   integer, intent(in) :: cond_bmin, bmin 
   real(dp), intent(in) :: pulse ! energy of the pulse
   !dipole(3, 3, nbnd, nkpt), enk(nband, nkpt)
   real(dp), intent(out) :: enk(:,:)
   complex(dp), intent(out) :: dipole(:,:,:,:) 
   ! local variables
   integer :: mem_max, step_max, nk_pool, qst_pool, qend_pool, nq_pool, nstep_k, ikc
   integer :: step_q, nstep_q, istep_k, kst, kend, step_k, n, icount, iter
   integer :: numk, numq, numb, ik, iq, it, im, ib, jb, i, bmax, ndim(4), neph
   real(dp) :: xk(3), xq(3), xkq(3), wcut, e_thr, ph_n, el_fkq, dt1, dt2, scat
   
   real(dp) :: ek(el%nb), wq(ph%nm), ekq(el%nb)
   complex(dp) :: gpol(ph%nm), ukq(el%nb, el%nb), uq(ph%nm, ph%nm)
   logical,  allocatable :: ltable(:,:,:)
   real(dp), allocatable :: g2(:,:,:)
   complex(dp), allocatable :: uk(:,:,:), gkq(:,:,:), g_kerp(:,:,:,:,:,:)
   complex(dp), allocatable :: vk_mat(:,:,:,:), vkq_mat(:,:,:)
   integer, pointer :: ik_loc(:)
   complex(dp) :: ttrans(3), trans_outer(3,3)
   integer :: il, ix, iy
   real(dp), allocatable:: enk_all(:,:)
   complex(dp), allocatable :: dipole_local(:,:,:,:) 

   call start_clock('calc_trans_dipole')
   enk = 0.0_dp;  dipole = 0.0_dp;  ndim = shape(dipole)

   !#. of selected bands, #.temperature, #. of kpts, number of q-points
   numb = ndim(3);  numk = ndim(4);  numq = size(qset,2)
   bmax = bmin + numb -1
   !array sanity check
   if(size(kpts,2).ne.numk .or. size(enk,1).ne.numb .or. size(enk,2) .ne. numk) &
      call errore('calc_selfel','arguments dimension mismatch.',1)
   if(bmax > elph%nb) &
      call errore('calc_selfel','bmin and imsgm mismatch',1)
   if(elph%nb .ne. el%nb .or. elph%na .ne. ph%na) &
      call errore('calc_selfel','array dimension mismatch',1)

   allocate(enk_all(el%nb,numk))
   enk_all = 0.0_dp

   e_thr = delta_smear*3.0_dp  !exp(-3*3) ~ 10^-4
   wcut  = phfreq_cutoff
   neph = el%nb * el%nb * ph%na * 3
   !max memory per process for g_kerp, hard-coded here, not sure what is the optimal value
   mem_max  = 1024*1024*1024/2 !8GB
   step_max = mem_max / ( neph * elph%max_nrp )
   
   ik_loc => null()
   if(numk > 2*step_max*npool) then
      !interleave distribution of kpoints over pools, for better load balance
      call distribute_points(numk, nk_pool, ik_loc, .true.)
      !call mp_split_pools(numk, kst_pool, kend_pool, nk_pool)
      qst_pool = 1;     qend_pool = numq;    nq_pool = numq
   else
      call distribute_points(numk, nk_pool, ik_loc, .false.)
      !kst_pool = 1;     kend_pool = numk;    nk_pool = numk
      !distrbite q among pools
      call mp_split_pools(numq, qst_pool, qend_pool, nq_pool)
   endif
   !
   call set_progress_step(nk_pool, step_k, nstep_k, step_max)
   if(nstep_k >= 8) then
      !tracking progess on k only
      step_q = nq_pool;    nstep_q = 1
   else
      !track step on k and q combined
      call set_progress_step(nq_pool, step_q, nstep_q)
   endif
   !output debug info
   !write(stdout,'(5x, 3(a,i8,1x))') &
   !   'step_k:', step_k, 'nstep_k:', nstep_k, 'step_max:', step_max
   
   !allocate work sapce
   allocate( uk(el%nb, el%nb, step_k), &
      g_kerp(3, elph%max_nrp, elph%nb, elph%nb, elph%na, step_k) )
   allocate( vk_mat(3, el%nb, el%nb, step_k))

   if(ionode) call progressbar_init('Transition dipole:')
   icount = 0
   do istep_k = 1, nstep_k
      if (ionode) write(*,*) istep_k, nstep_k
      kst = (istep_k-1)*step_k + 1
      kend = min(kst+step_k-1, nk_pool)
      
      uk = czero;  g_kerp = czero
      vk_mat = czero
!$omp parallel do schedule(guided) default(shared) private(n, ikc, ik, xk, ek)
      do n = kst, kend
         ikc = n - kst + 1
         !electronic wavefunction at ik
         ik = ik_loc(n)
         xk = kpts(:, ik)
         call solve_eigenvalue_vector(el, xk, ek, uk(:,:,ikc))
         enk_all(:,ik) = ek
         enk(:,ik) = ek(bmin:bmax)
         call eph_fourier_el(elph, xk, g_kerp(:,:,:,:,:,ikc))
         ! vk mat is on all bands
         call solve_band_velocity_mat(el, xk, vk_mat(:,:,:,ikc)) 
      enddo
!$omp end parallel do

      ! 
!$omp parallel default(shared) private(iq, xq, wq, uq, gpol, i, n, ik, ikc, xkq, &
!$omp& ekq, ukq, ltable, gkq, im, ib, jb, dt1, dt2, it, ph_n, el_fkq, scat, iter, &
!$omp vkq_mat, ttrans, trans_outer, dipole_local) 
      allocate( ltable(el%nb, numb, ph%nm) )
      allocate( gkq(el%nb, el%nb, ph%nm) )
      allocate( vkq_mat(3, el%nb, el%nb))
      allocate( dipole_local(3,3,numb,numk))
      dipole_local = czero
!$omp do schedule(guided)
      do iq = qst_pool, qend_pool
         xq = qset(1:3, iq)
         !get phonon frequcies and eigen-displacement at iq
         call solve_phonon_modes(ph, xq, wq, uq)
         !get e-ph matrix elements in wannier gauge and cartesian coord.
         gpol = czero
         if(elph%lpol) call eph_wan_longrange(elph%pol, xq, gpol, uq)
         
         do n = kst, kend
            ikc = n - kst + 1
            !electronic wavefunction at ik
            ik = ik_loc(n)
            xkq = kpts(:,ik) + xq
            !get electronic wavefunction at ikq
            call solve_eigenvalue_vector(el, xkq, ekq, ukq)
            !check if there is any scattering channel perserve energy conservation
            call match_table_pulse(ekq, enk(:,ik), wq, pulse, wcut, e_thr, ltable)
            if(.not. any(ltable)) cycle
            !
            !!!here
            call eph_fourier_elph(elph, xq, g_kerp(:,:,:,:,:,ikc), gkq)
            !transfor to phonon mode and bloch gauge.
            call eph_transform_fast(elph, uq, uk(:,:,ikc), ukq, gkq, gpol)

            !compute |g|^2
            vkq_mat(:,:,:) = czero
            call solve_band_velocity_mat(el, xkq, vkq_mat(:,:,:))

            do im = 1, ph%nm
            ! compute imsgm for bands between bmin:(bmin+numb-1)
            do  i = 1, numb
            do jb = 1, el%nb
               if(.not. ltable(jb,i,im)) cycle
               ib = i + bmin - 1

               if ((ib >= cond_bmin) .and. (jb < cond_bmin)) then
                  ! conduction band 
                  ! dt1 is emission, correspond to N+1
                  dt1 = gauss( delta_smear, enk(i,ik) - ekq(jb) + wq(im) - pulse )
                  dt2 = gauss( delta_smear, enk(i,ik) - ekq(jb) - wq(im) - pulse )
               else if ((jb >= cond_bmin) .and. (ib < cond_bmin)) then
                  !cycle
                  ! valence band 
                  ! dt1 is emission
                  dt1 = gauss( delta_smear, -enk(i,ik) + ekq(jb) + wq(im) - pulse )
                  dt2 = gauss( delta_smear, -enk(i,ik) + ekq(jb) - wq(im) - pulse )
               else
                  cycle
               endif
               ! compute T, gkq in the order of jb, ib, im
               ttrans(:) = czero
               do il = 1, el%nb 
                  if ((ib >= cond_bmin) .and. (jb < cond_bmin)) then
                     ttrans(:) = ttrans + conjg(gkq(il, ib, im)) * vkq_mat(:,il,jb) / (pulse - (ekq(il) - ekq(jb)))
                     ttrans(:) = ttrans + conjg(gkq(jb, il, im)) * vk_mat(:,ib,il,ikc) / ((enk(i,ik) - enk_all(il,ik)) - pulse)
                  else if ((jb >= cond_bmin) .and. (ib < cond_bmin)) then
                     ! gkq(mk+q, lk, im), v l n k omegalk - omega nk
                     ! gkq(nk, lk+q, im) = gkq*(lk+q, nk, im)
                     ttrans(:) = ttrans + gkq(jb, il, im) * vk_mat(:,il,ib,ikc) / (pulse - (enk_all(il,ik) - enk(i,ik)))
                     ! gkq(lk+q, nk, im), v mlk+q omegamk+q - omega lk+q
                     ttrans(:) = ttrans + gkq(il, ib, im) * vkq_mat(:,jb,il) / ((ekq(jb) - ekq(il)) - pulse)
                  endif
               enddo
               ph_n = bose(tmpr, wq(im)) ! N_q\mu
               scat = ((ph_n+1)*dt1 + ph_n*dt2) * qwt(iq)

               trans_outer = czero
               do ix = 1, 3
                  do iy = 1, 3
                     trans_outer(ix, iy) = ttrans(ix) * conjg(ttrans(iy))
                  enddo
               enddo
            
               !NOTE: the order of jb and ib are important! Tr[g^\dagger g]
               !scat = scat * pi * qwt(iq) * g2(jb, ib, im) / (2.0_dp*wq(im))
!!NB: only one thread can update this shared variable at a time. negligible overhead.
!               do ix = 1, 3
!                  do iy = 1, 3
!!$omp atomic update
!                     dipole(ix, iy, i, ik) = dipole(ix, iy, i, ik) + scat*trans_outer(ix, iy)
!                  enddo
!               enddo
               dipole_local(:, :, i, ik) = dipole_local(:,:, i, ik) + scat*trans_outer(:,:)
            enddo; enddo; enddo
         enddo

         !track progress
!!$omp atomic update
!         icount = icount + 1
!         iter = icount
!         if( mod(iter, step_q*nstep_k).eq.0 .or. iter.eq.(nq_pool*nstep_k) ) then
!!$omp critical (calc_selfel_progress)
!            write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/(nq_pool*nstep_k), '%'
!!$omp end critical (calc_selfel_progress)
!         endif

      enddo
!$omp end do

!$omp critical
      dipole = dipole + dipole_local
!$omp end critical

      deallocate(gkq, ltable)
      deallocate(vkq_mat)
      deallocate(dipole_local)
!$omp end parallel
   enddo
   deallocate(uk, g_kerp)
   deallocate(vk_mat)
   call mp_sum(dipole, inter_pool_comm)
   if(numk .ne. nk_pool) call mp_sum(enk, inter_pool_comm)
   !release space
   if( associated(ik_loc) ) deallocate( ik_loc )

   call stop_clock('calc_trans_dipole')
end subroutine calc_trans_dipole

end module transition_dipole
