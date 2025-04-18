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
!
! Maintenance:
!===============================================================================

!> Numerical calculation of the imaginary part of the off-shell (frequency dependent) electron-self energy in the Fan-Migdal approximation. 
!! Note: this module corresponds to the perturbo calculation mode `spectral-se` in the file `calc_selfenergy.f90`. 
module selfenergy
   use pert_const, only: dp, pi, czero, cone
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
      mp_split_pools, distribute_points
   use pert_param, only: delta_smear, phfreq_cutoff, prefix
   use pert_utils, only: abs2, bose, fermi, gauss
   use pert_output,only: progressbar_init, progressbar, set_progress_step, stopwatch
   use polar_correction, only: eph_wan_longrange
   use band_structure, only: electron_wann, solve_eigenvalue_vector
   use phonon_dispersion, only: lattice_ifc, solve_phonon_modes
   use elphon_coupling_matrix, only: elph_mat_wann, eph_fourier_el, eph_fourier_elph, &
      eph_transform_fast
   implicit none

   public :: selfenergy_fan
contains

!> Calculates the imaginary part of the off-shell electron self-energy in the Fan-Migdal approximation. 
subroutine selfenergy_fan &
      (elph, el, ph, kpts, egrid, qset, qwt, tmpr, ef, bmin, enk, selfe_i)
   implicit none
   type(electron_wann), intent(in) :: el     !< Wannier basis electronic structure for interpolation
   type(lattice_ifc), intent(in)   :: ph     !< Wannier basis inter-atomic force constants for interpolation
   type(elph_mat_wann), intent(in) :: elph   !< Wannier basis electron-phonon coupling for interpolation

   real(dp), intent(in) :: tmpr(:)           !< Array of temperatures at which to calculate self-energy
   real(dp), intent(in) :: ef(:)             !< Array of Fermi energies for each temperature
   !qwt: q-points (qset) weight. kpts: k-points in the current pool. kpts(3,:)
   real(dp), intent(in) :: kpts(:,:)         !< list of k-points at which to calculate self-energy. size (3 x nk) note: this may be
                                             !! only for the current pool, double check this. 
   real(dp), intent(in) :: qset(:,:)         !< list of q-points summed over in self-energy calculation. size (3 x nq)
   real(dp), intent(in) :: qwt(:)            !< Weight of each q-point. size (nq)
   real(dp), intent(in) :: egrid(:)          !< Frequency grid on which to calculate self-energy. size(nestep)
   integer, intent(in) :: bmin               !< Minimum band for which to calculate self-energy
   real(dp), intent(out) :: enk(:,:)         !< Band energy \f$\varepsilon_{nk}\f$. size (nb x nk)
   real(dp), intent(out) :: selfe_i(:,:,:,:) !< Imaginary part of electron self-energy, \f$ \Sigma_{n\mathbf{k}}(\omega,T) \f$. 
                                             !! size (nfreq,ntemp,nb,nk)
   
   ! local variables
   integer :: step_max, nk_pool, nq_pool, qst_pool, qend_pool, n
   integer :: step_q, nstep_q, istep_k, kst, kend, nstep_k, icount, iter
   integer :: numk, numq, ntmpr, numb, ik, iq, it, im, ib, jb, i, bmax, ndim(4)
   integer :: nmode, mem_max, nestep, ikc, iw, step_k, neph, nbd, nat, max_nrp
   real(dp) :: xk(3), xq(3), xkq(3), wcut, ph_n, el_fkq
   
   real(dp) :: ek(el%nb), wq(ph%nm), ekq(el%nb)
   complex(dp) :: gpol(ph%nm), ukq(el%nb, el%nb), uq(ph%nm, ph%nm)
   real(dp), allocatable :: g2(:,:,:), idt1(:), idt2(:), ipart(:)
   complex(dp), allocatable :: uk(:,:,:), gkq(:,:,:), g_kerp(:,:,:,:,:,:)
   integer, pointer :: ik_loc(:)

   call start_clock('selfenergy_fan')
   ndim = shape(selfe_i); nestep = ndim(1);  ntmpr = ndim(2);  numb = ndim(3)
   !#. of selected bands, #.temperature, #. of kpts, number of q-points
   numk = ndim(4); numq = size(qset,2); bmax = bmin + numb -1
   !array sanity check
   if(size(tmpr).ne.ntmpr .or. size(ef).ne.ntmpr .or. size(kpts,2).ne.numk .or. &
      size(egrid) .ne. nestep .or. size(enk,1) .ne. numb .or. size(enk,2) .ne. numk) &
      call errore('selfenergy_fan','arguments dimension mismatch.',1)
   if(bmax > elph%nb) &
      call errore('selfenergy_fan','bmin and spectral mismatch',1)
   if(elph%nb .ne. el%nb .or. elph%na .ne. ph%na) &
      call errore('selfenergy_fan','array dimension mismatch',1)
   
   !init
   enk = 0.0_dp;  selfe_i = 0.0_dp
   wcut  = phfreq_cutoff
   nmode = ph%nm
   nbd = elph%nb
   nat = elph%na
   neph = nat * nbd * nbd * 3
   max_nrp = elph%max_nrp

   !max memory per process for g_kerp, hard-coded here, not sure what is the optimal value
   mem_max  = 1024*1024*1024/4  !4GB
   step_max = mem_max / (neph * max_nrp)
   
   ik_loc => null()
   !if numk is larger, then k-parallel is more effecient
   ! otherwise use q-parallel to achieve better load balance.
   if(numk > 2*step_max*npool) then
      !interleave distribution of kpoints over pools, for better load balance
      call distribute_points(numk, nk_pool, ik_loc, .true.)
      !call mp_split_pools(numk, kst_pool, kend_pool, nk_pool)
      qst_pool = 1;  qend_pool = numq;  nq_pool = numq
   else
      call distribute_points(numk, nk_pool, ik_loc, .false.)
      !kst_pool = 1;     kend_pool = numk;    nk_pool = numk
      !distrbite q among pools
      call mp_split_pools(numq, qst_pool, qend_pool, nq_pool)
   endif
   
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
   allocate( uk(nbd, nbd, step_k), g_kerp(3, max_nrp, nbd, nbd, nat, step_k) )

   if(ionode) call progressbar_init('Selfenergy:')
   icount = 0
   do istep_k = 1, nstep_k
      kst = (istep_k-1)*step_k + 1
      kend = min(kst+step_k-1, nk_pool)
      
      uk = czero;  g_kerp = czero
!$omp parallel do schedule(guided) default(shared) private(n, ikc, ik, xk, ek)
      do n = kst, kend
         ikc = n - kst + 1
         !electronic wavefunction at ik
         ik = ik_loc(n)
         xk = kpts(:, ik)
         call solve_eigenvalue_vector(el, xk, ek, uk(:,:,ikc))
         enk(:,ik) = ek(bmin:bmax)
         call eph_fourier_el(elph, xk, g_kerp(:,:,:,:,:,ikc))
      enddo
!$omp end parallel do
 
      !
!$omp parallel default(shared) private(iq, xq, wq, uq, gpol, i, n, ik, ikc, xkq, ekq, &
!$omp& ukq, gkq, g2, im, ib, jb, it, ph_n, el_fkq, iw, idt1, idt2, ipart, iter)
      allocate( idt1(nestep), idt2(nestep)) 
      allocate( gkq(el%nb, el%nb, nmode), g2(el%nb, el%nb, nmode), ipart(nestep) )
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
            !
            call eph_fourier_elph(elph, xq, g_kerp(:,:,:,:,:,ikc), gkq)
            !transfor to phonon mode and bloch gauge.
            call eph_transform_fast(elph, uq, uk(:,:,ikc), ukq, gkq, gpol)
            !compute |g|^2
            g2 = abs2( gkq )

            do im = 1, nmode
               if( wq(im) < wcut ) cycle
            ! compute imsgm for bands between bmin:(bmin+numb-1)
            do  i = 1, numb
            do jb = 1, el%nb
               ib = i + bmin - 1
               
               idt1(:) = egrid(:) + ( enk(i,ik) - ekq(jb) + wq(im) )
               idt2(:) = egrid(:) + ( enk(i,ik) - ekq(jb) - wq(im) )

               idt1(:) = - pi * gauss( delta_smear, idt1(:) )
               idt2(:) = - pi * gauss( delta_smear, idt2(:) )
               
               do it = 1, ntmpr
                  ph_n = bose(tmpr(it), wq(im)) ! N_q\mu
                  el_fkq = fermi( tmpr(it), ekq(jb)-ef(it) )
                  !rpart(:) = (el_fkq + ph_n)*rdt1(:) + (1.0_dp - el_fkq + ph_n)*rdt2(:)
                  !rpart(:) = rpart(:) * qwt(iq) * g2(jb, ib, im) / (2.0_dp*wq(im))
                  
                  ipart(:) = (el_fkq + ph_n)*idt1(:) + (1.0_dp - el_fkq + ph_n)*idt2(:)
                  !NOTE: the order of jb and ib are important! Tr[g^\dagger g]
                  ipart(:) = ipart(:) * ( qwt(iq) * g2(jb, ib, im) / (2.0_dp*wq(im)) )

!!NB: only one thread can update this shared variable at a time. negligible overhead.
                  do iw = 1, nestep
!$omp atomic update
                     selfe_i(iw,it,i,ik) = selfe_i(iw,it,i,ik)  + ipart(iw)
                  enddo
               enddo
            enddo; enddo; enddo
         enddo
            
         !track progress
!$omp atomic update
         icount = icount + 1
         iter = icount
         if( mod(iter, step_q*nstep_k).eq.0 .or. iter.eq.(nq_pool*nstep_k) ) then
!$omp critical (selfenergy_progress)
            write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/(nq_pool*nstep_k), '%'
!$omp end critical (selfenergy_progress)
         endif

      enddo
!$omp end do
      deallocate(gkq, g2, idt1, idt2, ipart)
!$omp end parallel
   enddo
   deallocate(uk, g_kerp)
   call mp_sum(selfe_i, inter_pool_comm)
   if(nk_pool .ne. numk)  call mp_sum(enk, inter_pool_comm)
   if( associated(ik_loc) ) deallocate( ik_loc )

   call stop_clock('selfenergy_fan')
end subroutine selfenergy_fan

end module selfenergy
