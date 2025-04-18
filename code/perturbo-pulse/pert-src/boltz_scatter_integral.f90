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
!  solve iterative linearized Boltzmann equation, and real-time dynamics
!
! Maintenance:
!===============================================================================

#include "cdyn_config.h"

#undef _OPENACC

module boltz_scatter_integral
   use iso_c_binding
   use omp_lib
   use kinds, only: i8b
   use pert_const, only: dp, pi
   use boltz_grid, only: grid, symcheck_absvbs_ongrid, symmetrize_absvbs_ongrid
   use boltz_utils, only: idx2band
   use pert_utils, only: bose, fermi, gauss
   use pert_param, only: delta_smear, delta_smear_ph, boltz_qdim, &
      adapt_smear_eph, adapt_smear_phph, phfreq_cutoff, phfreq_cutoff_ph, &
      symm
   use pert_output, only: stopwatch
   !syp - [TOCHECK]
   use boltz_scatter, only: num_kq_pairs, num_qq_pair, num_scatter_channels, scatter, scatter_channels, qgrid, &
      scatter_ph, ph3_num_scat, tscat_ph, tsmear_list, tsmear_ph_list
   use qe_mpi_mod, only: stdout, ionode, mp_barrier, mp_sum, inter_pool_comm, my_pool_id
   implicit none
   private

   public :: rates_scat_int, trans_scat_int, cdyna_scat_int
   public :: cphdyna_scat_int, ph3_scat_int
   public :: ph3_scat_int_interp2
   public :: ph3_scat_int_interp_half
   !syp - added
   public :: rates_phe_scat_int, rates_ph3_scat_int, &
      trans_ph3_scat_int, drag_scat_4e, drag_scat_4ph
contains

! compute scattering rate
subroutine rates_scat_int(kg, eptemp, efermi, rates)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, efermi
   real(dp), intent(out) :: rates(:,:) !rates(kg%numb, kg%numk)
   !local variables
   integer :: i, iq, ik, ikq, im, n, m, bands(3)
   integer(kind=i8b) :: ns
   real(dp) :: wq, enk, emkq, ph_n, fkq, nk2mkq, mkq2nk
#if !defined(STORE_G2DT)
   real(dp) :: g2, dt1, dt2
#else
   real(dp) :: g2_dt1, g2_dt2
#endif

   call stopwatch('rates_scat_int', 'start')
   rates(:,:) = 0.0_dp

#if !defined(STORE_G2DT)
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns, &
!$omp& im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, mkq2nk, nk2mkq, bands)
#else
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns, &
!$omp& im, n, m, wq, enk, emkq, g2_dt1, g2_dt2, ph_n, fkq, mkq2nk, nk2mkq, bands)
#endif
#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         !write three logical variables, whether the pointer variables has been assigned
         wq   = qgrid%freq(im, iq)
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)

         ! phonon occupation
         ! the constrain of 'wq > eps_acustic' is aready done in scat_setup.
         ph_n = bose(eptemp, wq) ! N_q\mu

#if !defined(STORE_G2DT)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2

         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq)
         if ( adapt_smear_eph ) then
             dt1 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(ns), enk - emkq + wq)
             ! delta(Enk - Emkq -wq)
             dt2 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(ns), enk - emkq - wq)
         else
             dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
             ! delta(Enk - Emkq -wq)
             dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
         endif
#else
         g2_dt1 = scatter_channels(ns)%g2_dt1
         g2_dt2 = scatter_channels(ns)%g2_dt2
#endif

         ! compute scattering rates
         fkq = fermi(eptemp, emkq-efermi)
#if !defined(STORE_G2DT)
         mkq2nk = g2 * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
#else
         mkq2nk = (g2_dt1 * (ph_n + fkq) + g2_dt2 * (ph_n + 1.0_dp - fkq))
#endif

! lock memory when updating this variable
!syp - why not just use reduction(+:rates) here? more efficient.
!$omp atomic
         rates(n, ik)  = rates(n, ik)  + mkq2nk

         !update mk+q with contribution from k
         fkq = fermi(eptemp, enk-efermi)
         !delta(Ekq-Ek+wq)=delta(Ek-Ekq-wq); delta(Ekq-Ek-wq)=delta(Ek-Ekq+wq)
#if !defined(STORE_G2DT)
         nk2mkq = g2 * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
#else
         nk2mkq = (g2_dt2 * (ph_n + fkq) + g2_dt1 * (ph_n + 1.0_dp - fkq))
#endif

!$omp atomic
         rates(m, ikq) = rates(m, ikq) + nk2mkq
      enddo  !; enddo; enddo
#if defined(SCAT_FWD)
   enddo ! i loop
#endif
!$omp end parallel do
   ! combine contributions from all the pools/processors.
   call mp_sum(rates, inter_pool_comm)
   call stopwatch('rates_scat_int', 'stop')
end subroutine rates_scat_int

subroutine rates_phe_scat_int(kg, eptemp, efermi, qg, pecol)
   use boltz_utils, only: kpt2num, num2kpt, kpts_minus, num2ipt
   implicit none
   type(grid), intent(in) :: kg, qg
   real(dp), intent(in) :: eptemp, efermi
   ! e-ph collision term computed from given distribution function dist(:,:)
   real(dp), intent(out) :: pecol(qg%numb, qg%nk)
   !local variables
   integer :: i, iq, inq, ik, ikq, im, n, m, bands(3), iq_test
   integer(kind=i8b) :: ns
   real(dp):: wq, enk, emkq, ph_n, ph_nn, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq, imph, imphn
   real(dp) :: testq(3), wnq
   integer :: iqnums(3)
   logical :: symmornot
   real(dp) :: weight_ratio, smear_tmp, ediff_abs, ediff_ems 
   real(dp) :: factor, dscat

!#if !defined(STORE_G2DT)
   real(dp) :: g2, dt1, dt2
!#else
   real(dp) :: g2_dt1, g2_dt2
!#endif
   call start_clock('rates_phe_scat_int')
   ! kelly - using phonon_grid qgrid is a terrible idea here. too confusing
   ! kelly's conclusion is that the formula is correct
   ! we are doing half of k k' range, so twice the amount goes to each q
   ! only need to copy the occupation of q to that of -q at the end
   ! which can be done outside the loop
   ! only works when doing half of k k' range!!
   ! this should be effectively the same as updating -q in the half-range loop
   ! with some extra factor in place - worked out the factor to be 2
   ! we here are importing the whole BZ for qgrids and using type grid, qg 

   pecol(:,:) = 0.0E0_dp
   weight_ratio = kg%kweight/qg%kweight
   factor = qg%kweight * pi * 4.0_dp
    !!!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, &
    !!!$omp&  ns, im, n, m, wq, bands, enk, emkq, g2, dt1, dt2, ph_n, ph_nn, fnk, fmkq, &
    !!!$omp&  fabs, fem, mkq2nk, nk2mkq, imph, imphn, inq, wnq) 
!#if defined(_OPENACC)
!!OpenACC acceleration
!!$acc kernels
!#else
! OpenMP acceleration
!$omp parallel do schedule(guided) &
!#if defined(CDYN_USE_ARRAY_REDUCE)
!!$omp&  reduction(+:epcol) &
!#endif
!$omp&  default(shared) private(i, iq, ik, ikq, ns, bands, im, n, m, wq, enk, emkq, &
!#if !defined(STORE_G2DT)
!$omp&  g2, dt1, dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq, &
!$omp&  wnq, ph_nn, imph, imphn, smear_tmp, ediff_abs, ediff_ems, dscat, iqnums, inq)
!#else
!!$omp&  g2_dt1, g2_dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq)
!#endif
!#endif
#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
      iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      ! loop over all scatter-channels, using the scat_index to lookup
      ! the corresponding (iq, ik, ikq) values
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
         iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif
         !sypeng: unsure: maybe need to be modified for gpu 
         inq = kpt2num(-1.0_dp*num2kpt(iq-1, qg%ndim), qg%ndim) + 1
         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         wq   = qg%enk(im, iq)
         wnq  = qg%enk(im, inq)
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)

         fnk  = fermi(eptemp, enk-efermi) 
         fmkq = fermi(eptemp, emkq-efermi) 
         !if ( abs(enk - emkq) < phfreq_cutoff ) cycle
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
        !g2 = scat(i)%eph_g2(ns)
         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! syp - 
         !       1. here is 4 pi, not 2 pi as Eq. (34) in Bernardi et al. 2016
         !       2. we use q process: k + q -> k' and -q process: k' + (-q) -> k
         !          to tranverse the whole (k, k') space, because in set_kq_pairs
         !          we just consider half number of k' for given k (for sure,the 
         !          k is for full BZ there).
         if ( adapt_smear_eph ) then
            smear_tmp = tsmear_list(i)
         else
            smear_tmp = delta_smear
         endif
         dscat = factor * scatter_channels(ns)%eph_g2
         dt1 = gauss(smear_tmp, enk - emkq + wq)
         ! delta(Enk - Emkq -wq) : for -q : k' + (-q) -> k
         dt2 = gauss(smear_tmp, enk - emkq - wnq)
         fabs = fnk - fmkq
         fem  = fmkq - fnk
         imph  =  dscat * dt1 * fabs * weight_ratio !(w0g1*fabs - w0g2*fem)
         imphn =  dscat * dt2 * fem  * weight_ratio !(w0g1*fabs - w0g2*fem)
        !ediff_abs = abs(enk - emkq + wq)
        !ediff_ems = abs(enk - emkq - wnq)
        !dscat = factor * scatter_channels(ns)%eph_g2
        !if (ediff_abs < 3.0_dp*smear_tmp) then 
        !   dt1 = gauss(smear_tmp, ediff_abs)
        !   fabs = fnk - fmkq
        !   imph = dscat* dt1 * fabs * weight_ratio
        !else
        !   imph = 0.0_dp
        !endif
        !if (ediff_ems < 3.0_dp*smear_tmp) then 
        !   dt2 = gauss(smear_tmp, ediff_ems)
        !   fem  = fmkq - fnk
        !   imphn = dscat* dt2 * fem * weight_ratio
        !else
        !   imphn = 0.0_dp
        !endif

!$omp atomic
         pecol(im, iq) = pecol(im, iq) + imph
         !
         iqnums =  num2ipt(qg%kpt(iq), qg%ndim)
         if ( product(iqnums - boltz_qdim/2) .eq. 0 ) then
            cycle
         else
!$omp atomic
            pecol(im, inq) = pecol(im, inq) + imphn
         endif
      enddo   !; enddo; enddo ! ns loop
#if defined(SCAT_FWD)
   enddo ! i loop
#endif
!$omp end parallel do
   call mp_sum(pecol, inter_pool_comm)

  !!syp - symmetry of equivalent q points
  !call symcheck_absvbs_ongrid(qg, pecol, symmornot)
  !if (.not. symmornot) then
  !   if (ionode) write(stdout, *) '*** will symmetrize ph-e scattering rate'
  !   call symmetrize_absvbs_ongrid(qg, pecol)
  !endif

   call stop_clock('rates_phe_scat_int')
end subroutine rates_phe_scat_int

subroutine rates_ph3_scat_int(qg, eptemp, ph3col)
   use boltz_utils, only: kpts_plus_invert
   implicit none
   type(grid), intent(in) :: qg
   real(dp), intent(in) :: eptemp
   real(dp), intent(inout) :: ph3col(qg%numb, qg%nk)

   !local variables
   integer(kind=i8b) :: is
   integer :: ivec, iq0, iq1, iq2, it, im0, im1, im2, bands(3)
   real(dp) :: wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat
   real(dp) :: occ_sum_0, occ_sum_1, occ_sum_2, factor
   real(dp) :: qweight, e_thr
   integer :: scat_index
   real(dp), allocatable :: ph3col_tmp(:,:)

   call start_clock('rates_ph3_scat_int')
   qweight = 1.0_dp / (boltz_qdim(1) * boltz_qdim(2) * boltz_qdim(3))
   e_thr = delta_smear_ph*3.0_dp

   if(symm)then
      allocate(ph3col_tmp(qg%numb, qg%nk_irr))
   else
      allocate(ph3col_tmp(qg%numb, qg%nk))
   endif

   !syp - From Marco's review paper, the factor is 36 by a correction of Eq. (32) which is off by a factor of 2
   !    - Here I use a different counting prcess as Kelly, I will adopt irreducible points for q_0
   !    - But will loop over full Brillouin zone for q_1 and q_2 correspondingly. 
   !    - For scattering process, I will respect the Eq. (32) including both emission and absorption process for 
   !    - each q_0 which is different with Kelly. So the finally factor is 36, not 72 in Kelly's code. 
   !    - (Though kelly code is also 36 here, but she times another "2" elsewhere to make the final factor 72.)
   factor = 9.0E0_dp * 4.0E0_dp * pi * qweight 

   ph3col_tmp(:,:) = 0.0E0_dp
!$omp parallel do schedule(guided) default(shared) private(is, iq0, iq1, iq2, bands, im0, im1, im2, &
!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, scat_index,  dscat) &
!$omp& reduction(+:ph3col_tmp)
   do is = 1, ph3_num_scat
      ! key step, to keep g2 the same 
      g2 = tscat_ph(is)%g2
      ! ignore small values
      if ( g2 < 1.0E-30_dp ) then
         cycle
      endif
      scat_index = tscat_ph(is)%scat_index
      iq0 = scatter_ph(1, scat_index)
      iq1 = scatter_ph(2, scat_index)
      iq2 = kpts_plus_invert(qg%kpt(iq0), qg%kpt(iq1), qg%ndim) + 1
      ! band index is defined the same, here we use qg, but qg_fine%numb is the same 
      bands = idx2band( tscat_ph(is)%bands_index, qg%numb )
      im2 = bands(1); im0 = bands(2); im1 = bands(3)
         
      ! eigenvalue
      wq0 = qg%enk(im0, iq0); wq1 = qg%enk(im1, iq1); wq2 = qg%enk(im2, iq2)
        

      !--------------Someone's comments, I keep it here for reference----------------
      ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
      !> d3q  add the factor of 18pi/hbar^2/36, 0.5^3 is included in g2, hbar is omitted in atomic units.
      !> tdep add prefactor of hbar*pi/2, 1/8 is preincluded in g2
      !------------------------------------------------------------------------------
         
      ! phonon occupation
      ocq0 = bose(eptemp, wq0); ocq1 = bose(eptemp, wq1); ocq2 = bose(eptemp, wq2)
      
      !syp - I will split processes into absorption and emission process
      !    - this is a little different with electron strategy, we do it strickly
      !    - for each i_scatt, either absorption or emission, not both.
      ! 1. absorption 
      if (abs(wq0 + wq1 - wq2) < e_thr)then  
         if ( adapt_smear_phph ) then
            e_delta = gauss(tsmear_ph_list(is), wq1 + wq0 - wq2)
         else
            e_delta = gauss(delta_smear_ph, wq1 + wq0 - wq2)
         endif

         ! syp - [TOCHECK] @ 20241006: I don't know whether we need to keep 
         !       this here, which happens in ph3_scat_int by Kelly
         !! halving equal ik and iq
         !if ((im0 .eq. im1) .and. (iq0 .eq. iq1)) then
         !   e_delta = e_delta/2.0_dp
         !endif

         !syp - the factor of 2 is for two absorption processes
         dscat = g2 * e_delta * factor * 2.0_dp * (ocq1 - ocq2)
      
      ! 3. emission 
      elseif(abs(wq0 - wq1 - wq2) < e_thr) then  
         if ( adapt_smear_phph ) then
            e_delta = gauss(tsmear_ph_list(is), wq0 - wq1 - wq2)
         else
            e_delta = gauss(delta_smear_ph, wq0 - wq1 - wq2)
         endif
         ! syp - [TOCHECK] @ 20241006: I don't know whether we need to keep 
         !       this here, which happens in ph3_scat_int by Kelly
         !! halving equal ik and iq
         !if ((im0 .eq. im1) .and. (iq0 .eq. iq1)) then
         !   e_delta = e_delta/2.0_dp
         !endif
         dscat = g2 * e_delta * factor * (1.d0 + ocq1 + ocq2)
      endif
      ph3col_tmp(im0,qg%kpt2ir(iq0)) = ph3col_tmp(im0,qg%kpt2ir(iq0)) + dscat
   enddo ! is loop
!$omp end parallel do

   call mp_sum(ph3col_tmp, inter_pool_comm)

   !map irreducible points to full BZ data
   if(symm)then
       do iq0 = 1, qg%nk
          ph3col(:, iq0) = ph3col_tmp(:, qg%kpt2ir(iq0))
       enddo
   else
      ph3col = ph3col_tmp
   endif

   !call mp_barrier(inter_pool_comm)
   call stop_clock('rates_ph3_scat_int')
end subroutine rates_ph3_scat_int

subroutine trans_ph3_scat_int(qg, eptemp, mfd, ph3col)
   use pert_data, only: symop_inv_cart
   use boltz_utils, only: kpts_plus_invert
   implicit none
   type(grid), intent(in) :: qg
   real(dp), intent(in) :: eptemp
   real(dp), intent(in) :: mfd(:,:,:) !(3 or 6, nband, nk)
   real(dp), intent(out) :: ph3col(:,:,:)!(3 or 6, qg%numb, qg%nk)

   !local variables
   integer(kind=i8b) :: is
   integer :: iq0, iq1, iq2, it, im0, im1, im2, bands(3), icomp, ncomp
   real(dp) :: wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat
   real(dp) :: factor, qweight, e_thr
   integer :: scat_index
   real(dp), allocatable :: ph3col_tmp(:,:,:)

   call start_clock('trans_ph3_scat_int')

   !sanity check
   ncomp = size(mfd, 1)
   if( ncomp .ne. size(ph3col, 1) ) &
      call errore('trans_ph3_scat_int',"mismatch dimension in input arguments!",1)

   qweight = 1.0_dp / (boltz_qdim(1) * boltz_qdim(2) * boltz_qdim(3))
   e_thr = delta_smear_ph*3.0_dp

   !syp - From Marco's review paper, the factor is 36 by a correction of Eq. (32) which is off by a factor of 2
   !    - Here I use a different counting prcess as Kelly, I will adopt irreducible points for q_0
   !    - But will loop over full Brillouin zone for q_1 and q_2 correspondingly. 
   !    - For scattering process, I will respect the Eq. (32) including both emission and absorption process for 
   !    - each q_0 which is different with Kelly. So the finally factor is 36, not 72 in Kelly's code. 
   !    - (Though kelly code is also 36 here, but she times another "2" elsewhere to make the final factor 72.)
   factor = 9.0E0_dp * 4.0E0_dp * pi * qweight 
   if(symm)then
       allocate(ph3col_tmp(ncomp, qg%numb, qg%nk_irr))
   else
       allocate(ph3col_tmp(ncomp, qg%numb, qg%nk))
   endif
   ph3col(:,:,:) = 0.0E0_dp
   ph3col_tmp(:,:,:) = 0.0E0_dp 
!$omp parallel do schedule(guided) default(shared) private(is, iq0, iq1, iq2, bands, im0, im1, im2, &
!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat, scat_index,  icomp) &
!$omp& reduction(+:ph3col_tmp)
   do is = 1, ph3_num_scat
      g2 = tscat_ph(is)%g2
      ! ignore small values
      if ( g2 < 1.0E-30_dp ) then
         cycle
      endif
      scat_index = tscat_ph(is)%scat_index
      iq0 = scatter_ph(1, scat_index)
      iq1 = scatter_ph(2, scat_index)
      iq2 = kpts_plus_invert(qg%kpt(iq0), qg%kpt(iq1), qg%ndim) + 1
      ! band index is defined the same, here we use qg, but qg_fine%numb is the same 
      bands = idx2band( tscat_ph(is)%bands_index, qg%numb )
      im2 = bands(1); im0 = bands(2); im1 = bands(3)
         
      ! eigenvalue
      wq0 = qg%enk(im0, iq0); wq1 = qg%enk(im1, iq1); wq2 = qg%enk(im2, iq2)
         
      ! phonon occupation
      ocq0 = bose(eptemp, wq0); ocq1 = bose(eptemp, wq1); ocq2 = bose(eptemp, wq2)

      !syp - I will split processes into absorption and emission process
      !    - this is a little different with electron strategy, we do it strickly
      !    - for each i_scatt, either absorption or emission, not both.
      ! 1. absorption 
      if (abs(wq0 + wq1 - wq2) < e_thr)then  
         if ( adapt_smear_phph ) then
            e_delta = gauss(tsmear_ph_list(is), wq1 + wq0 - wq2)
         else
            e_delta = gauss(delta_smear_ph, wq1 + wq0 - wq2)
         endif
         dscat = g2 * e_delta * factor * 2.0_dp * (ocq1 - ocq2)
         do icomp = 1, ncomp 
            ph3col_tmp(icomp, im0, qg%kpt2ir(iq0)) = ph3col_tmp(icomp, im0, qg%kpt2ir(iq0)) + dscat * &
               (mfd(icomp, im2, iq2) - mfd(icomp, im1, iq1))
         enddo

      ! 2. emission 
      elseif(abs(wq0 - wq1 - wq2) < e_thr) then  
         if ( adapt_smear_phph ) then
            e_delta = gauss(tsmear_ph_list(is), wq0 - wq1 - wq2)
         else
            e_delta = gauss(delta_smear_ph, wq0 - wq1 - wq2)
         endif
         dscat = g2 * e_delta * factor*(1.d0 + ocq1 + ocq2)
         do icomp = 1, ncomp 
            ph3col_tmp(icomp, im0, qg%kpt2ir(iq0)) = ph3col_tmp(icomp, im0, qg%kpt2ir(iq0)) + dscat * &
               (mfd(icomp, im2, iq2) + mfd(icomp, im1, iq1))
         enddo

      endif
   enddo ! is loop
!$omp end parallel do

   call mp_sum(ph3col_tmp, inter_pool_comm)

   if(symm)then
      !map irreducible points to full BZ data
      !be careful, the first dimension (Cartesian coordinates) should be rorated accordingly
      do iq0 = 1, qg%nk
         do it = 1, qg%numb
             if (ncomp == 3) then
                ph3col(:, it, iq0) = &
                   matmul(ph3col_tmp(:, it, qg%kpt2ir(iq0)), symop_inv_cart(:,:,qg%symopindex(iq0)) )
             elseif (ncomp == 6) then
                ph3col(1:3, it, iq0) = &
                   matmul(ph3col_tmp(1:3, it, qg%kpt2ir(iq0)), symop_inv_cart(:,:,qg%symopindex(iq0)) )
                ph3col(4:6, it, iq0) = &
                   matmul(ph3col_tmp(4:6, it, qg%kpt2ir(iq0)), symop_inv_cart(:,:,qg%symopindex(iq0)) )
             else
                call errore('trans_ph3_scat_int',"ncomp value does not support, only 3 or 6!",1)
             endif
         enddo
      enddo
   else
      ph3col = ph3col_tmp
   endif

   call stop_clock('trans_ph3_scat_int')
end subroutine trans_ph3_scat_int

! iterative process for tranport calculation
subroutine drag_scat_4ph(kg, qg, eptemp, efermi, mfd_F, drag_4ph)
   use boltz_utils, only: kpt2num, num2kpt, kpts_minus, num2ipt
   implicit none
   type(grid), intent(in) :: kg, qg
   real(dp), intent(in) :: efermi, eptemp, mfd_F(:,:,:) !mfd(6, kg%numb, kg%nk)
   real(dp), intent(out) :: drag_4ph(:,:,:)  !epint(6, qg%numb, qg%nk)
   !local variables
   integer :: iqnums(3)
   integer(kind=i8b) :: ns
   integer :: i, iq, ik, ikq, im, n, m, ic, bands(3), ncomp, inq 
   real(dp) :: wq, wnq, enk, emkq, g2, dt1, dt2, ph_n, ph_nn, nk2mkq, mkq2nk, fnk, fmkq
   
   call stopwatch('drag_scat_4ph','start')

   ncomp = size(mfd_F, 1)

   !sanity check
   if( size(mfd_F, 1) .ne. 6 .or. &
      size(mfd_F, 2) .ne. kg%numb .or. &
      size(mfd_F, 3) .ne. kg%nk ) &
      call errore('drag_scat_ph',"mismatch dimension in input arguments!",1)
   
   if(size(drag_4ph, 1) .ne. 6 .or. &
      size(drag_4ph, 2) .ne. qg%numb .or. &
      size(drag_4ph, 3) .ne. qg%nk ) &
      call errore('drag_scat_ph',"mismatch dimension in input arguments!",1)

   drag_4ph(:,:,:) = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns, inq, &
!$omp&  im, n, m, wq, wnq, enk, emkq, g2, dt1, dt2, ph_n, ph_nn, fnk, fmkq, mkq2nk, nk2mkq, ic, bands, iqnums)
#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
      iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      ! loop over all scatter-channels, using the scat_index to lookup
      ! the corresponding (iq, ik, ikq) values
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
         iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif
         !sypeng: unsure: maybe need to be modified for gpu 
         inq = kpt2num(-1.0_dp*num2kpt(iq-1, qg%ndim), qg%ndim) + 1
         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! eigenvalue
         wq   = qg%enk(im, iq)
         wnq  = qg%enk(im, inq)

         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2
         
         if ( adapt_smear_eph ) then
            dt1 = kg%kweight * pi * 4.0_dp * gauss(tsmear_list(i), enk - emkq + wq)
            dt2 = kg%kweight * pi * 4.0_dp * gauss(tsmear_list(i), enk - emkq - wnq)
         else
            ! delta(Enk - Emkq + wq)
            dt1 = kg%kweight * pi * 4.0_dp * gauss(delta_smear, enk - emkq + wq)
            ! delta(Enk - Emkq -wq)
            dt2 = kg%kweight * pi * 4.0_dp * gauss(delta_smear, enk - emkq - wnq)
         endif

         ! phonon occupation
         ph_n = bose(eptemp, wq) ! N_q\mu
         ph_nn = bose(eptemp, wnq) ! N_q\mu

         ! Fermi function
         fnk = fermi(eptemp, enk-efermi)
         fmkq = fermi(eptemp, emkq-efermi) 

         do ic = 1, ncomp
            !update nk with contribution from k+q
            mkq2nk = g2 * dt1 * (fnk-fmkq) * (mfd_F(ic,n,ik) - mfd_F(ic,m,ikq))
!$omp atomic
            drag_4ph(ic, im, iq)  = drag_4ph(ic, im, iq)  - mkq2nk
            
            iqnums =  num2ipt(qg%kpt(iq), qg%ndim)
            if ( product(iqnums - boltz_qdim/2) .eq. 0 ) then
               cycle
            else
               !syp - could be optimized to nk2mkq = mkq2nk / dt1 * dt2
               nk2mkq = g2 * dt2 * (fmkq-fnk) * (mfd_F(ic,m,ikq) - mfd_F(ic,n,ik))
!$omp atomic
               drag_4ph(ic, im, inq) = drag_4ph(ic, im, inq) - nk2mkq
            endif

         enddo
      enddo
#if defined(SCAT_FWD)
   enddo
#endif
!$omp end parallel do
   !combine contributions from all the pools/processors.
   call mp_sum(drag_4ph, inter_pool_comm)
   call stopwatch('drag_scat_4ph','stop')
end subroutine drag_scat_4ph

! iterative process for tranport calculation
subroutine drag_scat_4e(kg, qg, eptemp, efermi, mfd_N, drag_4e)
   use boltz_utils, only: kpt2num, num2kpt, kpts_minus
   implicit none
   type(grid), intent(in) :: kg, qg
   real(dp), intent(in) :: efermi, eptemp, mfd_N(:,:,:) !mfd(6, qg%numb, qg%nk)
   real(dp), intent(out) :: drag_4e(:,:,:)  !epint(6, kg%numb, kg%nk)
   !local variables
   integer(kind=i8b) :: ns
   integer :: i, iq, ik, ikq, im, n, m, ic, bands(3), ncomp, inq 
   real(dp) :: wq, wnq, enk, emkq, g2, dt1, dt2, ph_n, ph_nn, nk2mkq, mkq2nk, fkq
   
   call stopwatch('drag_scat_4e','start')

   ncomp = size(mfd_N, 1)

   !sanity check
   if( size(mfd_N, 1) .ne. 6 .or. &
      size(mfd_N, 2) .ne. qg%numb .or. &
      size(mfd_N, 3) .ne. qg%nk ) then
      write(stdout, *) 'mfd_N.shape:', size(mfd_N, 1), size(mfd_N, 2), size(mfd_N, 3)
      write(stdout, *) 'qg%numb, qg%nk:', qg%numb, qg%nk
      call mp_barrier(inter_pool_comm)
      call errore('drag_scat_4e',"mfd_N mismatch dimension in input arguments!",1)
   endif
   
   if(size(drag_4e, 1) .ne. 6 .or. &
      size(drag_4e, 2) .ne. kg%numb .or. &
      size(drag_4e, 3) .ne. kg%nk ) then
      write(stdout, *) 'drag_4e.shape:', size(drag_4e, 1), size(drag_4e, 2), size(drag_4e, 3)
      write(stdout, *) 'kg%numb, kg%nk:', kg%numb, kg%nk
      call mp_barrier(inter_pool_comm)
      call errore('drag_scat_4e',"drag_4e mismatch dimension in input arguments!",1)
   endif

   drag_4e(:,:,:) = 0.0_dp
   !here we use openmp parallel do, use shared by default, and use private for private variables
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns, inq, &
!$omp&  im, n, m, wq, wnq, enk, emkq, g2, dt1, dt2, ph_n, ph_nn, fkq, mkq2nk, nk2mkq, ic, bands)
#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
      iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      ! loop over all scatter-channels, using the scat_index to lookup
      ! the corresponding (iq, ik, ikq) values
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
         iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif
         !sypeng: unsure: maybe need to be modified for gpu 
         inq = kpt2num(-1.0_dp*num2kpt(iq-1, qg%ndim), qg%ndim) + 1
         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! eigenvalue
         wq   = qg%enk(im, iq)
         wnq  = qg%enk(im, inq)

         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2
         
         if ( adapt_smear_eph ) then
            dt1 = qg%kweight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq + wq)
            dt2 = qg%kweight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq - wnq)
         else
            ! delta(Enk - Emkq + wq) 
            dt1 = qg%kweight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
            ! delta(Enk - Emkq -wq)
            dt2 = qg%kweight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wnq)
         endif

         ! phonon occupation
         ph_n = bose(eptemp, wq) ! N_q\mu
         ph_nn = bose(eptemp, wnq) ! N_q\mu

         do ic = 1, ncomp
            !update nk with contribution from k+q
            fkq = fermi(eptemp, emkq-efermi)
           !if(abs(fkq)>1E-10)write(stdout, *) 'drag_4e: fkq', fkq, 'emkq', emkq
            mkq2nk = g2 * mfd_N(ic, im, iq) * ( - dt1*(ph_n + fkq) + dt2*(ph_nn + 1.0_dp - fkq))
!$omp atomic
            drag_4e(ic, n, ik)  = drag_4e(ic, n, ik)  + mkq2nk

            ! update mk+q with contribution from k
            fkq = fermi(eptemp, enk-efermi)
            nk2mkq = g2 * mfd_N(ic, im, inq) * ( - dt2*(ph_nn + fkq) + dt1*(ph_n + 1.0_dp - fkq))
!$omp atomic
            drag_4e(ic, m, ikq) = drag_4e(ic, m, ikq) + nk2mkq

         enddo
      enddo
#if defined(SCAT_FWD)
   enddo
#endif
!$omp end parallel do
   !combine contributions from all the pools/processors.
   call mp_sum(drag_4e, inter_pool_comm)
   call stopwatch('drag_scat_4e','stop')
end subroutine drag_scat_4e

! iterative process for tranport calculation
subroutine trans_scat_int(kg, eptemp, efermi, mfd, epint)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: efermi, eptemp, mfd(:,:,:) !mfd(3 or 6, kg%numb, kg%nk)
   real(dp), intent(out) :: epint(:,:,:)  !epint(3 or 6, kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, im, n, m, ic, bands(3), ncomp
   integer(kind=i8b) :: ns
   real(dp) :: wq, enk, emkq, ph_n, nk2mkq, mkq2nk, fkq
#if !defined(STORE_G2DT)
   real(dp) :: g2, dt1, dt2
#else
   real(dp) :: g2_dt1, g2_dt2
#endif

   call stopwatch('trans_scat_int','start')
   !sanity check
   ncomp = size(mfd, 1)
   if( ncomp .ne. size(epint, 1) ) &
      call errore('trans_scat_int',"mismatch dimension in input arguments!",1)

   epint(:,:,:) = 0.0_dp

#if !defined(STORE_G2DT)
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns,&
!$omp&  im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, mkq2nk, nk2mkq, ic, bands)
#else
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns,&
!$omp&  im, n, m, wq, enk, emkq, g2_dt1, g2_dt2, ph_n, fkq, mkq2nk, nk2mkq, ic, bands)
#endif
#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! eigenvalue
         wq   = qgrid%freq(im, iq)
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)

         ! phonon occupation
         ! the constrain of 'wq > eps_acustic' is aready done in scat_setup.
         ph_n = bose(eptemp, wq) ! N_q\mu

#if !defined(STORE_G2DT)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2

         if ( adapt_smear_eph ) then
            ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
            ! delta(Enk - Emkq + wq)
            dt1 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(ns), enk - emkq + wq)
            ! delta(Enk - Emkq -wq)
            dt2 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(ns), enk - emkq - wq)
         else
            ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
            ! delta(Enk - Emkq + wq)
            dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
            ! delta(Enk - Emkq -wq)
            dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
         endif
#else
         g2_dt1 = scatter_channels(ns)%g2_dt1
         g2_dt2 = scatter_channels(ns)%g2_dt2
#endif

         do ic = 1, ncomp
            !update nk with contribution from k+q
            fkq = fermi(eptemp, emkq - efermi)
#if !defined(STORE_G2DT)
            mkq2nk = g2 * mfd(ic, m, ikq) * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
#else
            mkq2nk = mfd(ic, m, ikq) * (g2_dt1 * (ph_n + fkq) + g2_dt2 * (ph_n + 1.0_dp - fkq))
#endif

!$omp atomic
            epint(ic, n, ik)  = epint(ic, n, ik)  + mkq2nk

            ! update mk+q with contribution from k
            fkq = fermi(eptemp, enk-efermi);
#if !defined(STORE_G2DT)
            nk2mkq = g2 * mfd(ic, n, ik) * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
#else
            nk2mkq = mfd(ic, n, ik) * (g2_dt2 * (ph_n + fkq) + g2_dt1 * (ph_n + 1.0_dp - fkq))
#endif

!$omp atomic
            epint(ic, m, ikq) = epint(ic, m, ikq) + nk2mkq
         enddo
      enddo
#if defined(SCAT_FWD)
   enddo
#endif
!$omp end parallel do
   !combine contributions from all the pools/processors.
   call mp_sum(epint, inter_pool_comm)
   !call mp_barrier(inter_pool_comm)
   call stopwatch('trans_scat_int','stop')
end subroutine trans_scat_int

! Real time dynamics, Refer to eq.17. 18 in Eur. Phys. J. B (2016) 89: 239
! Bernardi, First-principles dynamics of electrons and phonons
subroutine cdyna_scat_int(kg, eptemp, dist, epcol)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, dist(kg%numb, kg%nk)
   ! e-ph collision term computed from given distribution function dist(:,:)
   real(dp), intent(out) :: epcol(kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, im, n, m, bands(3)
   integer(kind=i8b) :: ns
   real(dp):: wq, enk, emkq, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq
#if !defined(STORE_G2DT)
   real(dp) :: g2, dt1, dt2
#else
   real(dp) :: g2_dt1, g2_dt2
#endif

   epcol(:,:) = 0.0_dp

#if defined(_OPENACC)
! OpenACC acceleration
!$acc kernels
#else
! OpenMP acceleration
!$omp parallel do schedule(guided) &
#if defined(CDYN_USE_ARRAY_REDUCE)
!$omp&  reduction(+:epcol) &
#endif
!$omp&  default(shared) private(i, iq, ik, ikq, ns, bands, im, n, m, wq, enk, emkq, &
#if !defined(STORE_G2DT)
!$omp&  g2, dt1, dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq)
#else
!$omp&  g2_dt1, g2_dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq)
#endif
#endif
#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      ! loop over all scatter-channels, using the scat_index to lookup
      ! the corresponding (iq, ik, ikq) values
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! phonon occupation
         wq   = qgrid%freq(im, iq)
         ph_n = bose(eptemp, wq) ! N_q\mu

         ! electron distributaion
         fnk  = dist(n, ik)
         fmkq = dist(m, ikq)
         ! scattering process:  formula 1
         !fabs = fnk*(one - fmkq)*wgq - fmkq*(one - fnk)*(one + wgq)
         !fem  = fnk*(one - fmkq)*(one + wgq) - fmkq*(one - fnk)*wgq
         ! formula 2
         fabs = fnk*(ph_n + fmkq) - fmkq*(ph_n + 1.0_dp)
         fem  = fnk*(ph_n + 1.0_dp ) - fmkq*(ph_n + fnk)

         ! update nk with contribution from k+q
#if defined(STORE_G2DT)
         ! Use the precomputed -g2*dt1 and -g2*dt2 terms
         mkq2nk = -(scatter_channels(ns)%g2_dt1 * fabs + &
                    scatter_channels(ns)%g2_dt2 * fem)
#else
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2

         ! eigenvalue
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)
         if ( adapt_smear_eph ) then
            ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
            ! delta(Enk - Emkq + wq)
            dt1 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(ns), enk - emkq + wq)
            ! delta(Enk - Emkq -wq)
            dt2 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(ns), enk - emkq - wq)
         else
            ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
            ! delta(Enk - Emkq + wq)
            dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
            ! delta(Enk - Emkq -wq)
            dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
         endif

         mkq2nk = - g2 * (dt1*fabs + dt2*fem)
#endif
         !update nk with contribution from k+q
#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
         epcol(n, ik) = epcol(n, ik) + mkq2nk
         ! the inverse: update mk+q with contribution from k
         nk2mkq = - mkq2nk
#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
         epcol(m, ikq) = epcol(m, ikq) + nk2mkq
      enddo   !; enddo; enddo ! ns loop
#if defined(SCAT_FWD)
   enddo ! i loop
#endif
#if defined(_OPENACC)
!$acc end kernels
#else
!$omp end parallel do
#endif
   ! combine contributions from all the pools/processes.
   call mp_sum(epcol, inter_pool_comm)
!   call mp_barrier(inter_pool_comm)
end subroutine cdyna_scat_int

! Real time dynamics, Refer to eq.17. 18 in Eur. Phys. J. B (2016) 89: 239
! Bernardi, First-principles dynamics of electrons and phonons
subroutine cphdyna_scat_int(kg, qg, dist, epcol, disp, pecol)
   use boltz_utils, only: kpt2num, num2kpt, kpts_minus, num2ipt
   implicit none
   type(grid), intent(in) :: kg, qg
   real(dp), intent(in) :: dist(kg%numb, kg%nk), disp(qg%numb, qg%nk)
   ! e-ph collision term computed from given distribution function dist(:,:)
   real(dp), intent(out) :: epcol(kg%numb, kg%nk), pecol(qg%numb, qg%nk)
   !local variables
   integer :: i, iq, inq, ik, ikq, im, n, m, bands(3)
   integer(kind=i8b) :: ns
   integer :: iqnums(3)
   real(dp):: wq, wnq, enk, emkq, ph_n, ph_nn, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq
   real(dp) :: imph, imphn
   real(dp) :: weight_ratio, smear_tmp, ediff_abs, ediff_ems 
   real(dp) :: factor, dscat

!#if !defined(STORE_G2DT)
   real(dp) :: g2, dt1, dt2
!#else
   real(dp) :: g2_dt1, g2_dt2
!#endif

   call stopwatch('cphdyna_scat_int','start')
   epcol(:,:) = 0.0_dp
   pecol(:,:) = 0.0_dp
   weight_ratio = kg%kweight/qg%kweight
   factor = qg%kweight * pi * 2.0_dp
   
!#if defined(_OPENACC)
!!OpenACC acceleration
!!$acc kernels
!#else
! OpenMP acceleration
!$omp parallel do schedule(guided) &
!#if defined(CDYN_USE_ARRAY_REDUCE)
!!$omp&  reduction(+:epcol) &
!#endif
!$omp&  default(shared) private(i, iq, ik, ikq, ns, bands, im, n, m, wq, enk, emkq, &
!#if !defined(STORE_G2DT)
!$omp&  g2, dt1, dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq, &
!$omp&  wnq, ph_nn, imph, imphn, smear_tmp, ediff_abs, ediff_ems, dscat, iqnums, inq)
!#else
!!$omp&  g2_dt1, g2_dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq)
!#endif
!#endif
#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
      iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      ! loop over all scatter-channels, using the scat_index to lookup
      ! the corresponding (iq, ik, ikq) values
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
         iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif
         !sypeng: unsure: maybe need to be modified for gpu 
         inq = kpt2num(-1.0_dp*num2kpt(iq-1, qg%ndim), qg%ndim) + 1

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! phonon occupation
         wq   = qg%enk(im, iq)
         wnq  = qg%enk(im, inq)
         ph_n = disp(im, iq) ! N_q\mu
         ph_nn = disp(im, inq) ! N_q\mu
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)

         ! electron distributaion
         fnk  = dist(n, ik)
         fmkq = dist(m, ikq)
         if ( adapt_smear_eph ) then
            smear_tmp = tsmear_list(i)
         else
            smear_tmp = delta_smear
         endif
         ediff_abs = abs(enk - emkq + wq)
         ediff_ems = abs(enk - emkq - wnq)
         dscat = factor * scatter_channels(ns)%eph_g2
         !
         if (ediff_abs < 3.0_dp*smear_tmp) then 
            dt1 = gauss(smear_tmp, ediff_abs)
            ! formula 2
            ph_n = disp(im, iq) ! N_q\mu
            fabs = fnk*(ph_n + fmkq) - fmkq*(ph_n + 1.0_dp)
            imph = - dscat* dt1 * fabs !(w0g1*fabs - w0g2*fem)
         else
            imph = 0.0_dp
         endif
         if (ediff_ems < 3.0_dp*smear_tmp) then 
            dt2 = gauss(smear_tmp, ediff_ems)
            ph_nn = disp(im, inq) ! N_q\mu
            fem  = fnk*(ph_nn + 1.0_dp ) - fmkq*(ph_nn + fnk)
            imphn = dscat* dt2 * fem !(w0g1*fabs - w0g2*fem)
         else
            imphn = 0.0_dp
         endif
        !mkq2nk = - g2 * (dt1*fabs + dt2*fem)
         mkq2nk = imph - imphn
         nk2mkq = - mkq2nk
         imph = imph * weight_ratio 
         imphn = imphn * weight_ratio 
!#endif
         !update nk with contribution from k+q
!#if defined(_OPENACC)
!!$acc atomic
!#else
!#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
!#endif
!#endif
         epcol(n, ik) = epcol(n, ik) + mkq2nk
         ! the inverse: update mk+q with contribution from k
        !nk2mkq = - mkq2nk
!#if defined(_OPENACC)
!!$acc atomic
!#else
!#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
!#endif
!#endif
         epcol(m, ikq) = epcol(m, ikq) + nk2mkq
!$omp atomic
         pecol(im, iq) = pecol(im, iq) + imph
         ! this is only for 2D. for 3D, need to add another loop
         iqnums =  num2ipt(qg%kpt(iq), qg%ndim)
         if ( product(iqnums - boltz_qdim/2) .eq. 0 ) then
            cycle
         else
!$omp atomic
            pecol(im, inq) = pecol(im, inq) + imphn
         endif
      enddo   !; enddo; enddo ! ns loop
#if defined(SCAT_FWD)
   enddo ! i loop
#endif
!#if defined(_OPENACC)
!!$acc end kernels
!#else
!$omp end parallel do
!#endif
   ! combine contributions from all the pools/processes.
   call mp_sum(epcol, inter_pool_comm)
   call mp_sum(pecol, inter_pool_comm)
   call stopwatch('cphdyna_scat_int','stop')
end subroutine cphdyna_scat_int

subroutine ph3_scat_int(qg, disp, ph3col, scat_flag)
   use boltz_utils, only: kpts_plus_invert
   implicit none
   type(grid), intent(in) :: qg
   real(dp), intent(in) :: disp(qg%numb, qg%nk)
   real(dp), intent(out) :: ph3col(qg%numb, qg%nk)
   logical, intent(in), optional :: scat_flag
   !local variables
   logical :: scat_flag_loc
   integer(kind=i8b) :: is
   integer :: ivec, iq0, iq1, iq2, it, im0, im1, im2, bands(3)
   real(dp) :: wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat
   real(dp) :: occ_sum_0, occ_sum_1, occ_sum_2, factor
   real(dp) :: qweight, e_thr
   integer :: scat_index
   real(dp), allocatable :: ph3col_tmp(:,:)
   call start_clock('phph_scatint')
   scat_flag_loc = .false.
   if ( present(scat_flag) ) then
      if ( scat_flag ) scat_flag_loc = scat_flag
   endif
   
   qweight = 1.0_dp / (boltz_qdim(1) * boltz_qdim(2) * boltz_qdim(3))
   e_thr = delta_smear_ph*3.0_dp
   ! this constant means that the output is Im Simga only
   ! to convert to tau^-1 needs to multiply results by 2/hbar
   factor = 9.0E0_dp * 2.0E0_dp * pi * qweight 
   if(symm)then
      allocate(ph3col_tmp(qg%numb, qg%nk_irr))
   else
      allocate(ph3col_tmp(qg%numb, qg%nk))
   endif
   ph3col(:,:) = 0.0E0_dp
   ph3col_tmp(:,:) = 0.0E0_dp
!$omp parallel do schedule (guided) default(shared) private(is, iq0, iq1, iq2, bands, im0, im1, im2, &
!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat, occ_sum_0, occ_sum_1, occ_sum_2, &
!$omp& scat_index) &
!$omp& reduction(+:ph3col_tmp)
   do is = 1, ph3_num_scat
      ! key step, to keep g2 the same 
      g2 = tscat_ph(is)%g2
      ! ignore small values
      if ( g2 < 1.0E-30_dp ) then
         cycle
      endif
      scat_index = tscat_ph(is)%scat_index
      iq0 = scatter_ph(1, scat_index)
      iq1 = scatter_ph(2, scat_index)
      iq2 = kpts_plus_invert(qg%kpt(iq0), qg%kpt(iq1), qg%ndim) + 1
      ! band index is defined the same, here we use qg, but qg_fine%numb is the same 
      bands = idx2band( tscat_ph(is)%bands_index, qg%numb )
      im2 = bands(1); im0 = bands(2); im1 = bands(3)
         
      ! eigenvalue
      wq0 = qg%enk(im0, iq0); wq1 = qg%enk(im1, iq1); wq2 = qg%enk(im2, iq2)
      ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
      !> d3q  add the factor of 18pi/hbar^2/36, 0.5^3 is included in g2, hbar is omitted in atomic units.
      !> tdep add prefactor of hbar*pi/2, 1/8 is preincluded in g2
      ! delta(wq + wq1 - wq2)
        
      ! halving equal ik and iq
      !if ((im0 .eq. im1) .and. (iq0 .eq. iq1)) then
      !   e_delta = e_delta/2.0_dp
      !endif
      ! phonon occupation
      ocq0 = disp(im0, iq0)
      ocq1 = disp(im1, iq1)
      ocq2 = disp(im2, iq2)
      !syp - I will split processes into absorption and emission process
      !    - this is a little different with electron strategy, we do it strickly
      !    - for each i_scatt, either absorption or emission, not both.
      ! 1. absorption 
      if (abs(wq0 + wq1 - wq2) < e_thr)then  
         if ( adapt_smear_phph ) then
            e_delta = gauss(tsmear_ph_list(is), wq1 + wq0 - wq2)
         else
            e_delta = gauss(delta_smear_ph, wq1 + wq0 - wq2)
         endif

         !syp - the factor of 2 is for two absorption processes
         if (scat_flag_loc) then
            dscat = g2 * e_delta * factor * 2.0_dp * (ocq1 - ocq2)
         else
            dscat = g2 * e_delta * factor * 2.0d0 * (ocq2 * (ocq0 + ocq1 + 1.0_dp) - ocq0 * ocq1)
         endif
      
      ! 3. emission 
      elseif(abs(wq0 - wq1 - wq2) < e_thr) then  
         if ( adapt_smear_phph ) then
            e_delta = gauss(tsmear_ph_list(is), wq0 - wq1 - wq2)
         else
            e_delta = gauss(delta_smear_ph, wq0 - wq1 - wq2)
         endif
         if (scat_flag_loc) then
            dscat = g2 * e_delta * factor * (1.d0 + ocq1 + ocq2)
         else
            dscat = -g2 * e_delta * factor * 2.0d0 * (ocq2 * (ocq0 + ocq1 + 1.0_dp) - ocq0 * ocq1)
         endif
      endif

      ph3col_tmp(im0,qg%kpt2ir(iq0)) = ph3col_tmp(im0,qg%kpt2ir(iq0)) + dscat
   enddo ! is loop
!$omp end parallel do

   call mp_sum(ph3col_tmp, inter_pool_comm)

   !map irreducible points to full BZ data
   if(symm)then
       do iq0 = 1, qg%nk
          ph3col(:, iq0) = ph3col_tmp(:, qg%kpt2ir(iq0))
       enddo
   else
      ph3col = ph3col_tmp
   endif

   !call mp_sum(ph3col, inter_pool_comm)
   call stop_clock('phph_scatint')
end subroutine ph3_scat_int

subroutine ph3_scat_int_interp_half(qg, disp, ph3col, scat_flag)
   use boltz_utils, only: kpts_plus_invert, num2kpt, kpt2num, fold_k, &
      num2ipt
   implicit none
   type(grid), intent(in) :: qg
   real(dp), intent(in) :: disp(qg%numb, qg%nk)
   real(dp), intent(out) :: ph3col(qg%numb, qg%nk)
   logical, intent(in), optional :: scat_flag
   !local variables
   logical :: scat_flag_loc
   integer(kind=i8b) :: is
   integer :: ivec, iq0, iq1, iq2, it, im0, im1, im2, bands(3)
   real(dp) :: wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat
   real(dp) :: occ_sum_0, occ_sum_1, occ_sum_2, factor
   real(dp) :: qweight, e_thr
   integer :: scat_index
   real(dp), allocatable :: ph3col_tmp(:,:)
   integer :: xinc, yinc, zinc, iq0c, iq1c
   integer :: ipt0c(3), ipt1c(3), ipt0(3), ipt1(3)
   integer :: dim_vec(3), pos_vec(3)

   call start_clock('phph_scatint')

   dim_vec(1) = qg%ndim(2) * qg%ndim(3)
   dim_vec(2) = qg%ndim(3)
   dim_vec(3) = 1

   scat_flag_loc = .false.
   if ( present(scat_flag) ) then
      if ( scat_flag ) scat_flag_loc = scat_flag
   endif
   
   qweight = 1.0_dp / (boltz_qdim(1) * boltz_qdim(2) * boltz_qdim(3))
   e_thr = delta_smear_ph*3.0_dp

   ! this constant means that the output is Im Simga only
   ! to convert to tau^-1 needs to multiply results by 2/hbar
   factor = 9.0E0_dp * 2.0E0_dp * pi * qweight 
   if(symm)then
      allocate(ph3col_tmp(qg%numb, qg%nk_irr))
   else
      allocate(ph3col_tmp(qg%numb, qg%nk))
   endif
   ph3col(:,:) = 0.0E0_dp
   ph3col_tmp(:,:) = 0.0E0_dp
!$omp parallel do schedule (guided) default(shared) private(is, iq0, iq1, iq2, bands, im0, im1, im2, &
!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat, occ_sum_0, occ_sum_1, occ_sum_2, &
!$omp& scat_index, &
!$omp& iq1c, ipt1c, ipt1, xinc, yinc, zinc, pos_vec) &
!$omp& reduction(+:ph3col_tmp)
   do is = 1, ph3_num_scat
      ! key step, to keep g2 the same 
      g2 = tscat_ph(is)%g2
      ! ignore small values
      if ( g2 < 1.0E-30_dp ) then
         cycle
      endif
      scat_index = tscat_ph(is)%scat_index
      !iq0 = scatter_ph(1, scat_index)
      !iq1 = scatter_ph(2, scat_index)
      !iq2 = kpts_plus_invert(qg%kpt(iq0), qg%kpt(iq1), qg%ndim) + 1
      ! band index is defined the same, here we use qg, but qg_fine%numb is the same 

      bands = idx2band( tscat_ph(is)%bands_index, qg%numb )
      im2 = bands(1); im0 = bands(2); im1 = bands(3)

      iq0 = scatter_ph(1, scat_index)
      iq1c = scatter_ph(2, scat_index)
      ipt1c = num2ipt(qg%kpt(iq1c), qg%ndim)
      ! directly find corresponding indices
      ! see where the point is
      do xinc = -1, 1
      do yinc = -1, 1
      do zinc = -1, 1 
         pos_vec = (/xinc, yinc, zinc/)
         ! compute new iq1, iq2
         ipt1 = mod(ipt1c + qg%ndim + pos_vec, qg%ndim)
         iq1 = dot_product(ipt1, dim_vec) + 1
         iq2 = kpts_plus_invert(qg%kpt(iq0), qg%kpt(iq1), qg%ndim) + 1
         
         ! eigenvalue
         wq0 = qg%enk(im0, iq0); wq1 = qg%enk(im1, iq1); wq2 = qg%enk(im2, iq2)
         ! phonon occupation
         ocq0 = disp(im0, iq0)
         ocq1 = disp(im1, iq1)
         ocq2 = disp(im2, iq2)
         !syp - I will split processes into absorption and emission process
         !    - this is a little different with electron strategy, we do it strickly
         !    - for each i_scatt, either absorption or emission, not both.
         ! 1. absorption 
         if (abs(wq0 + wq1 - wq2) < e_thr)then  
            if ( adapt_smear_phph ) then
               e_delta = gauss(tsmear_ph_list(is), wq1 + wq0 - wq2)
            else
               e_delta = gauss(delta_smear_ph, wq1 + wq0 - wq2)
            endif

            !syp - the factor of 2 is for two absorption processes
            if (scat_flag_loc) then
               dscat = g2 * e_delta * factor * 1.0_dp * (ocq1 - ocq2)
            else
               dscat = g2 * e_delta * factor * 2.0_dp * (ocq2 * (ocq0 + ocq1 + 1.0_dp) - ocq0 * ocq1)
            endif

         elseif (abs(wq0 - wq1 + wq2) < e_thr) then  
            if ( adapt_smear_phph ) then
               e_delta = gauss(tsmear_ph_list(is), wq2 + wq0 - wq1)
            else
               e_delta = gauss(delta_smear_ph, wq2 + wq0 - wq1)
            endif

            !syp - the factor of 2 is for two absorption processes
            if (scat_flag_loc) then
               dscat = g2 * e_delta * factor * 1.0_dp * (ocq2 - ocq1)
            else
               dscat = g2 * e_delta * factor * 2.0_dp * (ocq1 * (ocq0 + ocq2 + 1.0_dp) - ocq0 * ocq2)
            endif
         
         ! 3. emission 
         elseif(abs(wq0 - wq1 - wq2) < e_thr) then  
            if ( adapt_smear_phph ) then
               e_delta = gauss(tsmear_ph_list(is), wq0 - wq1 - wq2)
            else
               e_delta = gauss(delta_smear_ph, wq0 - wq1 - wq2)
            endif
            if (scat_flag_loc) then
               dscat = g2 * e_delta * factor * (1.0_dp + ocq1 + ocq2)
            else
               dscat = -g2 * e_delta * factor * 2.0_dp * (ocq0 * (ocq1 + ocq2 + 1.0_dp) - ocq1 * ocq2)
            endif
         endif

         ph3col_tmp(im0,qg%kpt2ir(iq0)) = ph3col_tmp(im0,qg%kpt2ir(iq0)) + dscat
      enddo
      enddo
      enddo
   enddo ! is loop
!$omp end parallel do

   call mp_sum(ph3col_tmp, inter_pool_comm)

   !map irreducible points to full BZ data
   if(symm)then
       do iq0 = 1, qg%nk
          write(*,*) iq0, qg%kpt2ir(iq0)
          ph3col(:, iq0) = ph3col_tmp(:, qg%kpt2ir(iq0))
       enddo
   else
      ph3col = ph3col_tmp
   endif

   !call mp_sum(ph3col, inter_pool_comm)
   call stop_clock('phph_scatint')
end subroutine ph3_scat_int_interp_half

subroutine ph3_scat_int_interp2(qg, disp, ph3col, scat_flag)
   use boltz_utils, only: kpts_plus_invert, num2kpt, kpt2num, fold_k, &
      num2ipt
   use pert_utils, only : find_free_unit
   implicit none
   type(grid), intent(in) :: qg
   real(dp), intent(in) :: disp(qg%numb, qg%nk)
   real(dp), intent(out) :: ph3col(qg%numb, qg%nk)
   logical, intent(in), optional :: scat_flag
   !local variables
   logical :: scat_flag_loc
   integer(kind=i8b) :: is
   integer :: ivec, iq0, iq1, iq2, it, im0, im1, im2, bands(3)
   real(dp) :: wq0, wq1, wq2, wdiff, g2, e_delta, ocq0, ocq1, ocq2, dscat
   real(dp) :: occ_sum_0, occ_sum_1, occ_sum_2, factor
   real(dp) :: qweight
   real(dp), allocatable :: ph3col_tmp(:,:)
   integer :: xinc, yinc, zinc, iq0c, iq1c
   integer :: xinc1, yinc1, zinc1
   integer :: ipt0c(3), ipt1c(3), ipt0(3), ipt1(3)
   integer :: dim_vec(3), pos_vec(6)
   integer :: scat_index
   real(dp) :: fact
   
   character(len=120) :: fname, dset_name, msg

   call start_clock('phph_scatint')

   ! try using tscat_ph
   dim_vec(1) = qg%ndim(2) * qg%ndim(3)
   dim_vec(2) = qg%ndim(3)
   dim_vec(3) = 1
   scat_flag_loc = .false.
   if ( present(scat_flag) ) then
      if ( scat_flag ) scat_flag_loc = scat_flag
   endif
   qweight = 1.0_dp / (boltz_qdim(1) * boltz_qdim(2) * boltz_qdim(3))
   ! this constant means that the output is Im Simga only
   ! to convert to tau^-1 needs to multiply results by 2/hbar
   factor = 9.0E0_dp * 4.0E0_dp * pi * qweight 
   ph3col(:,:) = 0.0E0_dp
   fact = 1.0_dp
   allocate(ph3col_tmp(qg%numb, qg%nk))
!$omp parallel default(shared) private(is, iq0, iq1, iq2, bands, im0, im1, im2, &
!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat, occ_sum_0, occ_sum_1, occ_sum_2, &
!$omp& iq0c, iq1c, ipt0c, ipt1c, ipt0, ipt1, xinc, yinc, zinc, xinc1, yinc1, zinc1, scat_index, &
!$omp& pos_vec, ph3col_tmp, wdiff) 
   ph3col_tmp(:,:) = 0.0E0_dp
!$omp do schedule(guided)
   do is = 1, ph3_num_scat
      ! key step, to keep g2 the same 
      g2 = tscat_ph(is)%g2
      ! ignore small values
      if ( g2 < 1.0E-30_dp ) then
         cycle
      endif
      scat_index = tscat_ph(is)%scat_index
      iq0c = scatter_ph(1, scat_index)
      iq1c = scatter_ph(2, scat_index)
      ! band index is defined the same, here we use qg, but qg_fine%numb is the same 
      bands = idx2band( tscat_ph(is)%bands_index, qg%numb )
      im2 = bands(1); im0 = bands(2); im1 = bands(3)

      ipt0c = num2ipt(qg%kpt(iq0c), qg%ndim)
      ipt1c = num2ipt(qg%kpt(iq1c), qg%ndim)
      ! directly find corresponding indices
      ! see where the point is
      do xinc = 0, 1
      do yinc = 0, 1
      do zinc = 0, 1 
      do xinc1 = 0, 1 
      do yinc1 = 0, 1 
      do zinc1 = 0, 1 
         pos_vec = (/xinc, yinc, zinc, xinc1, yinc1, zinc1/)
         !fact = (0.5_dp) ** count(pos_vec .ne. 0)
         ipt0 = mod(ipt0c + qg%ndim + pos_vec(1:3), qg%ndim)
         ipt1 = mod(ipt1c + qg%ndim + pos_vec(4:6), qg%ndim)
         ! compute iq0, iq1, iq2
         iq0 = dot_product(ipt0, dim_vec) + 1
         iq1 = dot_product(ipt1, dim_vec) + 1
         iq2 = dot_product(mod(2*qg%ndim-ipt0-ipt1, qg%ndim), dim_vec) + 1
         
         ! eigenvalue
         wq0 = qg%enk(im0, iq0); wq1 = qg%enk(im1, iq1); wq2 = qg%enk(im2, iq2)
         wdiff = abs(- wq2 + wq1 + wq0)
         if (wdiff > 3.0_dp*delta_smear_ph) cycle

         if ( adapt_smear_phph ) then
            e_delta = gauss(tsmear_ph_list(is), wdiff)
         else
            e_delta = gauss(delta_smear_ph, wdiff)
         endif
         e_delta = 1.0_dp
         ! halving equal ik and iq
         if ((im0 .eq. im1) .and. (iq0 .eq. iq1)) then
            e_delta = e_delta/2.0_dp
         endif
         !! phonon occupation
         ocq0 = disp(im0, iq0)
         ocq1 = disp(im1, iq1)
         ocq2 = disp(im2, iq2)
         if ( scat_flag_loc ) then
            ! scattering process:
            dscat = g2 * e_delta * factor * fact
            occ_sum_0 =  dscat*(ocq1 - ocq2)
            occ_sum_1 =  dscat*(ocq0 - ocq2)
            occ_sum_2 =  dscat*(ocq0 + ocq1 + 1.0_dp)
            ph3col_tmp(im0, iq0) = ph3col_tmp(im0, iq0) + occ_sum_0
            ph3col_tmp(im1, iq1) = ph3col_tmp(im1, iq1) + occ_sum_1
            ph3col_tmp(im2, iq2) = ph3col_tmp(im2, iq2) + occ_sum_2
         else
            ! dynamics - another factor of 2
            dscat = 2.0_dp * g2 * e_delta * factor * fact
            occ_sum_0 = dscat * (ocq2 * (ocq0 + ocq1 + 1.0_dp) - ocq0 * ocq1)
            ph3col_tmp(im0, iq0) = ph3col_tmp(im0, iq0) + occ_sum_0
            ph3col_tmp(im1, iq1) = ph3col_tmp(im1, iq1) + occ_sum_0
            ph3col_tmp(im2, iq2) = ph3col_tmp(im2, iq2) - occ_sum_0
         endif
      enddo; enddo; enddo! iqf loop
      enddo; enddo; enddo! iqf loop
   enddo ! is loop
!$omp end do
do iq0 = 1, qg%nk
   do im0 = 1, qg%numb
!$omp atomic
      ph3col(im0, iq0) = ph3col(im0, iq0) + ph3col_tmp(im0, iq0)
   enddo
enddo
!$omp end parallel

   call mp_sum(ph3col, inter_pool_comm)
   call stop_clock('phph_scatint')
end subroutine ph3_scat_int_interp2
end module boltz_scatter_integral
