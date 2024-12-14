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

module boltz_scatter_integral
   use omp_lib
   use pert_const, only: dp, pi
   use boltz_grid, only: grid
   use boltz_utils, only: idx2band
   use pert_utils, only: bose, fermi, gauss
   use pert_param, only: delta_smear, delta_smear_ph, boltz_qdim, &
      adapt_smear_eph, adapt_smear_phph, phfreq_cutoff, phfreq_cutoff_ph
   use pert_output, only: stopwatch
   use boltz_scatter, only: num_kq_pair, num_qq_pair, scat, ph3_scat, tsmear_ph, qgrid, &
      scatter_e, scatter_ph, num_scat, ph3_num_scat, tscat, tscat_ph, tsmear_list, tsmear_ph_list
   use qe_mpi_mod, only: stdout, ionode, mp_barrier, mp_sum, inter_pool_comm, my_pool_id
   implicit none
   private 
   public :: rates_scat_int, trans_scat_int, cdyna_scat_int
   public :: cphdyna_scat_int, ph3_scat_int
   public :: ph3_scat_int_interp2
contains

! compute scattering rate
subroutine rates_scat_int(kg, eptemp, efermi, rates)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, efermi
   real(dp), intent(out) :: rates(:,:) !rates(kg%numb, kg%numk)
   !local variables
   integer :: i, iq, ik, ikq, ns, im, n, m, bands(3)
   real(dp) :: wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, nk2mkq, mkq2nk
   integer :: scat_index
 
   call stopwatch('rates_scat_int', 'start')
   rates(:,:) = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns, &
!$omp& im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, mkq2nk, nk2mkq, bands, scat_index)
   do i = 1, num_scat
      scat_index = tscat(i)%scat_index
      ik = scatter_e(1, scat_index)
      ikq = scatter_e(2, scat_index)
      iq = scatter_e(3, scat_index)
      !iq = scat(i)%iq;    ik = scat(i)%ik;    ikq = scat(i)%ikq
      ! map index to (mkq, nk, mu)
      bands = idx2band( tscat(i)%bands_index, kg%numb )
      m = bands(1);  n = bands(2);  im = bands(3)
      ! eigenvalue
      wq   = qgrid%freq(im, iq)
      enk  = kg%enk(n, ik)
      emkq = kg%enk(m, ikq)
      ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
      g2 = tscat(i)%g2
         
      if ( adapt_smear_eph ) then
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq + wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq - wq)
         !dt1 = qgrid%weight * pi * 2.0_dp * gauss(scat(i)%adsmear(ns), enk - emkq + wq)
         !dt2 = qgrid%weight * pi * 2.0_dp * gauss(scat(i)%adsmear(ns), enk - emkq - wq)
      else
         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq) 
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
      endif
      ! phonon occupation
      ! the constrain of 'wq > eps_acustic' is aready done in scat_setup.
      ph_n = bose(eptemp, wq) ! N_q\mu

      ! compute scattering rates
      fkq = fermi(eptemp, emkq-efermi)
      mkq2nk = g2 * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
! lock memory when updating this variable
!$omp atomic
      rates(n, ik)  = rates(n, ik)  + mkq2nk

      !update mk+q with contribution from k
      fkq = fermi(eptemp, enk-efermi)
      !delta(Ekq-Ek+wq)=delta(Ek-Ekq-wq); delta(Ekq-Ek-wq)=delta(Ek-Ekq+wq)
      nk2mkq = g2 * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
!$omp atomic
      rates(m, ikq) = rates(m, ikq) + nk2mkq
   enddo ! i loop
!$omp end parallel do
   ! combine contributions from all the pools/processors.
   call mp_sum(rates, inter_pool_comm)
!   call mp_barrier(inter_pool_comm)
   call stopwatch('rates_scat_int', 'stop')
end subroutine rates_scat_int

! iterative process for tranport calculation
subroutine trans_scat_int(kg, eptemp, efermi, mfd, epint)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: efermi, eptemp, mfd(:,:,:) !mfd(3 or 6, kg%numb, kg%nk)
   real(dp), intent(out) :: epint(:,:,:)  !epint(3 or 6, kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, ns, im, n, m, ic, bands(3), ncomp
   real(dp) :: wq, enk, emkq, g2, dt1, dt2, ph_n, nk2mkq, mkq2nk, fkq
   integer :: scat_index

   call stopwatch('trans_scat_int','start')
   !sanity check
   ncomp = size(mfd, 1)
   if( ncomp .ne. size(epint, 1) ) &
      call errore('trans_scat_int',"mismatch dimension in input arguments!",1)

   epint(:,:,:) = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns,&
!$omp&  im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, mkq2nk, nk2mkq, ic, bands, scat_index)
   do i = 1, num_scat
      scat_index = tscat(i)%scat_index
      ik = scatter_e(1, scat_index)
      ikq = scatter_e(2, scat_index)
      iq = scatter_e(3, scat_index)
      !iq = scat(i)%iq;   ik = scat(i)%ik;   ikq = scat(i)%ikq
      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      bands = idx2band( tscat(i)%bands_index, kg%numb )
      !   bands = idx2band( scat(i)%bands_index(ns), kg%numb )
      m = bands(1);  n = bands(2);  im = bands(3)
      ! eigenvalue
      wq   = qgrid%freq(im, iq)
      enk  = kg%enk(n, ik)
      emkq = kg%enk(m, ikq)
      ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
      g2 = tscat(i)%g2
         
      if ( adapt_smear_eph ) then
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq + wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq - wq)
        !dt1 = qgrid%weight * pi * 2.0_dp * gauss(scat(i)%adsmear(ns), enk - emkq + wq)
        !dt2 = qgrid%weight * pi * 2.0_dp * gauss(scat(i)%adsmear(ns), enk - emkq - wq)
      else
         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq) 
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
      endif
      ! phonon occupation
      ! the constrain of 'wq > eps_acustic' is aready done in scat_setup.
      ph_n = bose(eptemp, wq) ! N_q\mu

      do ic = 1, ncomp
         !update nk with contribution from k+q
         fkq = fermi(eptemp, emkq-efermi)
         mkq2nk = g2 * mfd(ic, m, ikq) * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
!$omp atomic
         epint(ic, n, ik)  = epint(ic, n, ik)  + mkq2nk
        
         ! update mk+q with contribution from k
         fkq = fermi(eptemp, enk-efermi)
         nk2mkq = g2 * mfd(ic, n, ik)  * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
!$omp atomic
         epint(ic, m, ikq) = epint(ic, m, ikq) + nk2mkq
      enddo
   enddo
!$omp end parallel do
   !combine contributions from all the pools/processors.
   call mp_sum(epint, inter_pool_comm)
   !call mp_barrier(inter_pool_comm)
   call stopwatch('trans_scat_int','stop')
end subroutine trans_scat_int

! Real time dynamics, Refer to eq.17. 18 in Eur. Phys. J. B (2016) 89: 239
! Bernardi, First-principles dynamics of electrons and phonons
subroutine cdyna_scat_int(kg, eptemp, dist, epcol)
   use boltz_utils, only: kpts_minus
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, dist(kg%numb, kg%nk)
   ! e-ph collision term computed from given distribution function dist(:,:)
   real(dp), intent(out) :: epcol(kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, ns, im, n, m, bands(3)
   real(dp):: wq, enk, emkq, g2, dt1, dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq
   integer :: scat_index
   epcol(:,:) = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, &
!$omp&  ns, im, n, m, wq, bands, enk, emkq, g2, dt1, dt2, ph_n, fnk, fmkq, &
!$omp&  fabs, fem, mkq2nk, nk2mkq, scat_index) 
   do i = 1, num_scat
      scat_index = tscat(i)%scat_index
      ik = scatter_e(1, scat_index)
      ikq = scatter_e(2, scat_index)
      iq = scatter_e(3, scat_index)
      ! map index to (mkq, nk, mu)
      bands = idx2band( tscat(i)%bands_index, kg%numb )
      !bands = idx2band( scat(i)%bands_index(ns), kg%numb )
      m = bands(1);  n = bands(2);  im = bands(3)
      ! eigenvalue
      wq   = qgrid%freq(im, iq)
      enk  = kg%enk(n, ik)
      emkq = kg%enk(m, ikq)
      ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
      g2 = tscat(i)%g2
      if ( adapt_smear_eph ) then
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq + wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq - wq)
      else
         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq) 
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
      endif

      ! phonon occupation
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
      mkq2nk = - g2 * (dt1*fabs + dt2*fem)
      !update nk with contribution from k+q
!$omp atomic
      epcol(n, ik) = epcol(n, ik) + mkq2nk
      ! the inverse: update mk+q with contribution from k
      nk2mkq = - mkq2nk
!$omp atomic
      epcol(m, ikq) = epcol(m, ikq) + nk2mkq
   enddo ! i loop
!$omp end parallel do
   ! combine contributions from all the pools/processes.
   call mp_sum(epcol, inter_pool_comm)
!   call mp_barrier(inter_pool_comm)
end subroutine cdyna_scat_int

subroutine cphdyna_scat_int(kg, qg, dist, epcol, disp, pecol)
   use boltz_utils, only: kpt2num, num2kpt, kpts_minus, num2ipt
   implicit none
   type(grid), intent(in) :: kg, qg
   real(dp), intent(in) :: dist(kg%numb, kg%nk), disp(qg%numb, qg%nk)
   ! e-ph collision term computed from given distribution function dist(:,:)
   real(dp), intent(out) :: epcol(kg%numb, kg%nk), pecol(qg%numb, qg%nk)
   !local variables
   integer :: i, iq, inq, ik, ikq, im, n, m, bands(3)
   real(dp):: wq, wnq, enk, emkq, g2, dt1, dt2, ph_n, ph_nn
   real(dp) :: fnk, fmkq, fabs, fem, mkq2nk, nk2mkq, imph, imphn
   integer :: scat_index, iqnums(3)
   real(dp) :: weight_ratio, smear_tmp, ediff_abs, ediff_ems
   real(dp) :: factor, dscat
   ! totalq is the number of total q points
   call start_clock('eph_scatint')
   epcol(:,:) = 0.0E0_dp; pecol(:,:) = 0.0E0_dp
   weight_ratio = kg%kweight/qg%kweight
   factor = qg%kweight * pi * 2.0_dp
!$omp parallel default(shared) private(i, iq, ik, ikq, &
!$omp& im, n, m, wq, bands, enk, emkq, g2, dt1, dt2, ph_n, ph_nn, fnk, fmkq, &
!$omp& fabs, fem, mkq2nk, nk2mkq, imph, imphn, inq, wnq, scat_index, iqnums, &
!$omp& ediff_abs, ediff_ems, smear_tmp, dscat) 
!$omp do schedule(guided)
   do i = 1, num_scat
      scat_index = tscat(i)%scat_index
      ik = scatter_e(1, scat_index)
      ikq = scatter_e(2, scat_index)
      ! correct iq index
      iq = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim) + 1
      inq = kpt2num(-1.0_dp*num2kpt(iq-1, qg%ndim), qg%ndim) + 1
      ! map index to (mkq, nk, mu)
      bands = idx2band( tscat(i)%bands_index, kg%numb )
      m = bands(1);  n = bands(2);  im = bands(3)
      ! cycle if m != n for graphene only
      ! eigenvalue
      ! kelly find frequency of q using qg
      wq   = qg%enk(im, iq)
      ! kelly find frequency of negative q
      wnq  = qg%enk(im, inq)
      ! test that qs are correct
      ! test if wnq and wq are the same
      enk  = kg%enk(n, ik)
      emkq = kg%enk(m, ikq)
      !if ( abs(enk - emkq) < phfreq_cutoff ) cycle
      ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
      g2 = tscat(i)%g2
      ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
      if ( adapt_smear_eph ) then
         smear_tmp = tsmear_list(i)
      else
         smear_tmp = delta_smear
      endif
      ediff_abs = abs(enk - emkq + wq)
      ediff_ems = abs(enk - emkq - wnq)
      dscat = factor * g2
      ! electron distributaion
      fnk  = dist(n, ik)
      fmkq = dist(m, ikq)
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

      ! update nk with contribution from k+q
      mkq2nk = imph - imphn
      nk2mkq = - mkq2nk
      imph = imph * weight_ratio 
      imphn = imphn * weight_ratio 
      !update nk with contribution from k+q
!$omp atomic
      epcol(n, ik) = epcol(n, ik) + mkq2nk
      ! the inverse: update mk+q with contribution from k
!$omp atomic
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
   enddo
!$omp end do
!$omp end parallel
   ! combine contributions from all the pools/processes.
   call mp_sum(epcol, inter_pool_comm)
   call mp_sum(pecol, inter_pool_comm)
!   call mp_barrier(inter_pool_comm)
   call stop_clock('eph_scatint')
end subroutine cphdyna_scat_int

!subroutine ph3_scat_int(qg, disp, ph3col, scat_flag)
!   use boltz_utils, only: kpts_plus_invert
!   implicit none
!   type(grid), intent(in) :: qg
!   real(dp), intent(in) :: disp(qg%numb, qg%nk)
!   real(dp), intent(out) :: ph3col(qg%numb, qg%nk)
!   logical, intent(in), optional :: scat_flag
!   !local variables
!   logical :: scat_flag_loc
!   integer :: ivec, is, iq0, iq1, iq2, it, im0, im1, im2, bands(3)
!   real(dp) :: wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat
!   real(dp) :: occ_sum_0, occ_sum_1, occ_sum_2, factor
!   real(dp) :: qweight
!   call start_clock('phph_scatint')
!   if ( present(scat_flag) ) then
!      if ( scat_flag ) scat_flag_loc = scat_flag
!   endif
!   qweight = 1.0_dp / (boltz_qdim(1) * boltz_qdim(2) * boltz_qdim(3))
!   ! this constant means that the output is Im Simga only
!   ! to convert to tau^-1 needs to multiply results by 2/hbar
!   factor = 9.0E0_dp * 4.0E0_dp * pi * qweight 
!   ph3col(:,:) = 0.0E0_dp
!!$omp parallel do schedule(guided) default(shared) private(is, iq0, iq1, iq2, it, bands, im0, im1, im2, &
!!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat, occ_sum_0, occ_sum_1, occ_sum_2) &
!!$omp& reduction(+:ph3col)
!   do is = 1, num_qq_pair
!      iq0  = ph3_scat(is)%ik;   iq1  = ph3_scat(is)%iq
!      !iq0  = ph3_scat(is)%ik;   iq1  = ph3_scat(is)%iq;    iq2 = ph3_scat(is)%ikq
!      iq2 = kpts_plus_invert(qg%kpt(iq0), qg%kpt(iq1), qg%ndim) + 1
!      do it = 1, ph3_scat(is)%nchl
!         
!         !if (ph3_scat(is)%bands_index(it) .eq. 0) cycle
!         
!         bands = idx2band( ph3_scat(is)%bands_index(it), qg%numb )
!         im2 = bands(1); im0 = bands(2); im1 = bands(3)
!         
!         g2 = ph3_scat(is)%eph_g2(it)
!         
!        ! eigenvalue
!         wq0 = qg%enk(im0, iq0); wq1 = qg%enk(im1, iq1); wq2 = qg%enk(im2, iq2)
!        ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
!        !> d3q  add the factor of 18pi/hbar^2/36, 0.5^3 is included in g2, hbar is omitted in atomic units.
!        !> tdep add prefactor of hbar*pi/2, 1/8 is preincluded in g2
!        ! delta(wq + wq1 - wq2)
!         if ( adapt_smear_phph ) then
!            e_delta = gauss(tsmear_ph(is)%adsmear(it), -wq2 + wq1 + wq0)
!            !e_delta = gauss(ph3_scat(is)%adsmear(it), -wq2 + wq1 + wq0)
!         else
!            e_delta = gauss(delta_smear_ph, -wq2 + wq1 + wq0)
!         endif
!        
!         ! halving equal ik and iq
!         if ((im0 .eq. im1) .and. (iq0 .eq. iq1)) then
!            e_delta = e_delta/2.0_dp
!         endif
!        ! phonon occupation
!         ocq0 = disp(im0, iq0)
!         ocq1 = disp(im1, iq1)
!         ocq2 = disp(im2, iq2)
!
!         if ( scat_flag_loc ) then
!            ! scattering process:
!            dscat = g2 * e_delta * factor
!            occ_sum_0 =  dscat*(ocq1 - ocq2)
!            occ_sum_1 =  dscat*(ocq0 - ocq2)
!            occ_sum_2 =  dscat*(ocq0 + ocq1 + 1.0_dp)
!         else
!            ! dynamics - another factor of 2
!            occ_sum_0 = 2.0_dp * g2 * e_delta * factor * &
!               (ocq2 * (ocq0 + ocq1 + 1.0_dp) - ocq0 * ocq1)
!            occ_sum_1 = occ_sum_0
!            occ_sum_2 = - occ_sum_0
!         endif
!         ph3col(im0, iq0) = ph3col(im0, iq0) + occ_sum_0
!         ph3col(im1, iq1) = ph3col(im1, iq1) + occ_sum_1
!         !if (rgrp_idx2 .eq. rgrp_idx) then
!         ph3col(im2, iq2) = ph3col(im2, iq2) + occ_sum_2
!         !endif
!       enddo! it loop
!   enddo ! is loop
!!$omp end parallel do
!
!   call mp_sum(ph3col, inter_pool_comm)
!   !call mp_barrier(inter_pool_comm)
!   call stop_clock('phph_scatint')
!end subroutine ph3_scat_int


subroutine ph3_scat_int(qg, disp, ph3col, scat_flag)
   use boltz_utils, only: kpts_plus_invert
   implicit none
   type(grid), intent(in) :: qg
   real(dp), intent(in) :: disp(qg%numb, qg%nk)
   real(dp), intent(out) :: ph3col(qg%numb, qg%nk)
   logical, intent(in), optional :: scat_flag
   !local variables
   logical :: scat_flag_loc
   integer :: ivec, is, iq0, iq1, iq2, it, im0, im1, im2, bands(3)
   real(dp) :: wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat
   real(dp) :: occ_sum_0, occ_sum_1, occ_sum_2, factor
   real(dp) :: qweight
   integer :: scat_index
   real(dp), allocatable :: ph3col_tmp(:,:)
   call start_clock('phph_scatint')
   if ( present(scat_flag) ) then
      if ( scat_flag ) scat_flag_loc = scat_flag
   endif
   qweight = 1.0_dp / (boltz_qdim(1) * boltz_qdim(2) * boltz_qdim(3))
   ! this constant means that the output is Im Simga only
   ! to convert to tau^-1 needs to multiply results by 2/hbar
   factor = 9.0E0_dp * 4.0E0_dp * pi * qweight 
   ph3col(:,:) = 0.0E0_dp
   allocate(ph3col_tmp(qg%numb, qg%nk))
!$omp parallel default(shared) private(is, iq0, iq1, iq2, bands, im0, im1, im2, &
!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat, occ_sum_0, occ_sum_1, occ_sum_2, &
!$omp& scat_index, ph3col_tmp) 
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
      if ( adapt_smear_phph ) then
         e_delta = gauss(tsmear_ph_list(is), -wq2 + wq1 + wq0)
      else
         e_delta = gauss(delta_smear_ph, -wq2 + wq1 + wq0)
      endif
        
      ! halving equal ik and iq
      if ((im0 .eq. im1) .and. (iq0 .eq. iq1)) then
         e_delta = e_delta/2.0_dp
      endif
      ! phonon occupation
      ocq0 = disp(im0, iq0)
      ocq1 = disp(im1, iq1)
      ocq2 = disp(im2, iq2)

      if ( scat_flag_loc ) then
         ! scattering process:
         dscat = g2 * e_delta * factor
         occ_sum_0 =  dscat*(ocq1 - ocq2)
         occ_sum_1 =  dscat*(ocq0 - ocq2)
         occ_sum_2 =  dscat*(ocq0 + ocq1 + 1.0_dp)
      else
         ! dynamics - another factor of 2
         dscat = 2.0_dp * g2 * e_delta * factor
         occ_sum_0 = dscat * (ocq2 * (ocq0 + ocq1 + 1.0_dp) - ocq0 * ocq1)
         occ_sum_1 = occ_sum_0
         occ_sum_2 = - occ_sum_0
      endif
      ph3col_tmp(im0, iq0) = ph3col_tmp(im0, iq0) + occ_sum_0
      ph3col_tmp(im1, iq1) = ph3col_tmp(im1, iq1) + occ_sum_1
      ph3col_tmp(im2, iq2) = ph3col_tmp(im2, iq2) + occ_sum_2
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
   !call mp_barrier(inter_pool_comm)
   call stop_clock('phph_scatint')
end subroutine ph3_scat_int

!subroutine ph3_scat_int_interp2(qg, disp, ph3col, scat_flag)
!   use boltz_utils, only: kpts_plus_invert, num2kpt, kpt2num, fold_k
!   implicit none
!   type(grid), intent(in) :: qg
!   real(dp), intent(in) :: disp(qg%numb, qg%nk)
!   real(dp), intent(out) :: ph3col(qg%numb, qg%nk)
!   logical, intent(in), optional :: scat_flag
!   !local variables
!   logical :: scat_flag_loc
!   integer :: ivec, is, iq0, iq1, iq2, it, im0, im1, im2, bands(3)
!   real(dp) :: wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat
!   real(dp) :: occ_sum_0, occ_sum_1, occ_sum_2, factor
!   real(dp) :: qweight
!   real(dp) :: fweights(3), q0c(3), q1c(3), q0(3), q1(3)
!   integer :: xinc, yinc, zinc, iq0c, iq1c
!   integer :: xinc1, yinc1, zinc1
!   call start_clock('phph_scatint')
!
!   fweights(1) = 1.0_dp/qg%ndim(1)
!   fweights(2) = 1.0_dp/qg%ndim(2)
!   fweights(3) = 1.0_dp/qg%ndim(3)
!   if ( present(scat_flag) ) then
!      if ( scat_flag ) scat_flag_loc = scat_flag
!   endif
!   qweight = 1.0_dp / (boltz_qdim(1) * boltz_qdim(2) * boltz_qdim(3))
!   !qweight = 1.0_dp / (qg_fine%ndim(1) * qg_fine%ndim(2) * qg_fine%ndim(3))
!   ! this constant means that the output is Im Simga only
!   ! to convert to tau^-1 needs to multiply results by 2/hbar
!   factor = 9.0E0_dp * 4.0E0_dp * pi * qweight 
!   ph3col(:,:) = 0.0E0_dp
!!$omp parallel do schedule(guided) default(shared) private(is, iq0, iq1, iq2, it, bands, im0, im1, im2, &
!!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat, occ_sum_0, occ_sum_1, occ_sum_2, &
!!$omp& iq0c, iq1c, q0c, q1c, q0, q1, xinc, yinc, zinc, xinc1, yinc1, zinc1)
!   do is = 1, num_qq_pair
!      iq0c = ph3_scat(is)%ik;   iq1c = ph3_scat(is)%iq!;   iq2 = ph3_scat(is)%ikq
!      ! find the corresponding iqfs
!      do it = 1, ph3_scat(is)%nchl
!         ! do a loop for iqf
!         ! band index is defined the same, here we use qg, but qg_fine%numb is the same 
!         bands = idx2band( ph3_scat(is)%bands_index(it), qg%numb )
!         im2 = bands(1); im0 = bands(2); im1 = bands(3)
!         ! key step, to keep g2 the same 
!         g2 = ph3_scat(is)%eph_g2(it)
!         if ( g2 == 0.0_dp ) cycle
!         q0c = num2kpt(qg%kpt(iq0c), qg%ndim)
!         q1c = num2kpt(qg%kpt(iq1c), qg%ndim)
!         !write(*,*) 'original= ', iq0, iq1, iq2
!         !write(*,*) 'num of points =', size_iqf+1
!         do xinc = 1, 2 
!         do yinc = 1, 2 
!         do zinc = 1, 2 
!         do xinc1 = 1, 2 
!         do yinc1 = 1, 2 
!         do zinc1 = 1, 2 
!            q0(1) = fold_k(q0c(1) + (xinc - 1) * fweights(1))
!            q0(2) = fold_k(q0c(2) + (yinc - 1) * fweights(2))
!            q0(3) = fold_k(q0c(3) + (zinc - 1) * fweights(3))
!            q1(1) = fold_k(q1c(1) + (xinc1 - 1) * fweights(1))
!            q1(2) = fold_k(q1c(2) + (yinc1 - 1) * fweights(2))
!            q1(3) = fold_k(q1c(3) + (zinc1 - 1) * fweights(3))
!            iq0 = kpt2num(q0, qg%ndim) + 1
!            iq1 = kpt2num(q1, qg%ndim) + 1
!            if (iq0 < 0) cycle
!            if (iq0 > qg%nk) cycle
!            if (iq1 < 0) cycle
!            if (iq1 > qg%nk) cycle
!            !iq2 = kpt2num(-q0-q1, qg%ndim) + 1
!            iq2 = kpts_plus_invert(qg%kpt(iq0), qg%kpt(iq1), qg%ndim) + 1
!            if (iq2 < 0) cycle
!            if (iq2 > qg%nk) cycle
!            !write(*,*) 'iqf_id = ', iq0, iq1, iq2
!         
!            ! eigenvalue
!            wq0 = qg%enk(im0, iq0); wq1 = qg%enk(im1, iq1); wq2 = qg%enk(im2, iq2)
!            !wq0 = qg%enk(im0, iq0); wq1 = qg%enk(im1, iq1); wq2 = qg%enk(im2, iq2)
!            ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
!            !> d3q  add the factor of 18pi/hbar^2/36, 0.5^3 is included in g2, hbar is omitted in atomic units.
!            !> tdep add prefactor of hbar*pi/2, 1/8 is preincluded in g2
!            ! delta(wq + wq1 - wq2)
!            if ( adapt_smear_phph ) then
!               e_delta = gauss(tsmear_ph(is)%adsmear(it), -wq2 + wq1 + wq0)
!               !e_delta = gauss(ph3_scat(is)%adsmear(it), -wq2 + wq1 + wq0)
!            else
!               e_delta = gauss(delta_smear_ph, -wq2 + wq1 + wq0)
!            endif
!        
!            ! halving equal ik and iq
!            if ((im0 .eq. im1) .and. (iq0 .eq. iq1)) then
!               e_delta = e_delta/2.0_dp
!            endif
!            ! phonon occupation
!            ocq0 = disp(im0, iq0)
!            ocq1 = disp(im1, iq1)
!            ocq2 = disp(im2, iq2)
!
!            if ( scat_flag_loc ) then
!               ! scattering process:
!               dscat = g2 * e_delta * factor
!               occ_sum_0 =  dscat*(ocq1 - ocq2)
!               occ_sum_1 =  dscat*(ocq0 - ocq2)
!               occ_sum_2 =  dscat*(ocq0 + ocq1 + 1.0_dp)
!            else
!               ! dynamics - another factor of 2
!               occ_sum_0 = 2.0_dp * g2 * e_delta * factor * &
!                  (ocq2 * (ocq0 + ocq1 + 1.0_dp) - ocq0 * ocq1)
!               occ_sum_1 = occ_sum_0
!               occ_sum_2 = - occ_sum_0
!            endif
!!$omp atomic
!            ph3col(im0, iq0) = ph3col(im0, iq0) + occ_sum_0
!!$omp atomic
!            ph3col(im1, iq1) = ph3col(im1, iq1) + occ_sum_1
!!$omp atomic
!            ph3col(im2, iq2) = ph3col(im2, iq2) + occ_sum_2
!         enddo; enddo; enddo! iqf loop
!         enddo; enddo; enddo! iqf loop
!      enddo! it loop
!   enddo ! is loop
!!$omp end parallel do
!
!   call mp_sum(ph3col, inter_pool_comm)
!   !call mp_barrier(inter_pool_comm)
!   call stop_clock('phph_scatint')
!end subroutine ph3_scat_int_interp2

subroutine ph3_scat_int_interp2(qg, disp, ph3col, scat_flag)
   use boltz_utils, only: kpts_plus_invert, num2kpt, kpt2num, fold_k, &
      num2ipt
   use hdf5_utils
   use pert_utils, only : find_free_unit
   implicit none
   type(grid), intent(in) :: qg
   real(dp), intent(in) :: disp(qg%numb, qg%nk)
   real(dp), intent(out) :: ph3col(qg%numb, qg%nk)
   logical, intent(in), optional :: scat_flag
   !local variables
   logical :: scat_flag_loc
   integer :: ivec, is, iq0, iq1, iq2, it, im0, im1, im2, bands(3)
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
   
   real(dp) :: dscat_save(729)
   real(dp) :: mean
   real(dp) :: averagec
   integer(HID_T) :: file_id, iout
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
!$omp parallel default(shared) private(is, iq0, iq1, iq2, bands, im0, im1, im2, &
!$omp& wq0, wq1, wq2, g2, e_delta, ocq0, ocq1, ocq2, dscat, occ_sum_0, occ_sum_1, occ_sum_2, &
!$omp& iq0c, iq1c, ipt0c, ipt1c, ipt0, ipt1, xinc, yinc, zinc, xinc1, yinc1, zinc1, scat_index, &
!$omp& pos_vec, ph3col_tmp, wdiff) 
   allocate(ph3col_tmp(qg%numb, qg%nk))
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
   deallocate(ph3col_tmp)
!$omp end parallel

   call mp_sum(ph3col, inter_pool_comm)
   call stop_clock('phph_scatint')
end subroutine ph3_scat_int_interp2


end module boltz_scatter_integral
