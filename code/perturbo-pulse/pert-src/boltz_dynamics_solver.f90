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
!  algorithms to solve the full BTE equation
!
! Maintenance:
!===============================================================================

!> Real-time solvers for the Boltzmann transport equation.
module boltz_dynamics_solver
   use pert_const, only: dp
   use boltz_grid, only: grid
   use boltz_scatter_integral, only: cdyna_scat_int, cphdyna_scat_int, ph3_scat_int, ph3_scat_int_interp2, &
      ph3_scat_int_interp_half
   use boltz_scatter_integral_tgt, only: cdyna_scat_int_tgt
   use boltz_external_field, only: calc_drift_efield
   use pert_param, only: boltz_efield, scat_impl, divide_grid 
   implicit none
   private
   
   public :: euler
   public :: runge_kutta_4th
   public :: euler_eph, euler_phph
   public :: runge_kutta_4th_eph, runge_kutta_4th_phph
contains


subroutine euler(kgrid, tstep, cdist, ndist, dertot, der_k, ept, calc_adv)
   implicit none
   type(grid), intent(in) :: kgrid
   ! time step, temperature for e-ph, and electric field in atomic unit
   real(dp), intent(in)  :: tstep, ept
   ! cdist: f( t_i ); ndist: f( t_(i+1) )
   real(dp), intent(in)  :: cdist(kgrid%numb, kgrid%nk)
   real(dp), intent(out) :: ndist(kgrid%numb, kgrid%nk)
   real(dp), intent(out) :: dertot(kgrid%numb, kgrid%nk)
   real(dp), allocatable, intent(inout) :: der_k(:,:) ! this array is allocated if calc_adv
   logical, intent(in)   :: calc_adv
   
   ! compute electron-phonon collision term and (if calc_adv_cur) advection term
   call calc_boltz_dertot(kgrid, ept, cdist, der_k, dertot, calc_adv)

   !update 
   ndist(:,:) = cdist(:,:) + dertot(:,:) * tstep

end subroutine euler


subroutine euler_eph(kgrid, qgrid, tstep, cdist, cndist, dertot, &
   phdist, phndist, phdertot)

   implicit none
   type(grid), intent(in) :: kgrid, qgrid
   ! time step, temperature for e-ph, and electric field in atomic unit
   real(dp), intent(in)  :: tstep
   ! cdist: f( t_i ); cndist: f( t_(i+1) ); phdist: n(t_i); phndist: n( t_(i+1) )
   real(dp), intent(in)  :: cdist(kgrid%numb, kgrid%nk)
   real(dp), intent(inout) :: phdist(qgrid%numb, qgrid%nk)
   real(dp), intent(out) :: cndist(kgrid%numb, kgrid%nk), phndist(qgrid%numb, qgrid%nk)
   real(dp), intent(out) :: dertot(kgrid%numb, kgrid%nk), phdertot(qgrid%numb, qgrid%nk)
   
   call start_clock('euler_eph')   
   dertot = 0.0E0_dp; phdertot(:,:) = 0.0E0_dp
   call cphdyna_scat_int(kgrid, qgrid, cdist, dertot, phdist, phdertot)
   !update
   cndist(:,:) = cdist(:,:) + dertot(:,:) * tstep
   phndist(:,:) = phdist(:,:) + phdertot(:,:) * tstep
   call stop_clock('euler_eph')
end subroutine euler_eph

subroutine euler_phph(qgrid, tstep, &
   phdist, phndist, phdertot, istep, phphstep)

   implicit none
   type(grid), intent(in) :: qgrid
   real(dp), intent(in) :: tstep
   real(dp), intent(inout) :: phdist(qgrid%numb, qgrid%nk)
   real(dp), intent(inout) :: phndist(qgrid%numb, qgrid%nk)
   real(dp), intent(inout) :: phdertot(qgrid%numb, qgrid%nk)
   integer, intent(in) :: istep, phphstep ! only update ph-ph collision term
   real(dp) :: newtstep
   call start_clock('euler_phph')
   newtstep = tstep * phphstep
   phdist(:,:) = phndist(:,:) ! crucial we calculate ph-ph based on the new
   !phonon distribution calculated from e-ph

   ! only update phdertot for every phphsteps
   if( mod(istep-1, phphstep) .eq. 0) then
      if (divide_grid) then
         !call ph3_scat_int_interp2(qgrid, phdist, phdertot)
         call ph3_scat_int_interp_half(qgrid, phdist, phdertot)
      else
         call ph3_scat_int(qgrid, phdist, phdertot)
      endif
   endif
   phndist(:,:) = phndist(:,:) + phdertot(:,:) * newtstep
   call stop_clock('euler_phph')
end subroutine euler_phph


subroutine runge_kutta_4th(kgrid, tstep, cdist, ndist, df_tmp, dertot, der_k, ept, calc_adv)
   implicit none
   type(grid), intent(in) :: kgrid
   ! time step, temperature for e-ph, and electric field in atomic unit
   real(dp), intent(in)  :: tstep, ept
   ! cdist: f( t_i ); ndist: f( t_(i+1) )
   real(dp), intent(in)  :: cdist(kgrid%numb, kgrid%nk)
   real(dp), intent(out) :: ndist(kgrid%numb, kgrid%nk)
   ! workspace, to avoid frequently allocate and deallocate large arrays,
   !   which could hurt efficiency and parallel performance severely.
   real(dp), intent(out) :: df_tmp(kgrid%numb, kgrid%nk), &
                            dertot(kgrid%numb, kgrid%nk)
   real(dp), allocatable, intent(inout) :: der_k(:,:) ! this array is allocated if calc_adv
   ! local variables
   logical, intent(in)   :: calc_adv
   real(dp) :: half_t, t_six, t_three
   
   half_t = tstep*0.5E0_dp;  t_three = tstep/3.E0_dp;  t_six = tstep/6.E0_dp
   !
   df_tmp = 0.0E0_dp; dertot = 0.0E0_dp
   ! compute k1
   ! compute electron-phonon collision term and (if calc_adv_cur) advection term
   call calc_boltz_dertot(kgrid, ept, cdist, der_k, dertot, calc_adv)
   ! collect contribution from k1
   ndist(:,:) = cdist(:,:) + dertot(:,:) * t_six

   ! compute k2
   df_tmp(:,:) = cdist(:,:) + half_t * dertot(:,:)
   call calc_boltz_dertot(kgrid, ept, df_tmp, der_k, dertot, calc_adv)
   ! collect contribution from k2
   ndist(:,:) = ndist(:,:) + dertot(:,:) * t_three

   ! compute k3
   df_tmp(:,:) = cdist(:,:) + half_t * dertot(:,:)
   call calc_boltz_dertot(kgrid, ept, df_tmp, der_k, dertot, calc_adv)
   ! collect contribution from k3
   ndist(:,:) = ndist(:,:) + dertot(:,:) * t_three

   ! compute k4
   df_tmp(:,:) = cdist(:,:) + tstep * dertot(:,:)
   call calc_boltz_dertot(kgrid, ept, df_tmp, der_k, dertot, calc_adv)
   ! collect contribution from k4
   ndist(:,:) = ndist(:,:) + dertot(:,:) * t_six

end subroutine runge_kutta_4th


subroutine runge_kutta_4th_eph(kgrid, qgrid, tstep, cdist, cndist, &
                df_tmp, dertot, phdist, phndist, dn_tmp, phdertot)

   implicit none
   type(grid), intent(in) :: kgrid, qgrid
   ! time step for e-ph, and electric field in atomic unit
   real(dp), intent(in)  :: tstep
   ! cdist: f( t_i ); cndist: f( t_(i+1) ); phdist: n(t_i); phndist: n( t_(i+1) )
   real(dp), intent(in)  :: cdist(kgrid%numb, kgrid%nk)
   real(dp), intent(inout) :: phdist(qgrid%numb, qgrid%nk)
   real(dp), intent(out) :: cndist(kgrid%numb, kgrid%nk), phndist(qgrid%numb, qgrid%nk)
   ! workspace, to avoid frequently allocate and deallocate large arrays,
   !   which could hurt efficiency and parallel performance severely.
   real(dp), intent(out) :: df_tmp(kgrid%numb, kgrid%nk), dertot(kgrid%numb, kgrid%nk), &
                            dn_tmp(qgrid%numb, qgrid%nk), phdertot(qgrid%numb, qgrid%nk)
   ! local variables
   real(dp) :: half_t, t_six, t_three
   
   ! record time
   call start_clock('rk4_eph')
   
   half_t = tstep*0.5E0_dp;  t_three = tstep/3.E0_dp;  t_six = tstep/6.E0_dp
   
   !- e-ph and ph-e considered at the same time
   df_tmp = 0.0E0_dp; dertot = 0.0E0_dp; dn_tmp = 0.0E0_dp; phdertot = 0.0E0_dp
   ! compute k1
   ! compute electron-phonon collision term
   call cphdyna_scat_int(kgrid, qgrid, cdist, dertot, phdist, phdertot)
   ! collect contribution from k1
   cndist(:,:) = cdist(:,:) + dertot(:,:) * t_six
   phndist(:,:) = phdist(:,:) + phdertot(:,:) * t_six

   
   ! compute k2
   df_tmp(:,:) = cdist(:,:) + half_t * dertot(:,:)
   dn_tmp(:,:) = phdist(:,:) + half_t * phdertot(:,:)
   call cphdyna_scat_int(kgrid, qgrid, df_tmp, dertot, dn_tmp, phdertot)
   ! collect contribution from k2
   cndist(:,:) = cndist(:,:) + dertot(:,:) * t_three
   phndist(:,:) = phndist(:,:) + phdertot(:,:) * t_three

   ! compute k3
   df_tmp(:,:) = cdist(:,:) + half_t * dertot(:,:)
   dn_tmp(:,:) = phdist(:,:) + half_t * phdertot(:,:)
   call cphdyna_scat_int(kgrid, qgrid, df_tmp, dertot, dn_tmp, phdertot)
   ! collect contribution from k3
   cndist(:,:) = cndist(:,:) + dertot(:,:) * t_three
   phndist(:,:) = phndist(:,:) + phdertot(:,:) * t_three

   ! compute k4
   df_tmp(:,:) = cdist(:,:) + tstep * dertot(:,:)
   dn_tmp(:,:) = phdist(:,:) + tstep * phdertot(:,:)
   call cphdyna_scat_int(kgrid, qgrid, df_tmp, dertot, dn_tmp, phdertot)
   ! collect contribution from k4
   cndist(:,:) = cndist(:,:) + dertot(:,:) * t_six
   phndist(:,:) = phndist(:,:) + phdertot(:,:) * t_six

   call stop_clock('rk4_eph')

end subroutine runge_kutta_4th_eph

subroutine runge_kutta_4th_phph(qgrid, tstep, &
                phdist, phndist, dn_tmp, phdertot, istep, phphstep)

   implicit none
   type(grid), intent(in) :: qgrid
   ! time step, temperature for e-ph, and electric field in atomic unit
   real(dp), intent(in)  :: tstep
   real(dp) :: newtstep
   ! cdist: f( t_i ); cndist: f( t_(i+1) ); phdist: n(t_i); phndist: n( t_(i+1) )
   real(dp), intent(inout) :: phdist(qgrid%numb, qgrid%nk)
   real(dp), intent(inout) :: phndist(qgrid%numb, qgrid%nk) ! has to be inout

   ! workspace, to avoid frequently allocate and deallocate large arrays,
   !   which could hurt efficiency and parallel performance severely.
   real(dp), intent(inout) :: dn_tmp(qgrid%numb, qgrid%nk), phdertot(qgrid%numb, qgrid%nk)
   real(dp) :: half_t, t_six, t_three
   
   integer, intent(in) :: istep, phphstep ! only update ph-ph collision term
   !every phphstep 
   newtstep = tstep * phphstep
   ! record time
   call start_clock('rk4_phph')
   
   phdist(:,:) = phndist(:,:) ! update the original distribution
   
   if( mod(istep-1, phphstep) .eq. 0) then   
      half_t = newtstep*0.5E0_dp;  t_three = newtstep/3.E0_dp;  t_six = newtstep/6.E0_dp
      
      !--ph-ph interaction--
      dn_tmp = 0.0E0_dp; phdertot=0.0E0_dp
      if (divide_grid) then
         !call ph3_scat_int_interp2(qgrid, phdist, phdertot)
         call ph3_scat_int_interp_half(qgrid, phdist, phdertot)

      else
         call ph3_scat_int(qgrid, phdist, phdertot)
      endif
      ! collect contribution from k1
      phndist(:,:) = phndist(:,:) + phdertot(:,:) * t_six
      !call errore('Good afternoon ', 'We stopped it intentionally, no worries!', 1)
   
   
      ! compute k2
      dn_tmp(:,:) = phdist(:,:) + half_t * phdertot(:,:)
      if (divide_grid) then
         !call ph3_scat_int_interp2(qgrid, dn_tmp, phdertot)
         call ph3_scat_int_interp_half(qgrid, phdist, phdertot)
      else
         call ph3_scat_int(qgrid, dn_tmp, phdertot)
      endif
      ! collect contribution from k2
      phndist(:,:) = phndist(:,:) + phdertot(:,:) * t_three
   
      ! compute k3
      dn_tmp(:,:) = phdist(:,:) + half_t * phdertot(:,:)
      if (divide_grid) then
         !call ph3_scat_int_interp2(qgrid, dn_tmp, phdertot)
         call ph3_scat_int_interp_half(qgrid, phdist, phdertot)
      else
         call ph3_scat_int(qgrid, dn_tmp, phdertot)
      endif
      ! collect contribution from k3
      phndist(:,:) = phndist(:,:) + phdertot(:,:) * t_three
   
      ! compute k4
      dn_tmp(:,:) = phdist(:,:) + newtstep * phdertot(:,:)
      if (divide_grid) then
         !call ph3_scat_int_interp2(qgrid, dn_tmp, phdertot)
         call ph3_scat_int_interp_half(qgrid, phdist, phdertot)
      else
         call ph3_scat_int(qgrid, dn_tmp, phdertot)
      endif
      ! collect contribution from k4
      phndist(:,:) = phndist(:,:) + phdertot(:,:) * t_six
   
   !else ! keep phdertot as input
   !   phndist(:,:) = phndist(:,:) + phdertot(:,:) * tstep

   endif

   call stop_clock('rk4_phph')

end subroutine runge_kutta_4th_phph

subroutine calc_boltz_dertot(kgrid, ept, dist, der_k, dertot, calc_adv)
   implicit none
   type(grid), intent(in) :: kgrid
   real(dp), intent(in)   :: ept
   real(dp), intent(in)  :: dist(kgrid%numb, kgrid%nk)
   real(dp), allocatable, intent(inout) :: der_k(:,:) ! this array is allocated if calc_adv
   real(dp), intent(out) :: dertot(kgrid%numb, kgrid%nk)
   logical, intent(in)   :: calc_adv
   ! local variables
   integer :: iband
   real(dp), allocatable :: dist_tmp(:,:)

   dertot = 0.0E0_dp

   call start_clock('dynamics_col.')

   if (trim(scat_impl) .eq. 'std') then
      call cdyna_scat_int(kgrid, ept, dist, dertot)
   else if (trim(scat_impl) .eq. 'tgt') then
      call cdyna_scat_int_tgt(kgrid, ept, dist, dertot)
   else
      call errore('carrier dynamics', 'unrecognized scat_impl value '// trim(scat_impl), 1)
   endif

   call stop_clock('dynamics_col.')

   if (calc_adv) then

      call start_clock('dynamics_adv.')

      call calc_drift_efield(kgrid, boltz_efield, dist, der_k)

      dertot(:,:) = dertot(:,:) + der_k(:,:)

      call stop_clock('dynamics_adv.')
   end if

end subroutine calc_boltz_dertot

end module boltz_dynamics_solver
