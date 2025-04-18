!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
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

module sundials
#if defined (__SUNDIALS)
   use, intrinsic :: iso_c_binding
   use kinds, only: dp
   use boltz_grid, only: grid
   use boltz_scatter_integral, only: cdyna_scat_int, cphdyna_scat_int, ph3_scat_int, &
      ph3_scat_int_interp_half
   use qe_mpi_mod, only: ionode, stdout
   use pert_param, only: divide_grid
   
   use fsundials_core_mod

   use farkode_mod
   use farkode_mristep_mod
   use farkode_arkstep_mod
   use farkode_erkstep_mod 
   use fnvector_serial_mod

   implicit none

   type(grid), save :: kgs, qgs
   real(dp), save :: temps
   integer(c_long), save :: knb, knk, qnb, qnk
   ! kelly changed save
   real(dp), allocatable :: cdist(:,:), phdist(:,:)
   real(dp), allocatable :: epcol(:,:), pecol(:,:), ph3col(:,:), pecol_sum(:,:)
   private :: create_rk4_butcher
   private :: f, ff, fs, f0, error_fun
   
   public :: vec_pert_to_sun, vec_sun_to_pert
   public :: init_sundials
   public :: fill_user_data
   public :: set_integrators_wrapper, sundials_stats_wrapper
   public :: sundials_evolve
   public :: sundials_reset
   public :: deallocate_sunarrays, free_sundials_memory
contains

subroutine vec_pert_to_sun(cdist, sundist, phdist)
!!! This subroutine reshapes carrier distribution and phonon distribution into one N-Vector
   !!! [(n1, k1), (n2, k1), ... , (n_nb, k1), (n1, k2) ..., (n_nb, k_nk),
   !!! (v1, q1), ... , (v_nm, q_nq) ]

   !!! one can also use this subroutine to reshape carrier and phonon collision
   !!! integrals (derivatives), since the dimension matches
   use, intrinsic :: iso_c_binding
   implicit none
   real(dp), intent(in) :: cdist(knb, knk)
   real(dp), intent(in), optional :: phdist(qnb, qnk)
   real(c_double), pointer, intent(inout) :: sundist(:)
   integer :: ik,ib,indx !counter
   call start_clock('pts') 
   indx=1
   do ik = 1,knk
      do ib = 1,knb
         sundist(indx) = cdist(ib,ik)
         indx=indx+1
      enddo
   enddo
   if ( present(phdist) ) then
      do ik = 1,qnk
         do ib = 1,qnb
            sundist(indx) = phdist(ib,ik)
            indx=indx+1
         enddo
      enddo
   endif
   call stop_clock('pts')

end subroutine vec_pert_to_sun

subroutine vec_sun_to_pert(cdist, sundist, phdist)
!
!   !!! This subroutine reshapes carrier distribution and phonon distribution from sundials to perturbo
   use, intrinsic :: iso_c_binding
   implicit none
   real(dp), intent(inout) :: cdist(kgs%numb, kgs%nk)
   real(dp), intent(inout), optional :: phdist(qgs%numb, qgs%nk)
   real(c_double), pointer, intent(in) :: sundist(:)
   integer :: ik,ib,indx !counter

   
   indx=1
   do ik = 1,knk
      do ib = 1,knb
         cdist(ib,ik) = sundist(indx)
         indx=indx+1
      enddo
   enddo

   if ( present(phdist) ) then
      do ik = 1,qnk
         do ib = 1,qnb
            phdist(ib,ik) = sundist(indx)
            indx=indx+1
         enddo
      enddo
   endif
end subroutine vec_sun_to_pert


subroutine fill_user_data(kg, eptemp, qg)
   !!! put kg, qg, eptemp into udata
   use, intrinsic :: iso_c_binding

   implicit none

   type(grid), intent(in) :: kg
   type(grid), intent(in), optional :: qg
   real(dp), intent(in) :: eptemp
   ! store user data as local static variables
   kgs=kg
   temps = eptemp
   knb = kg%numb
   knk = kg%nk
   qnb = 0
   qnk = 0
   ! allocate cdist and the derivatives
   allocate(cdist(knb,knk))
   allocate(epcol(knb,knk))
   if (present(qg) ) then
      qgs=qg
      qnb = qg%numb
      qnk = qg%nk
      allocate(phdist(qnb,qnk))
      allocate(pecol(qnb,qnk), ph3col(qnb, qnk), pecol_sum(qnb,qnk))
   endif
end subroutine fill_user_data

subroutine init_sundials(cdist0, svec_dist, sundist, ctx, phdist0)
   !!! this function initialize the vectors for the calculation
   use, intrinsic :: iso_c_binding
   
   implicit none
   real(dp), target, intent(in) :: cdist0(knb,knk)
   real(dp), target, intent(in), optional :: phdist0(qnb, qnk)
   type(N_Vector), pointer, intent(out) :: svec_dist
   real(c_double), pointer, intent(out) :: sundist(:)
   type(c_ptr), intent(out) :: ctx
   integer(c_int) :: retval
   integer(c_long) :: neq

   !calculate number of coefficients for population
   neq=knb*knk+qnb*qnk
   ! create sunContext
   retval = FSunContext_Create(SUN_COMM_NULL, ctx)
   if (retval /= 0) then
      call errore('sundials', 'Error in FSunContext_Create', 1)
   endif 
   svec_dist => FN_VNew_Serial(neq, ctx)
   if ( .not.associated(svec_dist)) then 
      call errore('sundials', 'Error creating the Nvector for data', 1)
   endif
   sundist => FN_VGetArrayPointer(svec_dist)
   if ( present(phdist0) ) then
      call vec_pert_to_sun(cdist0, sundist, phdist0)
   else
      call vec_pert_to_sun(cdist0, sundist)
   endif
   

end subroutine init_sundials


integer(c_int) function f(t, svec_dist, svec_dot, udata) result(retval) bind(C)
   !!! takes in sundials inputs for cphdyna_scat_int
   !!!  and output f for sundials
   use, intrinsic :: iso_c_binding

   implicit none
   
   ! user data
   type(c_ptr) :: udata

   real(c_double), value :: t ! the value is important here
   type(N_Vector) :: svec_dist
   type(N_Vector) :: svec_dot
   real(c_double), pointer :: sundist(:)
   real(c_double), pointer :: sundist_dot(:)

   ! associate pointers with nvectors for dist and derivative
   sundist => FN_VGetArrayPointer(svec_dist)
   sundist_dot => FN_VGetArrayPointer(svec_dot)
   
   ! convert sundist to cdist and phdist
   if (qnb > 0) then
      call vec_sun_to_pert(cdist, sundist, phdist)
   else
      call vec_sun_to_pert(cdist, sundist)
   endif
   ph3col(:,:) = 0.0_dp
   if (qnb > 0) then
      ! get derivatives from the collision integrals
      call cphdyna_scat_int(kgs, qgs, cdist, epcol, phdist, pecol)
      ! assign derivatives to the pointer associated with the sundial derivative vector
      if (divide_grid) then
         call ph3_scat_int_interp_half(qgs, phdist, ph3col)
         !ph3col = 0.0_dp
      else
         call ph3_scat_int(qgs, phdist, ph3col)
      endif
      pecol_sum(:,:) = ph3col(:,:) + pecol(:,:)
   else
      call cdyna_scat_int(kgs, temps, cdist, epcol)
   endif
   if (qnb > 0) then
      call vec_pert_to_sun(epcol, sundist_dot, pecol_sum)
   else
      call vec_pert_to_sun(epcol, sundist_dot)
   endif

   retval = 0

   return
end function f

integer(c_int) function ff(t, svec_dist, svec_dot, udata) result(retval) bind(C)
   !!! takes in sundials inputs for cphdyna_scat_int
   !!!  and output ff for sundials
   use, intrinsic :: iso_c_binding

   implicit none
   
   ! user data
   type(c_ptr) :: udata

   real(c_double), value :: t ! the value is important here
   type(N_Vector) :: svec_dist
   type(N_Vector) :: svec_dot
   real(c_double), pointer :: sundist(:)
   real(c_double), pointer :: sundist_dot(:)

   ! associate pointers with nvectors for dist and derivative
   sundist => FN_VGetArrayPointer(svec_dist)
   sundist_dot => FN_VGetArrayPointer(svec_dot)
   
   ! convert sundist to cdist and phdist
   if (qnb > 0) then
      call vec_sun_to_pert(cdist, sundist, phdist)
   else
      call vec_sun_to_pert(cdist, sundist)
   endif
   if (qnb > 0) then
      ! get derivatives from the collision integrals
      call cphdyna_scat_int(kgs, qgs, cdist, epcol, phdist, pecol)
      ! assign derivatives to the pointer associated with the sundial derivative vector
      !call ph3_scat_int(qgs, phdist, ph3col)
      !ph3col(:,:) = 0.0_dp
      !pecol_sum(:,:) = ph3col(:,:) + pecol(:,:)
   else
      call cdyna_scat_int(kgs, temps, cdist, epcol)
   endif
   if (qnb > 0) then
      ! assign derivatives to the pointer associated with the sundial derivative vector
      !call vec_pert_to_sun(epcol, sundist_dot, pecol_sum)
      call vec_pert_to_sun(epcol, sundist_dot, pecol)
   else
      call vec_pert_to_sun(epcol, sundist_dot)
   endif

   retval = 0

   return
end function ff


integer(c_int) function fs(t, svec_dist, svec_dot, udata) result(retval) bind(C)
   !!! takes in sundials inputs for ph3_scat_int
   !!!  and output fs for sundials
   use, intrinsic :: iso_c_binding

   implicit none
   
   ! user data
   type(c_ptr) :: udata

   real(c_double), value :: t ! the value is important here
   type(N_Vector) :: svec_dist
   type(N_Vector) :: svec_dot
   real(c_double), pointer :: sundist(:)
   real(c_double), pointer :: sundist_dot(:)
    
   ! associate pointers with nvectors for dist and derivative
   sundist => FN_VGetArrayPointer(svec_dist)
   sundist_dot => FN_VGetArrayPointer(svec_dot)
   
   ! convert sundist to cdist and phdist
   if (qnb > 0) then
      call vec_sun_to_pert(cdist, sundist, phdist)
   else
      call vec_sun_to_pert(cdist, sundist)
   endif
   ph3col(:,:) = 0.0_dp
   ! get derivatives from the collision integrals
   if (qnb > 0) then
      if (divide_grid) then
         call ph3_scat_int_interp_half(qgs, phdist, ph3col)
      else
         call ph3_scat_int(qgs, phdist, ph3col)
      endif
   endif

   ! no change in electron population for slow time step
   epcol(:,:) = 0.0_dp
   if (qnb > 0) then
      call vec_pert_to_sun(epcol, sundist_dot, ph3col)
   else
      call vec_pert_to_sun(epcol, sundist_dot)
   endif

   retval = 0

   return
end function fs


integer(c_int) function f0(t, svec_dist, svec_dot, udata) result(retval) bind(C)
   use, intrinsic :: iso_c_binding
   ! user data
   type(c_ptr) :: udata

   real(c_double), value :: t ! the value is important here
   type(N_Vector) :: svec_dist
   type(N_Vector) :: svec_dot
   
   call FN_VConst(0.0d0, svec_dot)

   retval = 0
   return
end function f0


integer(c_int) function error_fun(svec_dist, ewt, udata) result(retval) bind(C)
   !! user supplied error weights for vector y

   !! error is 
   use, intrinsic :: iso_c_binding
   implicit none

   type(c_ptr) :: udata
   type(N_Vector) :: svec_dist
   type(N_Vector) :: ewt
   real(c_double), pointer :: y(:)
   real(c_double), pointer :: yewt(:)
   integer(c_int) :: ylength
    
   ! associate pointers with nvectors
   y => FN_VGetArrayPointer(svec_dist)
   ylength = FN_VGetLength(svec_dist)
   yewt => FN_VGetArrayPointer(ewt)
   ! error weights are just of retol * sqrt(N/norm(y)), where N = length(y)
   yewt(:) = sqrt(ylength/norm2(y))*1.0d5
   !yewt(:) = 1.0d-4
   retval = 0
end function error_fun

subroutine set_tol_vec(abstol_c, abstol_ph, ctx, abvec)
   implicit none
   real(c_double), intent(in) :: abstol_c, abstol_ph
   type(N_Vector), pointer, intent(out) :: abvec
   type(c_ptr), intent(inout) :: ctx
   real(c_double), pointer :: abvec_val(:)
   integer(c_long) :: neq, neqc, neqph 
   neqc = knb * knk
   neqph = qnb * qnk
   neq = neqc + neqph
   abvec => FN_VNew_Serial(neq,ctx)
   abvec_val => FN_VGetArrayPointer(abvec)
   abvec_val(1:neqc) = abstol_c
   if ( neqph .gt. 0 ) abvec_val((neqc+1):neq) = abstol_ph

end subroutine set_tol_vec

subroutine create_rk4_butcher(BTf)
   implicit none
   
   !butcher tables
   real(c_double), allocatable :: Af(:,:), bf(:), cf(:), df(:)
   type(c_ptr), intent(out) :: BTf
   real(c_double) :: ONE = 1.0d0

   ! create rk4 butcher table
   allocate(Af(4,4))
   allocate(bf(4))
   allocate(cf(4))
   allocate(df(4))
   Af = 0.d0
   bf = 0.d0
   cf = 0.d0
   df = 0.d0
   Af(1,2) = 0.5d0
   Af(2,3) = 0.5d0
   Af(3,4) = ONE
   bf(1) = ONE/6.0d0
   bf(2) = ONE/3.0d0
   bf(3) = ONE/3.0d0
   bf(4) = ONE/6.0d0
   cf(2) = 0.5d0
   cf(3) = 0.5d0
   cf(4) = ONE
   BTf = FARKodeButcherTable_Create(4, 4, 0, cf, Af, bf, df)
   
end subroutine create_rk4_butcher


subroutine set_integrators_wrapper(sun_method, t0, hs, retol, &
   abstol_c, abstol_ph, mxsteps, svec_dist, ctx, &
   arkode_mem, inner_arkode_mem, inner_stepper)

   use, intrinsic :: iso_c_binding
   implicit none
   
   character(len=80), intent(in) :: sun_method
   real(dp), intent(in) :: hs
   real(c_double), intent(in) :: t0
   real(dp), intent(in) :: retol, abstol_c, abstol_ph
   integer(c_long), intent(in) :: mxsteps
   type(N_Vector), pointer, intent(inout) :: svec_dist
   type(c_ptr), intent(inout) :: ctx
   type(c_ptr), intent(inout) :: arkode_mem
   type(c_ptr), intent(inout) :: inner_arkode_mem, inner_stepper
   !local variables
   real(c_double) :: hs_c
   real(c_double) :: retol_c, abstol_c_c, abstol_ph_c
   real(c_double) :: septol_c
   hs_c = real(hs, c_double)
   retol_c = real(retol, c_double)
   abstol_c_c = real(abstol_c, c_double)
   abstol_ph_c = real(abstol_ph, c_double)
   septol_c = 10.0d0
   if (trim(sun_method) .eq. 'ark') then
      call set_ark_integrators(t0, retol_c, abstol_c_c, abstol_ph_c, mxsteps, &
         svec_dist, ctx, arkode_mem)
   elseif (trim(sun_method) .eq. 'erk') then
      call set_erk_integrators(t0, retol_c, abstol_c_c, abstol_ph_c, mxsteps, &
         svec_dist, ctx, arkode_mem)
   elseif (trim(sun_method) .eq. 'mri-hs') then
      call set_mri_integrators(t0, retol_c, abstol_c_c, abstol_ph_c, mxsteps, &
         svec_dist, ctx, arkode_mem, inner_arkode_mem, inner_stepper, hs_c)
   elseif (trim(sun_method) .eq. 'mri') then
      call set_mri_integrators(t0, retol_c, abstol_c_c, abstol_ph_c, mxsteps, &
         svec_dist, ctx, arkode_mem, inner_arkode_mem, inner_stepper)
   elseif (trim(sun_method) .eq. 'mri-septol') then
      call set_mri_integrators(t0, retol_c, abstol_c_c, abstol_ph_c, mxsteps, &
         svec_dist, ctx, arkode_mem, inner_arkode_mem, inner_stepper, septol=septol_c)
   else
      call errore('sundials', 'need to select ark, erk, mri, mri-hs or mri-septol for method_flag', 1)
   endif
   

end subroutine

subroutine set_ark_integrators(t0, retol, abstol_c, abstol_ph, mxsteps, &
         svec_dist, ctx, arkode_mem)
   !!! After initialization of svec_dist, create ark steps
   use, intrinsic :: iso_c_binding
   implicit none
   real(c_double), intent(in) :: t0
   real(c_double), intent(in) :: retol, abstol_c, abstol_ph
   integer(c_long), intent(in) :: mxsteps
   type(N_Vector), pointer, intent(inout) :: svec_dist
   type(c_ptr), intent(inout) :: ctx
   type(c_ptr), intent(inout) :: arkode_mem
   ! local variables
   type(N_Vector), pointer :: abvec
   integer(c_int) :: retval
   real(c_double) :: adapt_params(3)

   ! create fast integrater
   arkode_mem = FARKStepCreate(c_funloc(f), c_null_funptr, &
      t0, svec_dist, ctx)
   if (.not. c_associated(arkode_mem)) then
      call errore('sundials', 'Error: arkode_mem = NULL', 1)
   endif

   !retval = FARKStepSetAdaptivityMethod(arkode_mem, 4, 1, 0, adapt_params) 
   !if (retval /= 0) then
   !   call errore('sundials', 'Error in FERKStepSetAdaptivityMethod', 1)
   !endif
   
   ! associate user data 
   !retval = FARKStepSetUserData(arkode_mem, c_null_ptr)
   !if (retval /= 0) then
   !   call errore('sundials', 'Error in FARKStepSetUserData', 1)
   !endif
   ! set retol and abstol

   call set_tol_vec(abstol_c, abstol_ph, ctx, abvec)
   retval = FARKStepSVtolerances(arkode_mem, retol, abvec)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set tolerances', 1)
   endif

   retval = FARKStepSetFixedStepBounds(arkode_mem, 1.0d0, 1.0d0)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set ark fixed step bound', 1)
   endif

   retval = FARKStepSetAdaptivityMethod(arkode_mem, 0, 1, 1, adapt_params)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set ark adaptivity mehtod', 1)
   endif
   retval = FARKStepSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set ark interpolant type', 1)
   endif
   retval = FARKStepSetMaxNumSteps(arkode_mem, mxsteps)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set maximum number of steps', 1)
   endif
   
end subroutine
 

subroutine set_erk_integrators(t0, retol, abstol_c, abstol_ph, mxsteps, &
         svec_dist, ctx, arkode_mem)
   !!! After initialization of svec_dist, create erk steps
   use, intrinsic :: iso_c_binding
   implicit none
   real(c_double), intent(in) :: t0
   real(c_double), intent(in) :: retol, abstol_c, abstol_ph
   integer(c_long), intent(in) :: mxsteps
   type(N_Vector), pointer, intent(inout) :: svec_dist
   type(c_ptr), intent(inout) :: ctx
   type(c_ptr), intent(inout) :: arkode_mem
   ! local variables 
   integer(c_int) :: retval
   type(N_Vector), pointer :: abvec
   ! create rk4 butcher table
   type(c_ptr) :: BTf 
   real(c_double) :: adapt_params(3)
   ! create Erk step
   arkode_mem = FERKStepCreate(c_funloc(f), t0, svec_dist, ctx)
   if (.not. c_associated(arkode_mem)) then
      call errore('sundials', 'Error: arkode_mem = NULL', 1)
   end if
   ! set user
   ! associate user data 
   !endif
   
   call set_tol_vec(abstol_c, abstol_ph, ctx, abvec)
   retval = FERKStepSVtolerances(arkode_mem, retol, abvec)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set tolerances', 1)
   endif
   
end subroutine

subroutine set_mri_integrators(t0, retol, abstol_c, abstol_ph, mxsteps, &
         svec_dist, ctx, arkode_mem, inner_arkode_mem, inner_stepper, hs, septol)
   !!! After initialization of svec_dist, create fast integrater 
   !!! and set uptions
   use, intrinsic :: iso_c_binding
   implicit none
   real(c_double), intent(in) :: t0
   real(c_double), intent(in) :: retol, abstol_c, abstol_ph
   integer(c_long), intent(in) :: mxsteps
   type(N_Vector), pointer, intent(inout) :: svec_dist
   type(c_ptr), intent(inout) :: ctx
   type(c_ptr), intent(inout) :: inner_arkode_mem, arkode_mem
   type(c_ptr), intent(inout) :: inner_stepper
   real(c_double), intent(in), optional :: hs
   real(c_double), intent(in), optional :: septol
   ! local variables
   type(N_Vector), pointer :: abvec, abvec_outer
   integer(c_int) :: retval
   real(c_double) :: hf
   real(c_double), allocatable :: Af(:,:), bf(:), cf(:), df(:) ! Arrays for fast Butcher table, NOTE: must be in row-major order
   type(c_ptr) :: BTf, BTs
   real(c_double) :: adapt_params(3)
   type(c_ptr) :: mri_coup

   ! create fast integrater
   ! CHANGE BACK
   !inner_arkode_mem = FARKStepCreate(c_funloc(f), c_null_funptr, &
   
   inner_arkode_mem = FARKStepCreate(c_funloc(ff), c_null_funptr, &
      t0, svec_dist, ctx)
   if (.not. c_associated(inner_arkode_mem)) then
      call errore('sundials', 'Error: inner_arkode_mem = NULL', 1)
   end if

   retval = FARKStepSetAdaptivityMethod(inner_arkode_mem, 0, 1, 1, adapt_params)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set ark adaptivity mehtod', 1)
   endif
   retval = FARKStepSetMaxNumSteps(inner_arkode_mem, mxsteps)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set maximum number of steps', 1)
   endif
   retval = FARKStepSetInterpolantType(inner_arkode_mem, ARK_INTERP_LAGRANGE)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set ark interpolant type', 1)
   endif

   retval = FARKStepSetFixedStepBounds(inner_arkode_mem, 1.0d0, 1.0d0)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set ark fixed step bound', 1)
   endif

   call set_tol_vec(abstol_c, abstol_ph, ctx, abvec)

   retval = FARKStepSVtolerances(inner_arkode_mem, retol, abvec)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set fast step tolerances', 1)
   endif
   retval = FARKStepCreateMRIStepInnerStepper(inner_arkode_mem, inner_stepper)
   if (retval /= 0) then
      if (ionode) print *, retval
      call errore('sundials', 'Error: create MRI Step Inner Stepper', 1)
   endif
   ! create slow inegrator and set options
   arkode_mem = FMRIStepCreate(c_funloc(fs), c_null_funptr, t0, svec_dist, inner_stepper, ctx)
   if (.not. c_associated(arkode_mem)) then
      call errore('sudnials', 'Error: arkode_mem = NULL', 1)
   endif
   ! pass user_data to inner stepper

   ! set slow step size
   if (present(hs)) then
      retval = FMRIStepSetFixedStep(arkode_mem, hs)
      if (retval /= 0) then 
         call errore('sundials', 'Error in set slow step size', 1)
      endif
      retval = FARKStepSetMaxNumSteps(arkode_mem, mxsteps)
      if ( retval /= 0) then 
         call errore('sundials', 'Error in set maximum number of steps', 1)
      endif
   else
      if (present(septol)) then
         call set_tol_vec(septol*abstol_c, septol*abstol_ph, ctx, abvec_outer)
         retval = FARKodeSVtolerances(arkode_mem, septol*retol, abvec_outer)
      else
         retval = FARKodeSVtolerances(arkode_mem, retol, abvec)
      endif

      if ( retval /= 0) then 
         call errore('sundials', 'Error in set slow step tolerances', 1)
      endif
      retval = FARKStepSetFixedStepBounds(arkode_mem, 1.0d0, 1.0d0)
      if ( retval /= 0) then 
         call errore('sundials', 'Error in set ark fixed step bound', 1)
      endif
      retval = FARKStepSetAdaptivityMethod(arkode_mem, 0, 1, 1, adapt_params)
      if ( retval /= 0) then 
         call errore('sundials', 'Error in set ark adaptivity mehtod', 1)
      endif
   endif
   retval = FARKStepSetMaxNumSteps(arkode_mem, mxsteps)
   if ( retval /= 0) then 
      call errore('sundials', 'Error in set maximum number of steps', 1)
   endif

end subroutine

subroutine sundials_reset(sun_method, arkode_mem, tout, svec_dist)
   use, intrinsic :: iso_c_binding
   implicit none

   character(len=80), intent(in) :: sun_method
   type(c_ptr), intent(inout) :: arkode_mem
   real(c_double), intent(in) :: tout
   type(N_Vector), pointer, intent(inout) :: svec_dist
   !real(c_double), pointer :: sundist_reset(:)

   integer(c_int) :: retval

   !sundist_reset => FN_VGetArrayPointer(svec_dist_reset)

   if (trim(sun_method) .eq. 'ark') then
      retval = FARKStepReset(arkode_mem, tout, svec_dist)
      if ( retval /= 0 ) then
         call errore('sundials', "ARKStep Reset Error", 1)
      endif
   elseif (trim(sun_method) .eq. 'erk') then
      retval = FERKStepReset(arkode_mem, tout, svec_dist)
      if (retval /= 0) then
         write(*,*) 'retval', retval
         call errore('sundials', "ERKStep Reset Error", 1)
      endif
   elseif (trim(sun_method) .eq. 'mri') then
      retval = FMRIStepReset(arkode_mem, tout, svec_dist)
      if (retval /= 0) then
         call errore('sundials', "MRIStep Reset Error", 1)
      endif
   elseif (trim(sun_method) .eq. 'mri-hs') then
      retval = FMRIStepReset(arkode_mem, tout, svec_dist)
      if (retval /= 0) then
         call errore('sundials', "MRIStep Reset Error", 1)
      endif
   elseif (trim(sun_method) .eq. 'mri-septol') then
      retval = FMRIStepReset(arkode_mem, tout, svec_dist)
      if (retval /= 0) then
         call errore('sundials', "MRIStep Reset Error", 1)
      endif
  
   else
      call errore('sundials', "please select mri or ark for method_flag", 1)
   endif 

end subroutine sundials_reset
 
subroutine sundials_evolve(arkode_mem, tout, dtout, tf, svec_dist, tcur, sun_method)
   use, intrinsic :: iso_c_binding

   implicit none

   type(c_ptr), intent(inout) :: arkode_mem
   real(c_double) :: tcur(1)
   real(c_double), intent(inout) :: tout
   real(c_double), intent(in) :: dtout, tf
   type(N_Vector), pointer, intent(inout) :: svec_dist
   integer(c_int) :: retval
   character(len=80), intent(in) :: sun_method

   if (trim(sun_method) .eq. 'ark') then
      retval = FARKStepEvolve(arkode_mem, tout, svec_dist, tcur, ARK_NORMAL)
      if ( retval /= 0 ) then
         call errore('sundials', "ARKStep Evolve Error", 1)
      endif
   elseif (trim(sun_method) .eq. 'erk') then
      retval = FERKStepEvolve(arkode_mem, tout, svec_dist, tcur, ARK_NORMAL)
      if (retval /= 0) then
         write(*,*) 'retval', retval
         call errore('sundials', "ERKStep Evolve Error", 1)
      endif
   elseif (trim(sun_method) .eq. 'mri') then
      retval = FMRIStepEvolve(arkode_mem, tout, svec_dist, tcur, ARK_NORMAL)
      if (retval /= 0) then
         call errore('sundials', "MRIStep Evolve Error", 1)
      endif
   elseif (trim(sun_method) .eq. 'mri-hs') then
      retval = FMRIStepEvolve(arkode_mem, tout, svec_dist, tcur, ARK_NORMAL)
      if (retval /= 0) then
         call errore('sundials', "MRIStep Evolve Error", 1)
      endif
   elseif (trim(sun_method) .eq. 'mri-septol') then
      retval = FMRIStepEvolve(arkode_mem, tout, svec_dist, tcur, ARK_NORMAL)
      if (retval /= 0) then
         call errore('sundials', "MRIStep Evolve Error", 1)
      endif
  
   else
      call errore('sundials', "please select mri or ark for method_flag", 1)
   endif 
   tout = min(tout + dtout, tf)
end subroutine

subroutine sundials_stats_wrapper(sun_method, diagfile, arkode_mem, inner_arkode_mem)

   use, intrinsic :: iso_c_binding
   implicit none
   character(len=80), intent(in) :: sun_method
   type(c_ptr), intent(in) :: diagfile
   type(c_ptr), intent(in) :: arkode_mem
   type(c_ptr), intent(in) :: inner_arkode_mem
   integer(c_int)  :: ierr          ! error flag

   if (trim(sun_method) .eq. 'ark') then
      call ARKStepStats(arkode_mem)
      !ierr = FARKStepSetDiagnostics(arkode_mem, diagfile)
   elseif (trim(sun_method) .eq. 'erk') then
      call ERKStepStats(arkode_mem)
      
      !ierr = FERKStepSetDiagnostics(arkode_mem, diagfile)
   elseif (trim(sun_method) .eq. 'mri') then
      call MRIStepStats(arkode_mem, inner_arkode_mem)
   elseif (trim(sun_method) .eq. 'mri-hs') then
      call MRIStepStats(arkode_mem, inner_arkode_mem)
   elseif (trim(sun_method) .eq. 'mri-septol') then
      call MRIStepStats(arkode_mem, inner_arkode_mem)
      !ierr = FMRIStepSetDiagnostics(arkode_mem, diagfile)
      !ierr = FARKStepPrintAllStats(inner_arkode_mem, diagfile, SUN_OUTPUTFORMAT_TABLE)
   endif
end subroutine

subroutine MRIStepStats(arkode_mem, inner_arkode_mem)
 
   !======= Inclusions ===========
   use, intrinsic :: iso_c_binding
 
   !======= Declarations =========
   implicit none
    
   type(c_ptr), intent(in) :: arkode_mem ! solver memory structure
   type(c_ptr), intent(in) :: inner_arkode_mem ! fast solver memory structure
 
   integer(c_int)  :: ierr          ! error flag
 
   integer(c_long) :: nssteps(1)     ! num steps
   integer(c_long) :: nfsteps(1)      ! num steps attempted
   integer(c_long) :: nfe(1), nfe_inner(1)        ! num explicit function evals
   integer(c_long) :: nfi(1), nfi_inner(1)        ! num implicit function evals
   integer(c_long) :: nfevals(1), nfevals_inner(1)    ! num function evals
 
   real(c_double)  :: hlast(1)      ! last step size
   real(c_double)  :: tcur(1)       ! internal time reached
 
   ierr = FARKStepGetNumSteps(inner_arkode_mem, nfsteps)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumSteps', 1)
   endif

   ierr = FMRIStepGetNumRhsEvals(inner_arkode_mem, nfe_inner, nfi_inner)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumRhsEvals', 1)
   endif
   
   nfevals_inner=nfe_inner+nfi_inner

   ierr = FMRIStepGetNumSteps(arkode_mem, nssteps) !, nfsteps)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FMRIStepGetNumSteps', 1)
   endif

   ierr = FMRIStepGetNumRhsEvals(arkode_mem, nfe, nfi)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumRhsEvals', 1)
   endif
   
   nfevals=nfe+nfi
   ierr = FMRIStepGetLastStep(arkode_mem, hlast)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetLastStep', 1)
   endif
   
   ierr = FARKStepGetCurrentTime(arkode_mem, tcur)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetCurrentTime', 1)
   endif
   write(stdout, '(5x, a)') ' '
   write(stdout, '(5x, a)') ' General Solver Stats:'
   write(stdout, '(5x, a, i9)') 'Total slow steps taken    =', nssteps
   write(stdout, '(5x, a, i9)') 'Total internal steps taken (fast)  =', nfsteps
   write(stdout, '(5x, a, i9)') 'Total rhs function calls      =', nfevals
   write(stdout, '(5x, a, i9)') 'Total internal rhs function calls      =', nfevals_inner
   write(stdout, '(5x, a, es12.5)') 'Last internal step size       =', hlast
   write(stdout, '(5x, a, es12.5)') 'Current internal time         =', tcur
 
   return
 
end subroutine


subroutine ARKStepStats(arkode_mem)
 
   !======= Inclusions ===========
   use, intrinsic :: iso_c_binding
 
   !======= Declarations =========
   implicit none
    
   type(c_ptr), intent(in) :: arkode_mem ! solver memory structure
 
   integer(c_int)  :: ierr          ! error flag
 
   integer(c_long) :: nsteps(1)     ! num steps
   integer(c_long) :: nst_a(1)      ! num steps attempted
   integer(c_long) :: nfe(1)        ! num explicit function evals
   integer(c_long) :: nfi(1)        ! num implicit function evals
   integer(c_long) :: nfevals(1)    ! num function evals
   integer(c_long) :: netfails(1)   ! num error test fails
 
   real(c_double)  :: hinused(1)    ! initial step size
   real(c_double)  :: hlast(1)      ! last step size
   real(c_double)  :: hcur(1)       ! step size for next step
   real(c_double)  :: tcur(1)       ! internal time reached
 
 
   ierr = FARKStepGetNumSteps(arkode_mem, nsteps)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumSteps', 1)
   endif

   ierr = FARKStepGetNumStepAttempts(arkode_mem, nst_a)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumStepAttempts', 1)
   endif
   ierr = FARKStepGetNumRhsEvals(arkode_mem, nfe, nfi)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumRhsEvals', 1)
   endif
   nfevals=nfe+nfi
   ierr = FARKStepGetActualInitStep(arkode_mem, hinused)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetActualInitStep', 1)
   endif
   ierr = FARKStepGetLastStep(arkode_mem, hlast)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetLastStep', 1)
   endif
   ierr = FARKStepGetCurrentStep(arkode_mem, hcur)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetCurrentStep', 1)
   endif
   ierr = FARKStepGetCurrentTime(arkode_mem, tcur)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetCurrentTime', 1)
   endif
   ierr = FARKStepGetNumErrTestFails(arkode_mem, netfails)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumErrTestFails', 1)
   endif
   
   write(stdout, '(5x, a)') ' '
   write(stdout, '(5x, a)') ' General Solver Stats:'
   write(stdout, '(5x, a, i9)') 'Total steps taken    =', nsteps
   write(stdout, '(5x, a, i9)') 'Total steps attempts =', nst_a
   write(stdout, '(5x, a, i9)') 'Total rhs function calls =', nfevals
   write(stdout, '(5x, a, i9)') 'Num error test failures  =', netfails
   write(stdout, '(5x, a, es12.5)') 'First internal step size  =', hinused
   write(stdout, '(5x, a, es12.5)') 'Next internal step size   =', hcur
   write(stdout, '(5x, a, es12.5)') 'Current internal time =', tcur
 
   return
 
end subroutine

subroutine ERKStepStats(arkode_mem)
 
   !======= Inclusions ===========
   use, intrinsic :: iso_c_binding
 
   !======= Declarations =========
   implicit none
    
   type(c_ptr), intent(in) :: arkode_mem ! solver memory structure
 
   integer(c_int)  :: ierr          ! error flag
 
   integer(c_long) :: nsteps(1)     ! num steps
   integer(c_long) :: nst_a(1)      ! num steps attempted
   integer(c_long) :: nfevals(1)    ! num function evals
   integer(c_long) :: netfails(1)   ! num error test fails
   
   real(c_double)  :: hinused(1)    ! initial step size
   real(c_double)  :: hlast(1)      ! last step size
   real(c_double)  :: hcur(1)       ! step size for next step
   real(c_double)  :: tcur(1)       ! internal time reached
   !type(N_Vector)  :: ele           ! local truncation error
   !real(c_double), pointer :: ele_ptr(:)
   real(c_double)  :: tolsfac(1)    ! user tolerance scaling

   !ele_ptr => FN_VGetArrayPointer(ele)

   ierr = FERKStepGetTolScaleFactor(arkode_mem, tolsfac)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetTolScaleFactor', 1)
   endif
   !ierr = FERKStepGetEstLocalErrors(arkode_mem, ele)
   !if (ierr /= 0) then
   !   call errore('sundials', 'Error in FARKStepGetEstLocalErrors', 1)
   !endif

   ierr = FERKStepGetNumSteps(arkode_mem, nsteps)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumSteps', 1)
   endif

   ierr = FERKStepGetNumStepAttempts(arkode_mem, nst_a)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumStepAttempts', 1)
   endif
   
   ierr = FERKStepGetNumRhsEvals(arkode_mem, nfevals)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumRhsEvals', 1)
   endif

   ierr = FERKStepGetActualInitStep(arkode_mem, hinused)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetActualInitStep', 1)
   endif
   
   ierr = FERKStepGetLastStep(arkode_mem, hlast)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetLastStep', 1)
   endif

   ierr = FERKStepGetCurrentStep(arkode_mem, hcur)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetCurrentStep', 1)
   endif

   ierr = FERKStepGetCurrentTime(arkode_mem, tcur)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetCurrentTime', 1)
   endif

   ierr = FERKStepGetNumErrTestFails(arkode_mem, netfails)
   if (ierr /= 0) then
      call errore('sundials', 'Error in FARKStepGetNumErrTestFails', 1)
   endif

   write(stdout, '(5x, a)') ' '
   write(stdout, '(5x, a)') ' General Solver Stats:'
   write(stdout, '(5x, a, i9)') 'Total steps taken    =', nsteps
   write(stdout, '(5x, a, i9)') 'Total steps attempts =', nst_a
   write(stdout, '(5x, a, i9)') 'Total rhs function calls =', nfevals
   write(stdout, '(5x, a, i9)') 'Num error test failures  =', netfails
   write(stdout, '(5x, a, es12.5)') 'Suggest Tolerance Scale Factor=', tolsfac
   write(stdout, '(5x, a, es12.5)') 'First internal step size  =', hinused
   write(stdout, '(5x, a, es12.5)') 'Last internal step size   =', hlast
   write(stdout, '(5x, a, es12.5)') 'Next internal step size   =', hcur
   write(stdout, '(5x, a, es12.5)') 'Current internal time =', tcur
   !write(stdout, '(5x, a, es12.5)'), 'Local error |ynp-yn|          =',norm2(ele_ptr)
 
   return
 
end subroutine

subroutine free_sundials_memory(sun_method, arkode_mem, ctx, &
         inner_arkode_mem, inner_stepper)
   use, intrinsic :: iso_c_binding
   implicit none
   character(len=80), intent(in) :: sun_method
   type(c_ptr), intent(inout) :: arkode_mem
   type(c_ptr), intent(inout) :: ctx
   type(c_ptr), intent(inout), optional :: inner_arkode_mem, inner_stepper
   ! local variables
   integer(c_int) :: ierr   ! error check for sundials

   if ((trim(sun_method) .eq. 'mri') .or. (trim(sun_method) .eq. 'mri-hs') &
         .or. (trim(sun_method) .eq. 'mri-septol' )) then
      call FARKStepFree(inner_arkode_mem)
      ierr = FMRIStepInnerStepper_Free(inner_stepper)
      if (ierr /= 0) then
         call errore('sundials', 'fail to free MRI inner stepper', 1)
      endif
      call FMRIStepFree(arkode_mem)
   endif

   if (trim(sun_method) .eq. 'ark' ) then
      call FARKStepFree(arkode_mem)
   endif
   if (trim(sun_method) .eq. 'erk' ) then
      call FERKStepFree(arkode_mem)
   endif

   ierr=FSUNContext_Free(ctx)
   if (ierr /= 0) then
      call errore('sundials', 'fail to free sundials context', 1)
   endif
end subroutine

subroutine deallocate_sunarrays()

   if(allocated(epcol)) deallocate(epcol)
   if(allocated(pecol)) deallocate(pecol)
   if(allocated(ph3col)) deallocate(ph3col)
   if(allocated(pecol_sum)) deallocate(pecol_sum)

end subroutine

#endif

end module sundials
