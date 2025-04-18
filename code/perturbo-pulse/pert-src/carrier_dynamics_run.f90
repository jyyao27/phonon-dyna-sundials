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

!> Ultrafast dynamics principal subroutine (calc_mode = 'dynamics-run')
!! Key steps: initialize the distribution, probagate for boltz_nstep time steps while
!! writing down the occupations into the HDF5 file.
subroutine carrier_dynamics_run()
   use, intrinsic :: iso_c_binding
#if defined (__SUNDIALS)
   use sundials, only : vec_pert_to_sun, vec_sun_to_pert, init_sundials, &
      fill_user_data, set_integrators_wrapper, sundials_stats_wrapper, & 
      sundials_evolve, deallocate_sunarrays, free_sundials_memory, sundials_reset
   use fsundials_core_mod

   use farkode_mod
   use farkode_mristep_mod
   use farkode_arkstep_mod
   use farkode_erkstep_mod 
   use fnvector_serial_mod
#endif
   use pert_const, only: dp, efield_tolerance, drift_accel_tolerance, bohr2ang, &
      timeunit
   use qe_mpi_mod, only: ionode, stdout
   use yaml_utils, only: ymlout
   use boltz_grid, only: grid, init_boltz_grid, init_boltz_qgrid, init_boltz_qgrid_vel_sym
   use boltz_grid_neighbors, only: set_grid_neighbors
   use pert_data,  only: epr_fid, qc_dim, kc_dim
   use pert_param, only: ntemper, temper, band_min, band_max, prefix, boltz_qdim, &
      boltz_nstep, time_step, output_nstep, boltz_efield, solver, &
      boltz_norm_dist, hole, & 
      boltz_acc_thr, boltz_nstep_min, boltz_init_dist, &
      scat_impl, pump_pulse, pump_pulse_fname, &
      nphstep, sun_method, dyna_phph, &
      adapt_smear_eph, adapt_smear_phph, &
      hs, sundials_mxsteps, retol, abstol_c, abstol_ph, &
      boltz_norm_energy, boltz_de, boltz_de_ph, divide_grid

   !
   use qe_mpi_mod, only: ionode, mp_barrier, mp_bcast, ionode_id, inter_pool_comm
   use boltz_dynamics_mod, only: output_dist, calc_carrier_population, &
      calc_mean_velocity_integral, get_energy_array, output_stop_dyna_message, &
      output_vel_accel, print_header_vel_accel, &
      calc_carrier_energy, &
      read_pump_pulse_data, pump_pulse_data, read_pump_pulse_snapshot, dump_pump_pulse_data
   use pert_output,only: progressbar_init, progressbar
   use boltz_dynamics_solver, only: runge_kutta_4th, euler, &
      runge_kutta_4th_eph, runge_kutta_4th_phph, &
      euler_eph, euler_phph
   !
   use boltz_scatter, only: boltz_scatter_setup, boltz_ph3scat_setup, boltz_ph3scat_setup_fast
   use boltz_scatter_integral_tgt, only: cdyna_scatter_target_setup
   use boltz_scatter_sizest
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, lattice_ifct, init_lattice_ifc, &
        init_lattice_ifct
   use hdf5_utils

   implicit none
   !local variables
   integer  :: i, nstep, istep, nestep, nestep_ph
   integer :: vel_accel_unit
   character(len=120) :: fname !< Name of the HDF5 files
   integer(HID_T) :: file_id   !< HDF5 file ID
   integer(HID_T) ::  group_id1!< HDF5 group ID for carriers
   integer(HID_T) ::  group_id2!< HDF5 group ID phonons
   logical :: calc_adv         !< calculate the advection part of BTE (if Efield is present)
   real(dp), allocatable:: dist0(:,:)   !< Carrier occupation function of the previous step
   real(dp), allocatable:: dist1(:,:)   !< Carrier occupation function of the current step
   real(dp), allocatable:: dist_t1(:,:) !< Temporary occupation
   real(dp), allocatable:: dist_t2(:,:) !< Temporary occupation
   real(dp), allocatable:: phdist0(:,:) !< Phonon occupation function of the previous step
   real(dp), allocatable:: phdist1(:,:) !< Phonon occupation function of the current step
   real(dp), allocatable:: phdist_t1(:,:) !< Temporary occupation
   real(dp), allocatable:: phdist_t2(:,:) !< Temporary occupation
   real(dp), allocatable:: phdist_pht2(:,:) !< ph-ph derivative

   real(dp), allocatable :: deriv(:,:) !< Full derivative of the occupation function including advection and collision terms
   real(dp), allocatable :: ene(:)     !< Arrays of energies on a grid
   real(dp), allocatable :: ene_ph(:)     !< Arrays of energies on a grid for phonons
   real(dp), allocatable :: popu(:)    !< Carrier population: occupation integrated over BZ as a function of energy
   real(dp), allocatable :: boltz_efield_tmp(:) !< Electric field
   real(dp), allocatable :: vel_time(:)         !< Array of the drift velocities over time
   real(dp), allocatable :: accel_time(:)      !< Array of the drift acceleration over time
   real(dp) :: cnum !< Carrier number
   real(dp) :: cnum0 !< Carrier number at time zero for normalization
   real(dp) :: enum0 !< Total energy at time zero for normalization
   real(dp) :: ecnum !< Total carrier energy
   real(dp) :: ecnum0 !< Total carrier energy at time zero for normalization
   real(dp) :: ephnum !< Total phonon energy
   real(dp) :: ephnum0 !< Total phonon energy at time zero for normalization
   real(dp) :: carrier_norm_factor !< Distribution normalization factor for phonons
   real(dp) :: phonon_norm_factor !< Distribution normalization factor for phonons
   real(dp) :: drift_vel !< Drift velocity at given time
   !
   type(grid) :: kg              !< k-grid
   type(lattice_ifc)   :: phon   !< harmonic phonon
   type(grid) :: qg              !< q-grid
   type(lattice_ifct)  :: phont  !< anharmonic phonon
   type(electron_wann) :: elec   !< Wannier
   ! pump pulse excitation
   type(pump_pulse_data) :: pump_data                  !< Structure to store pump pulse data
   real(dp) :: pump_cnum                               !< Additional carrier number at a given time step
   character(len=6), external :: int_to_char           !< Convert integer to character

#if defined (__SUNDIALS)
   !sundials
   type(N_Vector), pointer :: svec_dist ! solution vector
   type(c_ptr) :: ctx
   real(c_double), pointer :: sundist(:)
   real(c_double) :: t0
   integer(c_long) :: mxsteps
   type(c_ptr) :: inner_arkode_mem, arkode_mem
   type(c_ptr) :: inner_stepper

   real(c_double) :: tout, dtout, tf ! time output, time_step, final time
   real(c_double) :: tcur(1) ! current time
   real(c_double) :: hlast(1), hlast_inner(1) ! last steps
   integer(c_int) :: retval
   type(N_Vector), pointer :: svec_dist_reset
   real(c_double), pointer :: sundist_reset(:)
   integer(c_long) :: dist_length

   ctx = c_null_ptr
   inner_stepper = c_null_ptr
   inner_arkode_mem = c_null_ptr
   arkode_mem = c_null_ptr
#endif
   if (trim(solver) .eq. 'sundials') then
#if defined (__SUNDIALS)
      write(stdout,'(5x,a)') "Using SUNDIALS solvers"
#else
     call errore('carrier_dynamics_run', "SUNDIALS solvers are chosen but Perturbo not compiled with SUNDIALS", 1) 
#endif
   endif
   
   !only needs one temperature for phonon occupation N(mod, q)
   if(ntemper > 1) write(stdout,'(5x, a)') &
      "Warn (carrier_dynamics_run): only the first temperature is used."

   !init
   ! if normalizing both energy and distribution, do not exist boltz_norm_dist
   if ( boltz_norm_energy ) boltz_norm_dist = .false. 
   if ( (.not. dyna_phph ).and. boltz_norm_energy) &
      call errore('carrier_dynamics_run', "must use energy normalization with time-dependent phonon populations", 1)
   
   call init_lattice_ifc(epr_fid, qc_dim, phon)
   call init_electron_wann(epr_fid, kc_dim, elec)
   !setup k-grid
   call init_boltz_grid(kg, elec, band_min, band_max)
   
   !setup q-grid, and e-ph scattering channel (a.k.a k-q pair).
   call boltz_scatter_setup(kg, elec, phon, boltz_qdim)
   if (trim(scat_impl) .eq. 'tgt') then
      ! The target element oriented method needs to do additional setup
      call cdyna_scatter_target_setup(kg)
   endif
   !
   
   if( NORM2(boltz_efield(1:3)) > efield_tolerance) then
      calc_adv = .true.
      ! allocate space for df/dk
      allocate(deriv(kg%numb, kg%nk))
   else
      calc_adv = .false.
   end if

   ! if the external field is present, setup the neighbors to calculate df/dk
   ! if adaptive smearings are used, setup the neighbors to calculate group velocity
   if(calc_adv .or. adapt_smear_eph) then
      call set_grid_neighbors(kg)
   endif
   !allocate space for distribution function
   allocate( dist0(kg%numb, kg%nk), dist1(kg%numb, kg%nk), dist_t1(kg%numb, kg%nk) )

   ! read the pump pulse data from the HDF5 file
   if (pump_pulse) then
      call read_pump_pulse_data(pump_pulse_fname, pump_data)

      ! check if the pump pulse data is compatible with the current simulation
      if (pump_data%num_bands /= kg%numb .or. pump_data%num_kpoints /= kg%nk) then
         write(stdout,'(a, i9, a, i9)') 'Number of k-points in sim.: ', kg%nk, ' in pulse: ', pump_data%num_kpoints
         write(stdout,'(a, i5, a, i5)') 'Number of bands in sim.: ', kg%numb, ' in pulse: ', pump_data%num_bands
         call errore('carrier_dynamics_run', 'Pump pulse data is not compatible with the current simulation.', 1)
      endif

      ! check the time steps
      if (pump_data%num_steps > boltz_nstep) call errore('carrier_dynamics_run', 'Pump pulse has more time steps than the simulation.', 1)
      if (abs(pump_data%time_step-time_step) > 1.0E-6_dp .and. pump_data%num_steps > 1) &
         write(stdout,'(a)') 'WARNING: Pump pulse time step is different from simulation time step.'

   endif

   if (pump_pulse .and. boltz_norm_dist) then
      call errore('carrier_dynamics_run', 'Pump pulse (pump_pulse = .true.) ' // &
         'and normalization of the distribution (boltz_norm_dist = .true.) are not compatible.', 1)
   endif
  
   !!! swap from full_dyna_run all phdist_t1 to phdist_t2, and dist_t1 to dist_t2, vice versa
   ! allocate additional workspace if solver is rk4
   if(trim(solver) .eq. 'rk4') & 
      allocate( dist_t2(kg%numb, kg%nk) )

   if (ionode) then
      ! Output details about our dist0/dist1/dist_t1/dist_t2 allocations
      call output_dist_alloc_stats_to_yaml(dist0, merge(3, 4, trim(solver) .eq. 'rk4'))
      ! Once we have done setup we should have a good idea of how much memory is required
      call boltz_scatter_estimate_pools(64.0_dp, 39.0_dp)
   endif

   dist0 = 0.0E0_dp
   dist1 = 0.0E0_dp
   
   !initialize distribution function: 
   !   restart from previous run or start a new one.
   fname = trim(prefix)//"_cdyna.h5"
   
   call start_clock('init_dyna')
   if ( dyna_phph ) then
      call init_boltz_qgrid(qg, phon)
      write(stdout,'(5x,a)') "Running dynamics with anharmonic phonon interaction"
      if ( adapt_smear_phph ) then
         call set_grid_neighbors(qg)
         call init_boltz_qgrid_vel_sym(qg, phon, .True.) 
      endif
      call init_lattice_ifct(epr_fid, qc_dim, phont, .true.)
      !call init_lattice_ifct(epr_fid, qc_dim, phont)
      if (divide_grid) then
         call boltz_ph3scat_setup_fast(qg, phon, phont, divide_grid)
      else
         call boltz_ph3scat_setup(qg, phon, phont)
      endif
      !call boltz_ph3scat_setup(qg, phon, phont, boltz_qdim)
      allocate( phdist0(qg%numb,qg%nk), phdist1(qg%numb,qg%nk), phdist_t1(qg%numb, qg%nk) )
      phdist_t1 = 0.0_dp
      allocate( phdist_pht2(qg%numb,qg%nk)) ! ph-ph only derivative
      if(trim(solver) .eq. 'rk4') then
         allocate( phdist_t2(qg%numb, qg%nk) )
         phdist_t2 = 0.0E0_dp
      endif
      phdist_t1 = 0.0E0_dp
      phdist0 = 0.0E0_dp
      phdist1 = 0.0E0_dp
      phdist_pht2 = 0.0E0_dp ! ph-ph only derivative
      call cdyna_setup(fname, kg, file_id, group_id1, dist1, &
         group_id2, qg, phdist1)
   else
      ! get initial distribution
      call cdyna_setup(fname, kg, file_id, group_id1, dist1, group_id2)
   endif
   call stop_clock('init_dyna')

   if( ionode ) write(ymlout, '(/,a)') 'dynamics-run:'
   ! setup the energy grid
   if(boltz_norm_dist .or. (boltz_acc_thr > drift_accel_tolerance) .or. pump_pulse) then
      ene = get_energy_array(kg)
      nestep = size(ene)
   endif

   ! calculate the initial carrier concentration
   if(boltz_norm_dist .or. pump_pulse) then
      allocate(popu(nestep))
      call calc_carrier_population(kg, dist1, ene, popu, hole, cnum0)

   endif

   if (boltz_norm_energy) then
      ! sum up electron population
      ene = get_energy_array(kg)
      nestep = size(ene)
      allocate(popu(nestep))
      call calc_carrier_population(kg, dist1, ene, popu, hole, cnum0)
      call calc_carrier_energy(kg, dist1, ene, boltz_de, ecnum0, hole)

      ! sum up phonon population   
      nestep_ph = int( floor((maxval(qg%enk)-minval(qg%enk))/boltz_de_ph) ) + 3
      allocate(ene_ph(nestep_ph))
      do i = 1, nestep_ph
         ene_ph(i) = minval(qg%enk) + (i-2)*boltz_de_ph
      enddo
      call calc_carrier_energy(qg, phdist1, ene_ph, boltz_de_ph, ephnum0)
      enum0 = ecnum0 + ephnum0

   endif


   if (boltz_acc_thr > drift_accel_tolerance) then
      allocate(boltz_efield_tmp(3))

      ! if the efield is specified, take its value
      ! otherwise, set to (1.0, 1.0, 1.0)
      if( NORM2(boltz_efield(1:3)) < efield_tolerance) then
         boltz_efield_tmp = (/1.0_dp, 1.0_dp, 1.0_dp /)
      else
         boltz_efield_tmp(:) = boltz_efield(:)
      end if

      ! drift velocity and acceleration as a function of time
      allocate(vel_time(boltz_nstep))
      allocate(accel_time(boltz_nstep))

      vel_time = 0.0_dp
      accel_time = 0.0_dp

      if(ionode) call print_header_vel_accel(vel_accel_unit)

   endif

   !start dynamics simulation, write the initial step
   if(ionode) call progressbar_init('Carrier Dynamics:')
   nstep = boltz_nstep
   call start_clock('dynamics_tot')
   !
   ! sundials initialization - kelly yao
#if defined (__SUNDIALS)
   if (trim(solver) .eq. 'sundials') then
      if ( dyna_phph ) then
         call fill_user_data(kg, temper(1), qg) 
         call init_sundials(dist1, svec_dist, sundist, ctx, phdist1)
      else
         ! if not dyna_phph, check whether we are using MRI. MRI doesn't make sense if not
         if ( .not. dyna_phph) then
            !if (trim(sun_method) .eq. 'mri') then
            if ((trim(sun_method) .eq. 'mri') .or. (trim(sun_method) .eq. 'mri-hs') &
               .or. (trim(sun_method) .eq. 'mri-septol')) then
               call errore('carrier_dynamics_run', 'Does not support using MRI without phonon dynamics', 1)
            endif
         endif
         call fill_user_data(kg, temper(1)) 
         call init_sundials(dist1, svec_dist, sundist, ctx)
      endif
      ! set time
      t0 = 0.0d0
      tcur = t0
      tf = t0 + time_step*nstep
      dtout = time_step
      tout = t0 + dtout
      mxsteps = int(sundials_mxsteps, c_long)
      call set_integrators_wrapper(sun_method, t0, hs, retol, &
         abstol_c, abstol_ph, mxsteps, svec_dist, ctx, &
         arkode_mem, inner_arkode_mem, inner_stepper)
      ! set reset vector
      dist_length = kg%numb * kg%nk
      if ( dyna_phph ) then
         dist_length = dist_length + qg%numb * qg%nk
      endif
      !svec_dist_reset => FN_VNew_Serial(dist_length, ctx)
      !if ( .not.associated(svec_dist_reset)) then 
      !   call errore('sundials', 'Error creating the Nvector for data', 1)
      !endif
      !sundist_reset => FN_VGetArrayPointer(svec_dist_reset)
   endif
#endif
   !istep = 1
   do istep = 1, nstep
      dist0(:,:) = dist1(:,:)

      ! add the pump pulse to the distribution
      if (pump_pulse) then
         if (istep .le. pump_data%num_steps) then
            call read_pump_pulse_snapshot(pump_data, istep)
            dist0(:,:) = dist0(:,:) + pump_data%dist(:,:)

            call calc_carrier_population(kg, pump_data%dist, ene, popu, hole, pump_cnum)
            pump_data%cnum_time(istep) = pump_cnum

         endif
      endif
      if ( dyna_phph ) phdist0(:,:) = phdist1(:,:)

      if(ionode) write(*,*) 'time step :', istep
      !
      !get the f(t_i+1)
      if(trim(solver) .eq. 'rk4') then
         if ( dyna_phph ) then
         !       separate eph and phph, and keep phdist_pht2 constant unless mod(i-1, nphstep)
!       .eq. 0
         ! disable eph for now
            call runge_kutta_4th_eph &
               (kg, qg, time_step, dist0, dist1, dist_t2, dist_t1, &
               phdist0, phdist1, phdist_t2, phdist_t1)
          
            call runge_kutta_4th_phph &
               (qg, time_step, phdist0, phdist1, phdist_t2, phdist_pht2, &
               istep, nphstep)
         else
            
            call runge_kutta_4th &
               (kg, time_step, dist0, dist1, dist_t1, dist_t2, deriv, temper(1), calc_adv)
         endif

      elseif(trim(solver) .eq. 'euler') then
         if ( dyna_phph ) then
         !   call errore('full_dyna_setup','only rk4 method is supported at this moment', 5)
            call euler_eph(kg, qg, time_step, dist0, dist1, dist_t1, &
               phdist0, phdist1, phdist_t1) ! no dist_t2
          
            call euler_phph(qg, time_step, phdist0, phdist1, phdist_pht2, &
               istep, nphstep)
         else
            call euler(kg, time_step, dist0, dist1, dist_t1, deriv, temper(1), calc_adv)
         
         endif

#if defined (__SUNDIALS)
      elseif(trim(solver) .eq. 'sundials') then
         call sundials_evolve(arkode_mem, tout, time_step, tf, svec_dist, tcur, sun_method)
#endif
      else
         call errore('carrier_dynamics_run', 'only rk4 and forward euler &
           is supported at this moment',5)
      endif
      !
      ! normalize distribution
      if(boltz_norm_dist) then
         call start_clock('norm_dist')
#if defined (__SUNDIALS)
         if (trim(solver) .eq. 'sundials') then
            if ( dyna_phph ) then 
               call vec_sun_to_pert(dist1, sundist, phdist1)
            else
               call vec_sun_to_pert(dist1, sundist)
            endif
         endif
#endif
         call calc_carrier_population(kg, dist1, ene, popu, hole, cnum)

         dist1(:,:) = dist1(:,:) / cnum * cnum0
#if defined (__SUNDIALS)
         if (trim(solver) .eq. 'sundials') then
            ! fill reset vector
            sundist => FN_VGetArrayPointer(svec_dist)
            if ( dyna_phph) then
               call vec_pert_to_sun(dist1, sundist, phdist1)
            else
               call vec_pert_to_sun(dist1, sundist)
            endif 
            !if (ionode) write(*,*) sundist_reset(1)
            if (trim(sun_method) .eq. 'mri') then
               retval = FARKStepReset(inner_arkode_mem, tout-dtout, svec_dist)
               if (retval /= 0) then
                  call errore('sundials', "ARKStep Reset Error", 1)
               endif
               retval = FMRIStepReset(arkode_mem, tout-dtout, svec_dist)
               if (retval /= 0) then
                  call errore('sundials', "MRIStep Reset Error", 1)
               endif
            endif
   
            if (trim(sun_method) .eq. 'ark') then
               retval = FARKStepReset(arkode_mem, tout-dtout, svec_dist)
               if (retval /= 0) then
                  call errore('sundials', "ARKStep Reset Error", 1)
               endif
            endif
            !call sundials_reset(sun_method, arkode_mem, tout, svec_dist)!, sundist_reset)
         endif
#endif
         call stop_clock('norm_dist')
      endif

      if(boltz_norm_energy) then
         call start_clock('norm_energy')
#if defined (__SUNDIALS)
         if (trim(solver) .eq. 'sundials') then
            if ( dyna_phph ) then 
               call vec_sun_to_pert(dist1, sundist, phdist1)
            else
               call vec_sun_to_pert(dist1, sundist)
            endif
         endif
#endif
         call calc_carrier_population(kg, dist1, ene, popu, hole, cnum)
         dist1(:,:) = dist1(:,:) * cnum0 / cnum
         call calc_carrier_energy(kg, dist1, ene, boltz_de, ecnum, hole)
         call calc_carrier_energy(qg, phdist1, ene_ph, boltz_de_ph, ephnum)
         phdist1(:,:) = phdist1(:,:) * (enum0 - ecnum) / ephnum
         !call calc_carrier_energy(qg, phdist1, ene_ph, boltz_de_ph, ephnum)

#if defined (__SUNDIALS)
         if (trim(solver) .eq. 'sundials') then
            ! fill reset vector
            sundist => FN_VGetArrayPointer(svec_dist)
            if ( dyna_phph) then
               call vec_pert_to_sun(dist1, sundist, phdist1)
            else
               call vec_pert_to_sun(dist1, sundist)
            endif 
            !if (ionode) write(*,*) sundist_reset(1)
            !if (trim(sun_method) .eq. 'mri') then
            if ((trim(sun_method) .eq. 'mri') .or. (trim(sun_method) .eq. 'mri-hs' ) &
               .or. (trim(sun_method) .eq. 'mri-septol')) then
               retval = FARKStepReset(inner_arkode_mem, tout-dtout, svec_dist)
               if (retval /= 0) then
                  call errore('sundials', "ARKStep Reset Error", 1)
               endif
               retval = FMRIStepReset(arkode_mem, tout-dtout, svec_dist)
               if (retval /= 0) then
                  call errore('sundials', "MRIStep Reset Error", 1)
               endif
            endif
   
            if (trim(sun_method) .eq. 'ark') then
               retval = FARKStepReset(arkode_mem, tout-dtout, svec_dist)
               if (retval /= 0) then
                  call errore('sundials', "ARKStep Reset Error", 1)
               endif
            endif
            !call sundials_reset(sun_method, arkode_mem, tout, svec_dist)!, sundist_reset)
         endif
#endif

         call stop_clock('norm_energy')
      endif

      !output f(t_i+1) if needed.
      if( mod(istep, output_nstep) .eq. 0 ) then
         if(ionode) then
            call progressbar(istep, boltz_nstep)
#if defined (__SUNDIALS)
            if (trim(solver) .eq. 'sundials') then
               if ( dyna_phph ) then 
                  call vec_sun_to_pert(dist1, sundist, phdist1)
               else
                  call vec_sun_to_pert(dist1, sundist)
               endif
            endif
#endif
            call output_dist(kg, group_id1, dist1, (istep/output_nstep) )
            if ( dyna_phph ) then 
               call output_dist(qg, group_id2, phdist1, (istep/output_nstep) )
            endif
         endif
      endif

      ! if accel. threshold is specified, stop the simulation when accel.
      ! is smaller than boltz_acc_thr
      if (boltz_acc_thr > drift_accel_tolerance) then

         call start_clock('drift_accel')
         drift_vel = calc_mean_velocity_integral(kg, dist1, ene, boltz_efield_tmp, hole)
         vel_time(istep) = drift_vel
         
         if (istep > 1) then

            ! drift acceleration in cm/s2
            accel_time(istep) = (drift_vel - vel_time(istep - 1)) / &
                                     (time_step * timeunit * 1.0E15_dp)

            if(ionode) call output_vel_accel(vel_accel_unit, istep*time_step*(timeunit*1.0E15_dp), drift_vel, accel_time(istep))

         endif

         if (istep > 10) then

            if(all(abs(accel_time(istep-9:istep)) < boltz_acc_thr) &
               .and. istep >= boltz_nstep_min) then
               if(ionode) call output_stop_dyna_message(istep, drift_vel, accel_time(istep), vel_time, accel_time)
               exit
            endif

         endif

         call stop_clock('drift_accel')

      endif

      ! exit
      !if (istep >= nstep) exit
      !istep = istep + 1

   enddo

   ! close files
   if(ionode) then
#if defined (__SUNDIALS)
      if (trim(solver) .eq. 'sundials') then
         !call sundials_stats_wrapper(sun_method, diagfile, arkode_mem, inner_arkode_mem)
      endif
#endif
   endif
   !close file
   if(ionode) then
      call hdf_close_group(group_id1)
      if ( dyna_phph) then 
         call hdf_close_group(group_id2)
      endif
      call hdf_close_file(file_id)
      if (boltz_acc_thr > drift_accel_tolerance) close(vel_accel_unit)

      if (pump_pulse) then
         call dump_pump_pulse_data(pump_data)
         call hdf_close_group(pump_data%pump_dist_gid)
         call hdf_close_file(pump_data%pump_fid)
      endif
   endif

   
   ! release work space
   deallocate(dist0, dist1, dist_t1) 
   if ( dyna_phph ) deallocate(phdist0, phdist1, phdist_t1, phdist_pht2)
   if(trim(solver) .eq. 'rk4') then
      deallocate(dist_t2)
      if ( dyna_phph ) deallocate(phdist_t2)
   endif
   
   if(allocated(deriv)) deallocate(deriv)
   if(allocated(ene)) deallocate(ene)
   if(allocated(popu)) deallocate(popu)
   if(allocated(vel_time)) deallocate(vel_time)
   if(allocated(accel_time)) deallocate(accel_time)

#if defined (__SUNDIALS)
   if(trim(solver) .eq. 'sundials') then
      call FN_VDestroy_Serial(svec_dist)
      call free_sundials_memory(sun_method, arkode_mem, ctx, &
         inner_arkode_mem, inner_stepper)
      call deallocate_sunarrays()

      !call FSUNDIALSFileClose(diagfile)
   endif
#endif
   call stop_clock('dynamics_tot')

   return
1050 format(1x, f15.8, 2x, f15.8)

   contains
      !
   subroutine cdyna_setup(filename, kgrid, fid, gid, dist, &
                   qgid, qgrid, qdist)
      use pert_const, only: dp, timeunit, ryd2ev, ryd2mev
      use qe_mpi_mod, only: ionode, mp_barrier, mp_bcast, ionode_id, inter_pool_comm
      use boltz_dynamics_mod, only: output_dist, restart_dist, restart_dist_full, &
         init_dist_fermi, init_dist_lorentz, init_dist_gaussian, init_disp_bose, &
         init_dist_filled
      use pert_param, only: boltz_nstep, output_nstep, time_step, &
         boltz_init_dist, boltz_init_e0, boltz_init_smear, boltz_init_ampl, &
         dyna_phph, hole
      use boltz_grid, only: grid
      use hdf5_utils
      implicit none
      character(*), intent(in) :: filename
      type(grid), intent(in) :: kgrid
      type(grid), intent(in), optional :: qgrid
      integer(HID_T), intent(out) :: fid, gid, qgid
      real(dp), intent(out) :: dist( kgrid%numb, kgrid%nk )
      real(dp), intent(out), allocatable, optional :: qdist(:,:)
      !local variables
      logical :: has_file
      real(dp) :: t_step_fs, qtemp
      character(len=120) :: group_name, msg
      integer :: nrun, current_run
      character(len=6), external :: int_to_char
      real(dp) :: efield_norm2
      integer:: ib, ik
      real(dp) :: enk
   
      if ( present(qgrid) ) allocate(qdist( qgrid%numb, qgrid%nk) )
      gid = 0;  fid = 0; qgid = 0;  ! init
      qtemp = temper(1); !init phonon temperature
      !restart dynamics simulation from previous run
      !same indicator for full dyanmics
      if(trim(boltz_init_dist) .eq. 'restart') then
         inquire(file=trim(filename), exist=has_file)
         if(.not. has_file) call errore('cdyna_setup','missing '// trim(filename), 1)
         if(ionode) then
            call hdf_open_file(fid, trim(filename), status='OLD', action='READWRITE')
            call hdf_read_dataset(fid, 'num_runs', nrun)
            if ( dyna_phph ) then 
               call restart_dist_full(kgrid, qgrid, fid, nrun, dist, qdist)
            else
               call restart_dist(kgrid, fid, nrun, dist)
            endif
            current_run = nrun + 1
            !
            call hdf_update_dataset(fid, 'num_runs', current_run)
         endif
         call mp_bcast(dist, ionode_id, inter_pool_comm)
         if ( dyna_phph ) call mp_bcast(qdist, ionode_id, inter_pool_comm)
            !
      else
            ! for carrier initialization
         select case (boltz_init_dist)
         case('fermi')
            call init_dist_fermi(kgrid, dist, boltz_init_e0, boltz_init_smear)
         case ('lorentz')
            call init_dist_lorentz(kgrid, dist, boltz_init_e0, boltz_init_smear, hole)
         case('gaussian')
            call init_dist_gaussian(kgrid, dist, boltz_init_e0, boltz_init_smear, boltz_init_ampl, hole)
         case('filled')
            call init_dist_filled(kgrid, dist, boltz_init_e0)
         case default
            msg = "invalid boltz_init_dist, valid options are 'restart', 'fermi', 'lorentz', 'gaussian', 'filled'."
            call errore('carrier_dynamics_run', trim(msg), 1)
         end select
            
         if ( dyna_phph ) then 
            ! for phonon initialization, other functions could also be adopted
            call init_disp_bose(qgrid, qdist, qtemp) !qtemp may be changed to an separate input parameter later -xt
         endif

         if(ionode) then
            current_run = 1
            call hdf_open_file(fid, trim(filename), status='NEW')
            call hdf_write_dataset(fid, 'num_runs', current_run)
            !write kg%enk(:,:)
            call hdf_write_dataset(fid, 'band_structure_ryd', kgrid%enk)
            call hdf_write_attribute(fid, 'band_structure_ryd', 'ryd2ev', ryd2ev)

            if ( dyna_phph ) then 
               !write qg%enk(:,:)
               call hdf_write_dataset(fid, 'phonon_dispersion_ryd', qgrid%enk)
               call hdf_write_attribute(fid, 'phonon_dispersion_ryd', 'ryd2ev', ryd2ev)
            endif

         endif
      endif
      
      !record status of the current run
      if(ionode) then
         if (dyna_phph) then
            group_name = "carrier_dynamics_run_" // trim( int_to_char(current_run) )
         else
            group_name = "dynamics_run_" // trim( int_to_char(current_run) )
         endif
         call hdf_create_group(fid, trim(group_name))
         call hdf_open_group(fid, trim(group_name), gid)
         !
         t_step_fs = time_step * output_nstep * (timeunit*1.0E15_dp)
         call hdf_write_dataset(gid, 'time_step_fs', t_step_fs)
         call hdf_write_dataset(gid, 'num_steps',  0)
         
         ! nvfortran version 22.7 has a bug that prevents us from using the NORM2()
         ! intrinsic here - even though it works perfectly fine in the enclosing
         ! subroutine.  So, we implement it manually for the time being.  Hopefully
         ! at some point we can switch back to the original way of doing this.
         !
         ! if( NORM2(boltz_efield(1:3)) > efield_tolerance) call hdf_write_dataset(gid, 'efield',  &
         !
         efield_norm2 = SQRT(boltz_efield(1) * boltz_efield(1) + &
                             boltz_efield(2) * boltz_efield(2) + &
                             boltz_efield(3) * boltz_efield(3))

         !if(efield_norm2 > efield_tolerance) call hdf_write_dataset(gid, 'efield',  &
         if( NORM2(boltz_efield(1:3)) > efield_tolerance) call hdf_write_dataset(gid, 'efield',  &
                     boltz_efield / (bohr2ang*1E-8_dp) * ryd2ev )
         !
         if(current_run .eq. 1)  call output_dist(kgrid, gid, dist, 0)

     
         if ( dyna_phph ) then 
            group_name = "phonon_dynamics_run_" // trim( int_to_char(current_run) )
            call hdf_create_group(fid, trim(group_name))
            call hdf_open_group(fid, trim(group_name), qgid)
            !
            t_step_fs = time_step * output_nstep * (timeunit*1.0E15_dp)
            call hdf_write_dataset(qgid, 'time_step_fs', t_step_fs)
            call hdf_write_dataset(qgid, 'num_steps',  0)
            !
            if(current_run .eq. 1)  call output_dist(qgrid, qgid, qdist, 0)
         endif

      endif      
      !
      call mp_barrier(inter_pool_comm)
   end subroutine cdyna_setup
   !

   ! Output details regarding the memory allocation of the scatter_targets
   ! data structure.
   subroutine output_dist_alloc_stats_to_yaml(dist, count)
      implicit none
      real(dp), intent(in) :: dist(:,:)
      integer :: count

      integer(8) :: bytes_per_entry, total_bytes, subelems, total_subelems
 
     ! dist0/dist1/dist_t1/dist_t2 arrays

      bytes_per_entry = STORAGE_SIZE(dist(1,1)) / 8
      total_bytes = count * bytes_per_entry * SIZE(dist, kind=8)

      write(ymlout,'(3x, a)') 'dist arrays:'
      write(ymlout,'(6x, a, i5)') 'bytes per element: ', bytes_per_entry
      write(ymlout,'(6x, a, i12, a, i12)') 'elements: ', SIZE(dist, 1), ' x ', SIZE(dist, 2)
      write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
      write(ymlout,'(6x, a, i3/)') 'allocations: ', count

      call boltzscat_accum_cpugpu_bytes(total_bytes)
   end subroutine output_dist_alloc_stats_to_yaml

end subroutine carrier_dynamics_run
