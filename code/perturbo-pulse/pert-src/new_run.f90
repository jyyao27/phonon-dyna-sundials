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
   use pert_const, only: dp, efield_tolerance, drift_accel_tolerance, bohr2ang, &
      timeunit
   use qe_mpi_mod, only: ionode, stdout
   use yaml_utils, only: ymlout
   use boltz_grid, only: grid, init_boltz_grid
   use boltz_grid_neighbors, only: set_grid_neighbors
   use pert_data,  only: epr_fid, qc_dim, kc_dim
   use pert_param, only: ntemper, temper, band_min, band_max, prefix, boltz_qdim, &
      boltz_nstep, time_step, output_nstep, boltz_efield, solver, &
      boltz_norm_dist, hole, &
      boltz_acc_thr, boltz_nstep_min, boltz_init_dist, &
      scat_impl, pump_pulse, pump_pulse_fname
   !
   use boltz_dynamics_mod, only: output_dist, calc_carrier_population, &
      calc_mean_velocity_integral, get_energy_array, output_stop_dyna_message, &
      output_vel_accel, print_header_vel_accel, &
      read_pump_pulse_data, pump_pulse_data, read_pump_pulse_snapshot, dump_pump_pulse_data
   !
   use pert_output,only: progressbar_init, progressbar
   use boltz_dynamics_solver, only: runge_kutta_4th, euler
   !
   use boltz_scatter, only: boltz_scatter_setup
   use boltz_scatter_integral_tgt, only: scatter_target_setup
   use boltz_scatter_sizest
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   use hdf5_utils
   implicit none
   !local variables
   integer  :: i, nstep, istep, nestep
   integer :: vel_accel_unit
   character(len=120) :: fname !< Name of the HDF5 files
   integer(HID_T) :: file_id   !< HDF5 file ID
   integer(HID_T) :: group_id  !< HDF5 group ID
   logical :: calc_adv         !< calculate the advection part of BTE (if Efield is present)
   real(dp), allocatable:: dist0(:,:)   !< Carrier occupation function of the previous step
   real(dp), allocatable:: dist1(:,:)   !< Carrier occupation function of the current step
   real(dp), allocatable:: dist_t1(:,:) !< Temporary occupation
   real(dp), allocatable:: dist_t2(:,:) !< Temporary occupation
   real(dp), allocatable :: deriv(:,:)  !< Full derivative of the occupation function including advection and collision terms
   real(dp), allocatable :: ene(:)      !< Arrays of energies on a grid
   real(dp), allocatable :: popu(:)     !< Carrier population: occupation integrated over BZ as a function of energy
   real(dp), allocatable :: boltz_efield_tmp(:) !< Electric field
   real(dp), allocatable :: vel_time(:)         !< Array of the drift velocities over time
   real(dp), allocatable :: accel_time(:)       !< Array of the drift accelerations over time
   real(dp) :: cnum              !< Carrier number
   real(dp) :: cnum0             !< Carrier number at time zero for normalization
   real(dp) :: drift_vel         !< Drift velocity at given time
   ! objects
   type(grid) :: kg              !< k-grid object
   type(lattice_ifc)   :: phon   !< phonon
   type(electron_wann) :: elec   !< Wannier
   ! pump pulse excitation
   type(pump_pulse_data) :: pump_data                  !< Structure to store pump pulse data
   real(dp) :: pump_cnum                               !< Additional carrier number at a given time step
   !
   character(len=6), external :: int_to_char           !< Convert integer to character

   !only needs one temperature for phonon occupation N(mod, q)
   if(ntemper > 1) write(stdout,'(5x, a)') &
      "WARNING (carrier_dynamics_run): only the first temperature is used."
   !init
   call init_lattice_ifc(epr_fid, qc_dim, phon)
   call init_electron_wann(epr_fid, kc_dim, elec)
   !
   !setup k-grid
   call init_boltz_grid(kg, elec, band_min, band_max)
   !setup q-grid, and e-ph scattering channel (a.k.a k-q pair).
   call boltz_scatter_setup(kg, elec, phon, boltz_qdim)
   if (trim(scat_impl) .eq. 'tgt') then
      ! The target element oriented method needs to do additional setup
      call scatter_target_setup(kg,2)
   endif
   !
   if( NORM2(boltz_efield(1:3)) > efield_tolerance) then
      calc_adv = .true.
      ! allocate space for df/dk
      allocate(deriv(kg%numb, kg%nk))
   else
      calc_adv = .false.
   endif

   ! if the external field is present, setup the neighbors to calculate df/dk
   if(calc_adv) call set_grid_neighbors(kg)

   ! allocate space for distribution function
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

   ! allocate additional workspace if solver is rk4
   if(trim(solver) .eq. 'rk4') then
      allocate( dist_t2(kg%numb, kg%nk) )
   endif

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

   ! get initial distribution
   call cdyna_setup(fname, kg, file_id, group_id, dist1)

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
   nstep = (boltz_nstep / output_nstep) * output_nstep

   call start_clock('dynamics_tot')

   !
   do istep = 1, nstep
      dist0(:,:) = dist1(:,:)

      !
      ! add the pump pulse to the distribution
      if (pump_pulse) then
         if (istep .le. pump_data%num_steps) then
            call read_pump_pulse_snapshot(pump_data, istep)
            dist0(:,:) = dist0(:,:) + pump_data%dist(:,:)

            call calc_carrier_population(kg, pump_data%dist, ene, popu, hole, pump_cnum)
            pump_data%cnum_time(istep) = pump_cnum

         endif
      endif

      !
      ! get the f(t_i+1)
      if(trim(solver) .eq. 'rk4') then
         call runge_kutta_4th &
            (kg, time_step, dist0, dist1, dist_t1, dist_t2, deriv, temper(1), calc_adv)
      else
         call euler(kg, time_step, dist0, dist1, dist_t1, deriv, temper(1), calc_adv)
      endif

      !
      ! normalize distribution
      if(boltz_norm_dist) then
         call start_clock('norm_dist')
         call calc_carrier_population(kg, dist1, ene, popu, hole, cnum)

         dist1(:,:) = dist1(:,:) / cnum * cnum0
         call stop_clock('norm_dist')
      endif

      !
      !output f(t_i+1) if needed.
      if( mod(istep, output_nstep) .eq. 0 ) then
         if(ionode) then
            call progressbar(istep, boltz_nstep)
            call output_dist(kg, group_id, dist1, (istep/output_nstep) )
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

   enddo

   call stop_clock('dynamics_tot')

   ! close files
   if(ionode) then

      call hdf_close_group(group_id)
      call hdf_close_file(file_id)

      if (boltz_acc_thr > drift_accel_tolerance) close(vel_accel_unit)

      if (pump_pulse) then
         call dump_pump_pulse_data(pump_data)
         call hdf_close_group(pump_data%pump_dist_gid)
         call hdf_close_file(pump_data%pump_fid)
      endif

   endif

   ! release work space
   if(trim(solver) .eq. 'rk4') deallocate(dist_t1, dist_t2)

   deallocate(dist0, dist1)

   if(allocated(deriv)) deallocate(deriv)
   if(allocated(ene)) deallocate(ene)
   if(allocated(popu)) deallocate(popu)
   if(allocated(vel_time)) deallocate(vel_time)
   if(allocated(accel_time)) deallocate(accel_time)

   return

   contains

   !> Setup the carrier dynamics simulation:
   !! - read the initial distribution from the file
   !! - initialize the distribution function
   !! - write the initial distribution to the HDF5 file
   subroutine cdyna_setup(filename, kgrid, fid, gid, dist)
      use pert_const, only: dp, timeunit, ryd2ev
      use qe_mpi_mod, only: ionode, mp_barrier, mp_bcast, ionode_id, inter_pool_comm
      use boltz_dynamics_mod, only: output_dist, restart_dist, &
         init_dist_fermi, init_dist_lorentz, init_dist_gaussian
      use pert_param, only: boltz_nstep, output_nstep, time_step, &
         boltz_init_dist, boltz_init_e0, boltz_init_smear, boltz_init_ampl, &
         hole
      use boltz_grid, only: grid
      use hdf5_utils
      implicit none
      character(*), intent(in) :: filename                  !< Name of the HDF5 file: prefix_cdyna.h5
      type(grid), intent(in) :: kgrid                       !< k-grid object
      integer(HID_T), intent(out) :: fid                    !< HDF5 file ID for prefix_cdyna.h5
      integer(HID_T), intent(out) :: gid                    !< HDF5 group ID
      real(dp), intent(out) :: dist( kgrid%numb, kgrid%nk ) !< Distribution function
      !local variables
      logical :: has_file
      real(dp) :: t_step_fs                                 !< Time step in fs
      character(len=120) :: group_name, msg
      integer :: nrun                                       !< Number of dynamics runs
      integer :: current_run                                !< Current dynamics run number
      real(dp) :: efield_norm2                              !< Norm of the external electric field
      character(len=6), external :: int_to_char

      gid = 0;  fid = 0  ! init
      ! restart dynamics simulation from previous run
      if( (trim(boltz_init_dist) .eq. 'restart')) then
         inquire(file=trim(filename), exist=has_file)
         if(.not. has_file) call errore('cdyna_setup','missing '// trim(filename), 1)

         ! read only on the ionode
         if(ionode) then
            call hdf_open_file(fid, trim(filename), status='OLD', action='READWRITE')
            call hdf_read_dataset(fid, 'num_runs', nrun)

            call restart_dist(kgrid, fid, nrun, dist)

            ! increment the number of runs
            current_run = nrun + 1
            call hdf_update_dataset(fid, 'num_runs', current_run)
         endif
         ! distribute to all nodes
         call mp_bcast(dist, ionode_id, inter_pool_comm)

      ! initialize the distribution function with one of the functions below
      else
         select case (boltz_init_dist)
         case('fermi')
            call init_dist_fermi(kgrid, dist, boltz_init_e0, boltz_init_smear)
         case ('lorentz')
            call init_dist_lorentz(kgrid, dist, boltz_init_e0, boltz_init_smear, hole)
         case('gaussian')
            call init_dist_gaussian(kgrid, dist, boltz_init_e0, boltz_init_smear, boltz_init_ampl, hole)
         case default
            msg = "invalid boltz_init_dist, valid options are 'restart', 'fermi', 'lorentz', 'gaussian'."
            call errore('carrier_dynamics_run', trim(msg), 1)
         end select

         if(ionode) then
            current_run = 1
            call hdf_open_file(fid, trim(filename), status='NEW')
            call hdf_write_dataset(fid, 'num_runs', current_run)
            !write kg%enk(:,:)
            call hdf_write_dataset(fid, 'band_structure_ryd', kgrid%enk)
            call hdf_write_attribute(fid, 'band_structure_ryd', 'ryd2ev', ryd2ev)
         endif
      endif

      !record status of the current run
      if(ionode) then
         group_name = "dynamics_run_" // trim( int_to_char(current_run) )
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

         if(efield_norm2 > efield_tolerance) call hdf_write_dataset(gid, 'efield',  &
                     boltz_efield / (bohr2ang*1E-8_dp) * ryd2ev )

         !
         if(current_run .eq. 1)  call output_dist(kgrid, gid, dist, 0)
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

