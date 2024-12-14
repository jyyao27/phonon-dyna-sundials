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

module boltz_dynamics_mod
   use pert_const, only: dp, ryd2ev, timeunit, bohr2ang
   use pert_data,  only: alat
   use boltz_grid, only: grid
   use pert_utils, only: find_free_unit, bose
   use qe_mpi_mod, only: ionode, stdout, mp_split_pools, mp_sum, inter_pool_comm
   use hdf5_utils
   implicit none
   private
   character(len=6), external :: int_to_char

   public :: restart_dist, init_dist_lorentz, init_dist_gaussian, init_dist_fermi, &
      output_dist, read_dist, calc_carrier_population, output_carrier_population, &
      calc_mean_velocity, calc_mean_velocity_integral, get_energy_array, &
      output_stop_dyna_message, output_vel_accel, print_header_vel_accel, &
      init_disp_bose, restart_dist_full, calc_carrier_energy, init_dist_fermi_separate
contains

subroutine restart_dist(kgrid, file_id, last_run, dist)
   implicit none
   type(grid), intent(in) :: kgrid
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: last_run
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk)
   !
   integer :: last_t
   integer(HID_T) :: group_id
   character(len=120) :: group_name, dname
   
   group_name = "dynamics_run_" // trim( int_to_char(last_run) )
   call hdf_open_group(file_id, trim(group_name), group_id)
   !  
   if(.not. hdf_exists(group_id, 'num_steps')) call errore('restart_dist', &
      'Restart from previous run fails, the previous run is unsuccessful.', 1)
   
   call hdf_read_dataset(group_id, 'num_steps', last_t)
   dname = "snap_t_" // trim( int_to_char(last_t) )
   !
   if(.not. hdf_exists(group_id, trim(dname))) call errore('restart_dist', &
      'Restart from previous run fails, the previous run is unsuccessful.', 1)

   call read_dist(kgrid, group_id, dist, last_t)
   
   call hdf_close_group(group_id)
end subroutine restart_dist


subroutine restart_dist_full(kgrid, qgrid, file_id, last_run, dist, qdist)
   implicit none
   type(grid), intent(in) :: kgrid, qgrid
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: last_run
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk), qdist(qgrid%numb, qgrid%nk)
   !
   integer :: last_t
   integer(HID_T) :: group_id
   character(len=120) :: group_name, dname

   group_name = "carrier_dynamics_run_" // trim( int_to_char(last_run) )
   call hdf_open_group(file_id, trim(group_name), group_id)
   !  
   if(.not. hdf_exists(group_id, 'num_steps')) call errore('restart_dist', &
      'Restart from previous run fails, the previous run is unsuccessful.', 1)

   call hdf_read_dataset(group_id, 'num_steps', last_t)
   dname = "snap_t_" // trim( int_to_char(last_t) )
   !
   if(.not. hdf_exists(group_id, trim(dname))) call errore('restart_dist', &
      'Restart from previous run fails, the previous run is unsuccessful.', 1)

   call read_dist(kgrid, group_id, dist, last_t)

   call hdf_close_group(group_id)

   group_name = "phonon_dynamics_run_" // trim( int_to_char(last_run) )
   call hdf_open_group(file_id, trim(group_name), group_id)
   !  
   if(.not. hdf_exists(group_id, 'num_steps')) call errore('restart_dist', &
      'Restart from previous run fails, the previous run is unsuccessful.', 1)

   call hdf_read_dataset(group_id, 'num_steps', last_t)
   dname = "snap_t_" // trim( int_to_char(last_t) )
   !
   if(.not. hdf_exists(group_id, trim(dname))) call errore('restart_dist', &
      'Restart from previous run fails, the previous run is unsuccessful.', 1)

   call read_dist(qgrid, group_id, qdist, last_t)

   call hdf_close_group(group_id)   
end subroutine restart_dist_full


subroutine init_disp_bose(qgrid, qdist, qtemp)
   use pert_param, only: phfreq_cutoff_ph
   implicit none
   type(grid), intent(in) :: qgrid
   real(dp), intent(in) :: qtemp
   real(dp), intent(out) :: qdist(qgrid%numb, qgrid%nk)
   !local variables
   integer :: ik, ib
   real(dp) :: wq

   do ik = 1, qgrid%nk
      do ib = 1, qgrid%numb
         wq = qgrid%enk(ib, ik)
         qdist(ib, ik) = bose(qtemp, wq)
         if (wq < phfreq_cutoff_ph) qdist(ib,ik) = 0.0_dp
      enddo
   enddo
end subroutine init_disp_bose

subroutine init_dist_fermi_separate(kgrid, dist, midline, efe, efh, smear)
   use pert_utils, only: fermi
   implicit none
   type(grid), intent(in) :: kgrid
   real(dp), intent(in) :: efe, efh, midline, smear
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk)
   !local variables
   integer :: ik, ib
   real(dp) :: enk
   
   do ik = 1, kgrid%nk
      do ib = 1, kgrid%numb
         if (kgrid%enk(ib, ik) > midline) then
            ! for electrons
            enk = kgrid%enk(ib, ik) - efe
         else
            ! for holes
            enk = kgrid%enk(ib, ik) - efh
         endif
         dist(ib, ik) = fermi(smear, enk)
      enddo
   enddo
end subroutine init_dist_fermi_separate

subroutine init_dist_fermi(kgrid, dist, e0, smear)
   use pert_utils, only: fermi
   implicit none
   type(grid), intent(in) :: kgrid
   real(dp), intent(in) :: e0, smear
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk)
   !local variables
   integer :: ik, ib
   real(dp) :: enk
   
   do ik = 1, kgrid%nk
      do ib = 1, kgrid%numb
         enk = kgrid%enk(ib, ik) - e0
         dist(ib, ik) = fermi(smear, enk)
      enddo
   enddo
end subroutine init_dist_fermi


subroutine init_dist_lorentz(kgrid, dist, e0, smear, lhole)
   implicit none
   type(grid), intent(in) :: kgrid
   real(dp), intent(in) :: e0, smear
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk)
   !
   logical,  intent(in), optional :: lhole
   !local variables
   integer :: ik, ib
   real(dp) :: enk
   
   do ik = 1, kgrid%nk
      do ib = 1, kgrid%numb
         enk = (kgrid%enk(ib, ik) - e0) / smear
         dist(ib, ik) = 0.1_dp / (enk*enk + 1.0_dp)
      enddo
   enddo
   !in case of hole carrier
   if( present(lhole) ) then
      if(lhole)  dist(:,:) = 1.0E0_dp - dist(:,:)
   endif
end subroutine init_dist_lorentz

!> Initialize the electron occupation function with
!! the Guassian distribution around a given energy.
subroutine init_dist_gaussian(kgrid, dist, e0, smear, lhole)
   implicit none
   type(grid), intent(in) :: kgrid                     !< k-point grid
   real(dp), intent(in) :: e0                          !< Energy center of the Gaussian
   real(dp), intent(in) :: smear                       !< Guassian smearing
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk) !< Occupation function
   !
   logical,  intent(in), optional :: lhole
   !local variables
   integer :: ik, ib
   real(dp) :: enk
   
   do ik = 1, kgrid%nk
      do ib = 1, kgrid%numb
         enk = (kgrid%enk(ib, ik) - e0) / smear
         dist(ib, ik) = 0.1_dp*exp( -min(enk*enk, 200.0_dp) )
      enddo
   enddo
   !in case of hole carrier
   if( present(lhole) ) then
      if(lhole)  dist(:,:) = 1.0E0_dp - dist(:,:)
   endif
end subroutine init_dist_gaussian


subroutine output_dist(kgrid, group_id, dist, irec)
   implicit none
   type(grid), intent(in) :: kgrid
   integer(HID_T), intent(in) :: group_id
   real(dp), intent(in) :: dist(kgrid%numb, kgrid%nk)
   integer, intent(in) :: irec
   !
   character(len=120) :: dset_name
   
   dset_name = "snap_t_" // trim( int_to_char(irec) )
   call hdf_write_dataset(group_id, trim(dset_name), dist)
   call hdf_update_dataset(group_id, 'num_steps', irec)
end subroutine output_dist


subroutine output_stop_dyna_message(nstep, drift_vel, drift_accel, vel_time, accel_time)
   use qe_mpi_mod, only: stdout
   use yaml_utils, only: ymlout
   implicit none
   integer, intent(in) :: nstep
   real(dp), intent(in) :: vel_time(:), accel_time(:)
   real(dp), intent(in) :: drift_vel, drift_accel
   !
   integer :: it

   write(ymlout, '(/,3x,a,i7)') 'number of steps: ', nstep
   
   write(stdout,'(5x, a, /, 5x, a, i10, 1x, a, /, 5x, a, 1x, f18.8, /, 5x, a, 1x, f18.8)') &
      'Drift acceleration in 10 last steps is lower than the threshold.', &
      'Dynamics ended in', nstep, 'steps.', &
      'Drift velocity (cm/s):', drift_vel, &
      'Drift acceleration (cm/s2):', drift_accel


   write(ymlout, '(/,3x,a)') 'velocity:'
   do it = 1, nstep
      write(ymlout,'(6x,"-", 1x, es23.16)') vel_time(it)
   enddo

   write(ymlout, '(/,3x,a)') 'acceleration:'
   do it = 1, nstep
      write(ymlout,'(6x,"-", 1x, es23.16)') accel_time(it)
   enddo

end subroutine output_stop_dyna_message


subroutine print_header_vel_accel(vel_accel_unit)
   use pert_param, only: prefix
   use yaml_utils, only: ymlout
   implicit none
   integer, intent(out) :: vel_accel_unit
   !
   character(len=120) :: vel_accel_name

   ! open the file for drift vel. and accel.
   vel_accel_name = trim(prefix)//'_vel_accel.dat'

   open(newunit = vel_accel_unit, file=vel_accel_name, status='unknown',action='write')
   write(vel_accel_unit,'(a)') ' # time (fs)    drift velocity (cm/s)    drift acceleration (cm/s2)'

   write(ymlout, '(/,3x,a)') 'velocity units: cm/s'
   write(ymlout, '(/,3x,a)') 'acceleration units: cm/s2'

end subroutine print_header_vel_accel


subroutine output_vel_accel(vel_accel_unit, time_cur, drift_vel, drift_accel)
   implicit none
   integer, intent(in) :: vel_accel_unit
   real(dp), intent(in) :: time_cur, drift_vel, drift_accel
   
   write(vel_accel_unit,'(1x,f16.4,2x,es23.16,2x,es23.16)') &
      time_cur, drift_vel, drift_accel

end subroutine output_vel_accel


subroutine read_dist(kgrid, group_id, dist, irec)
   implicit none
   type(grid), intent(in) :: kgrid
   integer(HID_T), intent(in) :: group_id
   real(dp), intent(out)  :: dist(kgrid%numb, kgrid%nk)
   integer, intent(in)  :: irec
   !local variables
   character(len=120) :: dset_name
   
   dset_name = "snap_t_" // trim( int_to_char(irec) )
   
   call hdf_read_dataset(group_id, trim(dset_name), dist)
end subroutine read_dist


subroutine calc_carrier_population(kgrid, dist, ene, popu, lhole, cnum)
   use pert_param, only: boltz_de
   implicit none
   type(grid), intent(in) :: kgrid
   real(dp), intent(in)  :: dist(kgrid%numb, kgrid%nk)
   real(dp), intent(in)  :: ene(:)
   real(dp), intent(out) :: popu(:)
   !
   logical,  intent(in), optional :: lhole
   real(dp), intent(out), optional :: cnum
   
   logical :: lh
   integer :: nstep, it_st, it_end, it, ic, ie, ib
   real(dp) :: et(kgrid%numb, 4), fnk(kgrid%numb, 4), dtmp
   real(dp), allocatable :: pop_t(:)

   nstep = size(ene)
   popu(1:nstep) = 0.0_dp
   !
   if(present(lhole)) then
      lh = lhole
   else
      lh = .false.
   endif

   call mp_split_pools(kgrid%num_tet, it_st, it_end)
!$omp parallel default(shared) private(it, ic, et, fnk, pop_t, ie, ib, dtmp)
   allocate( pop_t(nstep) ); pop_t = 0.0_dp
!$omp do schedule(guided)
   do it = it_st, it_end
      do ic = 1, 4
         do ib = 1, kgrid%numb
            et(ib, ic)  = kgrid%enk(ib, kgrid%tetra(2,ic,it) )
            dtmp = dist(ib, kgrid%tetra(2,ic,it))
            fnk(ib, ic) = merge(1.0E0_dp - dtmp, dtmp, lh)
         enddo
      enddo
      call tetra_int(kgrid%numb, et, fnk, nstep, ene, pop_t)
   enddo
!$omp end do nowait

   do ie = 1, nstep
!$omp atomic
      popu(ie) = popu(ie) + pop_t(ie)
   enddo
   deallocate(pop_t)
!$omp end parallel

   !! proper weight, no spin degeneracy for now.  popu per unit cell.
   ! tweight  = V_T / V_G, 1 / (total number of tetrahedron)
   popu(1:nstep) = popu(1:nstep) * kgrid%tweight
   ! collect results from different pools
   call mp_sum( popu(1:nstep), inter_pool_comm )

   if( present(cnum) ) then
      cnum = sum(popu) * boltz_de
   endif

end subroutine calc_carrier_population


subroutine calc_carrier_energy(kgrid, dist, ene, de, total_e, lhole)
   implicit none
   type(grid), intent(in) :: kgrid
   real(dp), intent(in)  :: dist(kgrid%numb, kgrid%nk)
   real(dp), intent(in)  :: de
   real(dp), intent(in)  :: ene(:)
   real(dp), intent(out) :: total_e
   !
   logical,  intent(in), optional :: lhole
   
   logical :: lh
   integer :: nstep, it_st, it_end, it, ic, ie, ib
   real(dp) :: et(kgrid%numb, 4), fenk(kgrid%numb, 4), dtmp
   real(dp), allocatable :: popu(:)
   real(dp), allocatable :: pop_t(:)

   nstep = size(ene)
   allocate(popu(nstep))
   popu(1:nstep) = 0.0_dp
   total_e = 0.0_dp
   !
   if(present(lhole)) then
      lh = lhole
   else
      lh = .false.
   endif

   call mp_split_pools(kgrid%num_tet, it_st, it_end)
!$omp parallel default(shared) private(it, ic, et, fenk, pop_t, ie, ib, dtmp)
   allocate( pop_t(nstep) ); pop_t = 0.0_dp
!$omp do schedule(guided)
   do it = it_st, it_end
      do ic = 1, 4
         do ib = 1, kgrid%numb
            et(ib, ic)  = kgrid%enk(ib, kgrid%tetra(2,ic,it) )
            dtmp = dist(ib, kgrid%tetra(2,ic,it))
            fenk(ib, ic) = merge(1.0E0_dp - dtmp, dtmp, lh) * et(ib, ic)
         enddo
      enddo
      call tetra_int(kgrid%numb, et, fenk, nstep, ene, pop_t)
   enddo
!$omp end do nowait

   do ie = 1, nstep
!$omp atomic
      popu(ie) = popu(ie) + pop_t(ie)
   enddo
   deallocate(pop_t)
!$omp end parallel

   !! proper weight, no spin degeneracy for now.  popu per unit cell.
   ! tweight  = V_T / V_G, 1 / (total number of tetrahedron)
   popu(1:nstep) = popu(1:nstep) * kgrid%tweight
   ! collect results from different pools
   call mp_sum( popu(1:nstep), inter_pool_comm )

   total_e = sum(popu) * de

end subroutine calc_carrier_energy


function get_energy_array(kgrid) result(energy_array)
   use pert_param, only: boltz_emin, boltz_emax, boltz_de
   implicit none
   type(grid), intent(in) :: kgrid
   !
   real(dp), allocatable :: energy_array(:)
   !
   integer :: nestep, i
   real(dp) :: emin, emax

   emin = max(boltz_emin, minval(kgrid%enk))
   emax = min(boltz_emax, maxval(kgrid%enk))
   nestep = int( floor((emax-emin)/boltz_de) ) + 3

   allocate(energy_array(nestep))

   do i = 1, nestep
      energy_array(i) = emin + (i-2)*boltz_de
   enddo


end function get_energy_array


function calc_mean_velocity_integral(kgrid, dist, ene, efield, lhole) result(vel_integral)
   use pert_param, only: boltz_de
   implicit none
   type(grid), intent(in) :: kgrid
   real(dp), intent(in)  :: dist(kgrid%numb, kgrid%nk)
   real(dp), intent(in)  :: ene(:), efield(3)
   !
   real(dp) :: vel_integral
   !
   logical,  intent(in), optional :: lhole
   
   logical :: lh
   integer :: nstep, it_st, it_end, it, ic, ie, ib
   real(dp) :: cnum, estep
   real(dp) :: et(kgrid%numb, 4), fnk(kgrid%numb, 4), vnk(kgrid%numb, 4), evec(3), dtmp
   real(dp), allocatable :: vel_t(:), pop_t(:), popu(:), vel(:)
   
   if(present(lhole)) then
      lh = lhole
   else
      lh = .false.
   endif

   nstep = size(ene)

   estep = (ene(nstep) - ene(1)) / real(nstep - 1, dp)

   allocate(vel(nstep))
   allocate(popu(nstep))

   popu(:) = 0.0_dp
   vel(:) = 0.0_dp

   ! E/|E|
   evec(1:3) = efield(1:3) / sqrt(dot_product(efield(1:3),efield(1:3)))

   call mp_split_pools(kgrid%num_tet, it_st, it_end)
!$omp parallel default(shared) private(it, ic, et, fnk, vnk, pop_t, vel_t, ie, ib, dtmp)
   allocate( pop_t(nstep) ); pop_t = 0.0_dp
   allocate( vel_t(nstep) ); vel_t = 0.0_dp
!$omp do schedule(guided)
   do it = it_st, it_end
      do ic = 1, 4
         do ib = 1, kgrid%numb
            et(ib, ic) = kgrid%enk(ib, kgrid%tetra(2,ic,it) )
            dtmp = dist(ib, kgrid%tetra(2,ic,it) )
            fnk(ib, ic) = merge(1.0E0_dp - dtmp, dtmp, lh)
            vnk(ib, ic) = merge(1.0E0_dp - dtmp, dtmp, lh) * &
               dot_product(kgrid%vnk(:, ib, kgrid%tetra(2,ic,it)), evec(:))
         enddo
      enddo

      call tetra_int(kgrid%numb, et, fnk, nstep, ene, pop_t)
      call tetra_int(kgrid%numb, et, vnk, nstep, ene, vel_t)

   enddo
!$omp end do nowait

   do ie = 1, nstep
!$omp atomic
      popu(ie) = popu(ie) + pop_t(ie)
!$omp atomic
      vel(ie) = vel(ie) + vel_t(ie)
   enddo

   deallocate(vel_t)
   deallocate(pop_t)
!$omp end parallel

   !! proper weight, no spin degeneracy for now.  popu per unit cell.
   ! tweight  = V_T / V_G, 1 / (total number of tetrahedron)
   ! "alat*bohr2ang": convert from Ryd*Bohr*alat/hbar to Ryd*angstrom/hbar
   ! timeunit t0=hbar/Ryd, so the unit is equal to Angstrom/t0
   ! and then "1.0E-8_dp/timeunit", convert Angstrom/t0 to cm/s
   vel(1:nstep) = vel(1:nstep) * (kgrid%tweight * alat*bohr2ang*1.0E-8_dp/timeunit)

   !! proper weight, no spin degeneracy for now.  popu per unit cell.
   ! tweight  = V_T / V_G, 1 / (total number of tetrahedron)
   popu(1:nstep) = popu(1:nstep) * kgrid%tweight

   ! collect results from different pools
   call mp_sum( vel(1:nstep), inter_pool_comm )
   call mp_sum( popu(1:nstep), inter_pool_comm )

   ! number of carriers
   cnum  = sum(popu) * boltz_de

   vel_integral = sum(vel) * estep / cnum

   !write(*,*) 'INSIDE1', vel(4:6), sum(vel)

end function calc_mean_velocity_integral


subroutine calc_mean_velocity(kgrid, dist, ene, vel, efield, lhole)
   implicit none
   type(grid), intent(in) :: kgrid
   real(dp), intent(in)  :: dist(kgrid%numb, kgrid%nk)
   real(dp), intent(in)  :: ene(:), efield(3)
   real(dp), intent(out) :: vel(:) ! <v> * E/|E|, mean velocity along E-directio
   !
   logical,  intent(in), optional :: lhole
   
   logical :: lh
   integer :: nstep, it_st, it_end, it, ic, ie, ib
   real(dp) :: et(kgrid%numb, 4), fnk(kgrid%numb, 4), evec(3), dtmp
   real(dp), allocatable :: vel_t(:)
   
   if(present(lhole)) then
      lh = lhole
   else
      lh = .false.
   endif

   nstep = size(ene)
   vel(1:nstep) = 0.0_dp
   ! E/|E|
   evec(1:3) = efield(1:3) / sqrt(dot_product(efield(1:3),efield(1:3)))

   call mp_split_pools(kgrid%num_tet, it_st, it_end)
!$omp parallel default(shared) private(it, ic, et, fnk, vel_t, ie, ib, dtmp)
   allocate( vel_t(nstep) ); vel_t = 0.0_dp
!$omp do schedule(guided)
   do it = it_st, it_end
      do ic = 1, 4
         do ib = 1, kgrid%numb
            et(ib, ic) = kgrid%enk(ib, kgrid%tetra(2,ic,it) )
            dtmp = dist(ib, kgrid%tetra(2,ic,it) )
            fnk(ib,ic) = merge(1.0E0_dp - dtmp, dtmp, lh) * &
               dot_product(kgrid%vnk(:, ib, kgrid%tetra(2,ic,it)), evec(:))
         enddo
      enddo
      call tetra_int(kgrid%numb, et, fnk, nstep, ene, vel_t)
   enddo
!$omp end do nowait

   do ie = 1, nstep
!$omp atomic
      vel(ie) = vel(ie) + vel_t(ie)
   enddo
   deallocate(vel_t)
!$omp end parallel

   !! proper weight, no spin degeneracy for now.  popu per unit cell.
   ! tweight  = V_T / V_G, 1 / (total number of tetrahedron)
   ! "alat*bohr2ang": convert from Ryd*Bohr*alat/hbar to Ryd*angstrom/hbar
   ! timeunit t0=hbar/Ryd, so the unit is equal to Angstrom/t0
   ! and then "1.0E-8_dp/timeunit", convert Angstrom/t0 to cm/s
   vel(1:nstep) = vel(1:nstep) * (kgrid%tweight * alat*bohr2ang*1.0E-8_dp/timeunit)
   ! collect results from different pools
   call mp_sum( vel(1:nstep), inter_pool_comm )
end subroutine calc_mean_velocity


subroutine output_carrier_population(ene, popu, mean_vel, punit, irec, fname)
   implicit none
   real(dp), intent(in) :: ene(:), popu(:), mean_vel
   integer, intent(inout)  :: punit
   integer, intent(in)  :: irec
   character(len=*), intent(in), optional :: fname
   
   integer :: i, nsize
   real(dp) :: cnum, estep

   nsize = size(ene)
   estep = (ene(nsize)-ene(1)) / real(nsize-1, dp)
   if(irec .eq. 0) then
      if(present(fname)) then
         punit = find_free_unit()
         open(punit, file=fname, form='formatted', status='unknown',action='write')
      else
         close(punit)
         return
      endif
   else
      cnum = sum(popu(1:nsize))*estep
      write(punit, '(a)') ''
      write(punit, 1001)  irec, cnum, mean_vel/cnum
      do i = 1, nsize
         write(punit,'(f10.5, 2x, es23.16)') ene(i)*ryd2ev, popu(i)
      enddo
      write(punit, '(a)') ''
   endif

1001 format(1x,'#.record ', i8, 2x, ' #.carrier/u.c. ', es23.16, ' #.mean.vel. ', es23.16)
end subroutine output_carrier_population

end module boltz_dynamics_mod
