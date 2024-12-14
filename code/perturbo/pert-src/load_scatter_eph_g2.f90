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

subroutine load_scatter_eph_g2(nc_loc, kq_pair, num_kq_pair, scat, tsmear, kg, nmodes)
   use hdf5_utils
   use boltz_utils, only: idx2band
   use qe_mpi_mod, only: ionode, stdout
   use pert_param, only: prefix, tmp_dir, ph_mode_exclude_ranges, adapt_smear_eph
   implicit none
   integer, intent(in) :: nc_loc, num_kq_pair, nmodes
   type(channel), intent(inout) :: kq_pair(nc_loc)
   type(t_scatt), intent(out) :: scat(num_kq_pair)
   type(t_scatt_smear), allocatable, intent(out) :: tsmear(:) 
   type(grid), intent(in) :: kg
   !
   integer :: n, i, ik, nkq, j, tot_chnl, nchl, it, ib
   integer, allocatable :: ist(:), bnd_idx(:)
   real(dp), allocatable :: eph_g2(:)
   ! adaptive smear
   real(dp), allocatable :: adsmear(:)
   !
   logical :: has_file
   integer(HID_T) :: file_id
   character(len=120) :: fname, dset_name
   character(len=6), external :: int_to_char
   ! Exclude ph mode variables
   integer :: im, jm, im1, im2, iscat, ns
   integer :: bands(3)
   logical :: include_phonon_modes(nmodes)
   
   call start_clock("load_eph_g2")
   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_eph_g2_p" 
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   ! check if file exist
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('load_scatter_eph_g2', "Missing file: "//trim(fname), 1)
   ! open file
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')

   n = 0
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      nkq = kq_pair(i)%nkq
      !skip if nkq = 0
      if(nkq < 1) cycle
      
      allocate( ist( nkq ) )
      ist(1) = 0 
      do j = 1, nkq-1
         ist(j+1) = ist(j) + kq_pair(i)%nchl(j)
      enddo

      tot_chnl = sum( kq_pair(i)%nchl(:) )
      allocate( bnd_idx(tot_chnl), eph_g2(tot_chnl) )
      bnd_idx = 0;   eph_g2 = 0.0_dp
      if ( adapt_smear_eph ) then
         allocate( adsmear(tot_chnl) )
         adsmear = 0.0_dp
      endif
      
      !read bnd_idx and eph_g2 from hdf5 file
      dset_name = "bands_index_" // trim( int_to_char(i) )
      call hdf_read_dataset(file_id, trim(dset_name), bnd_idx)
      dset_name = "eph_g2_" // trim( int_to_char(i) )
      call hdf_read_dataset(file_id, trim(dset_name), eph_g2)
      if ( adapt_smear_eph ) then
         dset_name = "adapt_smear_" // trim( int_to_char(i) )
         call hdf_read_dataset(file_id, trim(dset_name), adsmear)
      endif

      do j = 1, nkq
         n = n + 1
         !
         scat(n)%ik = ik
         !scat(n)%iq = kq_pair(i)%iq(j)
         scat(n)%ikq = kq_pair(i)%ikq(j)
         !
         nchl = kq_pair(i)%nchl(j)
         scat(n)%nchl = nchl
         allocate( scat(n)%bands_index(nchl), scat(n)%eph_g2(nchl) )

         !if ( adapt_smear_eph ) allocate( scat(n)%adsmear(nchl) )
         if (adapt_smear_eph ) allocate( tsmear(n)%adsmear(nchl) )
         it = ist(j)
         do ib = 1, nchl
            it = it + 1
            scat(n)%bands_index(ib) = bnd_idx(it)
            scat(n)%eph_g2(ib) = eph_g2(it)
            if ( adapt_smear_eph ) then
               tsmear(n)%adsmear(ib) = adsmear(it)
               !scat(n)%adsmear(ib) = adsmear(it)
            endif
         enddo
      enddo
      deallocate(ist, bnd_idx, eph_g2)
      if ( adapt_smear_eph ) deallocate( adsmear )
      !
   enddo
   !

   call hdf_close_file(file_id)
   call mp_barrier( inter_pool_comm )
   call stop_clock("load_eph_g2")

   ! Set eph_g2 to zero for the phonon modes to exclude
   ! This is done in separate nested do loops to not to effect
   ! the efficiency for the case where no phonon modes have to be excluded.
   if( any(ph_mode_exclude_ranges(:, :) > 0) ) then
   
      call start_clock("exclude_ph_modes")

      write(stdout, '(/, 5x, a)') '>excluding phonon modes:'

      call interpret_ph_exclude_ranges(nmodes, include_phonon_modes)

      ! Iterate over the scat array
!$omp parallel do schedule(guided) default(shared) private(iscat, &
!$omp&  ns, im, bands)
      do iscat = 1, num_kq_pair
         do ns = 1, scat(iscat)%nchl
            bands = idx2band( scat(iscat)%bands_index(ns), kg%numb )
            im = bands(3)
            ! Set the g2 elements to zero for the phonon modes to exclude
            if(.not. include_phonon_modes(im)) &
               scat(iscat)%eph_g2(ns) = 0.0_dp
         end do
      end do
!$omp end parallel do

      call stop_clock("exclude_ph_modes")
   end if

end subroutine load_scatter_eph_g2

