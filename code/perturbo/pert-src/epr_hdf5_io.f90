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
!   read and write prefix_epr.hdf5 file
!
! Maintenance:
!===============================================================================

module epr_hdf5_io
   use qe_mpi_mod, only: stdout
   use kinds, only: dp
   use force_constant, only: lattice_ifc
   use force_constant_thirdorder, only: lattice_ifct
   use electron_wannier, only: electron_wann
   use elph_matrix_wannier, only: elph_mat_wann, eph_wannier
   use hdf5_utils
   implicit none
   character(len=6), external :: int_to_char

   public :: open_epr_file
   public :: write_force_constant, read_force_constant
   public :: write_force_constant_thirdorder, read_force_constant_thirdorder
   public :: write_electron_wannier, read_electron_wannier
   public :: write_elph_mat_wann, read_elph_mat_wann
   public :: write_elph_mat_wann_part, read_elph_mat_wann_part
   public :: write_ephmat_part, read_ephmat_part
contains

!> Open the file containing the e-ph elements in Wannier basis
!! called <prefix>_epr.h5
subroutine open_epr_file(prefix, epr_fid)
   implicit none
   character(len=*), intent(in) :: prefix !< prefix of the calculation
   integer(HID_T), intent(out) :: epr_fid !< HDF5-formatter ID of the file

   character(len=50) :: epr_filename
   logical :: has_epr_file, has_epwan_file

   epr_fid = 0

   ! <1> Look for the prefix_epr.h5 file
   epr_filename = trim( trim(prefix)//"_epr.h5" )

   inquire(file=epr_filename, exist=has_epr_file)

   if (.not. has_epr_file) then
      epr_filename = trim( trim(prefix)//"_epwan.h5" )

      ! <2> Look for the prefix_epwan.h5 file
      inquire(file=epr_filename, exist=has_epwan_file)

      if(has_epwan_file) then
         ! In the 3.0 release: call errore() instead of the warning
         write(stdout,'(/,5x,a)') repeat("!", 50)
         write(stdout,'(5x,a)') "Warning: "//trim(prefix)//"_epr.h5 file was not found in the directory."
         write(stdout,'(5x,a)') "However, "//trim(prefix)//"_epwan.h5 file was found. The "&
            //trim(epr_filename)//" will be used in the calculation."
         write(stdout,'(5x,a)') "Please note, <prefix>_epr.h5 name is now the standard name for the HDF5 file"
         write(stdout,'(5x,a)') "storing the e-ph elements in the Wannier basis."
         write(stdout,'(5x,a)') "The <prefix>_epwan.h5 filename is deprecated and will not be supported in the future release."
         write(stdout,'(5x,a,/)') repeat("!", 50)

      else
         call errore('perturbo','missing '// trim(prefix)// "_epr.h5 file!", 1)
      endif

   endif

   write(stdout,'(5x,a,/)') "The "//trim(epr_filename)//" HDF5 file was found and will be used in the calculation."
   call hdf_open_file(epr_fid, epr_filename, status='OLD', action='READ')

end subroutine open_epr_file


subroutine write_force_constant(file_id, phon)
   implicit none
   integer(HID_T), intent(in) :: file_id !< prefix of the calculation
   type(lattice_ifc), intent(in) :: phon !< HDF5-formatter ID of the file
   !local
   integer :: m, nelem
   integer(HID_T) :: group_id
   character(len=120) :: dset_name
   
   nelem = phon%na * (phon%na + 1) / 2
   !write force constant
   call hdf_create_group(file_id, 'force_constant')
   call hdf_open_group(file_id, 'force_constant', group_id)
   do m = 1, nelem
      dset_name = 'ifc'// trim( int_to_char(m) )
      call hdf_write_dataset( group_id, trim(dset_name), real(phon%phi(m)%ifc) )
   enddo
   call hdf_close_group(group_id)
end subroutine write_force_constant


subroutine write_force_constant_thirdorder(file_id, phont)
   implicit none
   integer(HID_T), intent(in) :: file_id !< prefix of the calculation
   type(lattice_ifct), intent(in) :: phont !< HDF5-formatter ID of the file
   !local
   integer :: m, nelem
   integer(HID_T) :: group_id
   character(len=120) :: dset_name

   !nelem = (phont%na**3+3*phont%na**2+2*phont%na)/6
   nelem = phont%na**3
   !write force constant
   call hdf_create_group(file_id, 'force_constant_thirdorder')
   call hdf_open_group(file_id, 'force_constant_thirdorder', group_id)
   do m = 1, nelem
      dset_name = 'ifct'// trim( int_to_char(m) )
      call hdf_write_dataset( group_id, trim(dset_name), real(phont%phit(m)%ifct) )
   enddo
   call hdf_close_group(group_id)
end subroutine write_force_constant_thirdorder


subroutine read_force_constant(file_id, phon)
   implicit none
   integer(HID_T), intent(in) :: file_id
   type(lattice_ifc), intent(inout) :: phon
   !local
   integer :: m, nelem, nr
   integer(HID_T) :: group_id
   real(dp), allocatable :: r_val(:,:,:)
   character(len=120) :: dset_name
   
   nelem = phon%na * (phon%na + 1) / 2
   allocate( r_val(3, 3, phon%max_nr) )
   !read force constant
   call hdf_open_group(file_id, 'force_constant', group_id)
   do m = 1, nelem
      nr = phon%phi(m)%ws_ph%nr
      !
      dset_name = 'ifc'// trim( int_to_char(m) )
      call hdf_read_dataset(group_id, trim(dset_name), r_val(:,:,1:nr))

      phon%phi(m)%ifc = cmplx(r_val(:,:,1:nr), 0.0_dp, kind=dp)
   enddo
   call hdf_close_group( group_id )
   
   deallocate( r_val )
end subroutine read_force_constant

subroutine read_force_constant_thirdorder(file_id, phont)
   implicit none
   integer(HID_T), intent(in) :: file_id
   type(lattice_ifct), intent(inout) :: phont
   !local
   integer :: m, nelem, nr
   integer(HID_T) :: group_id
   real(dp), allocatable :: r_val(:,:,:,:)
   character(len=120) :: dset_name

   !nelem = (phont%na**3+3*phont%na**2+2*phont%na)/6
   nelem = phont%na**3
   allocate( r_val(3, 3, 3, phont%max_nrt) )
   !read force constant
   call hdf_open_group(file_id, 'force_constant_thirdorder', group_id)
   do m = 1, nelem
      nr = phont%phit(m)%trip_ph%nrt
      !
      dset_name = 'ifct'// trim( int_to_char(m) )
      call hdf_read_dataset(group_id, trim(dset_name), r_val(:,:,:,1:nr))

      phont%phit(m)%ifct = cmplx(r_val(:,:,:,1:nr), 0.0_dp, kind=dp)
   enddo
   call hdf_close_group( group_id )

   deallocate( r_val )
end subroutine read_force_constant_thirdorder

subroutine write_electron_wannier(file_id, elec)
   implicit none
   integer(HID_T), intent(in) :: file_id
   type(electron_wann), intent(in) :: elec
   !local
   integer :: m, nelem
   integer(HID_T) :: group_id
   character(len=120) :: dset_name
   
   nelem = elec%nb * (elec%nb + 1) / 2
   !write electron hamiltonian
   call hdf_create_group(file_id, 'electron_wannier')
   call hdf_open_group(file_id, 'electron_wannier', group_id)
   do m = 1, nelem
      dset_name = 'hopping_r'// trim( int_to_char(m) )
      call hdf_write_dataset(group_id, trim(dset_name), real(elec%ham_r(m)%hop(:)) )
      dset_name = 'hopping_i'// trim( int_to_char(m) )
      call hdf_write_dataset(group_id, trim(dset_name), aimag(elec%ham_r(m)%hop(:)))
   enddo
   call hdf_close_group(group_id)
end subroutine write_electron_wannier


subroutine read_electron_wannier(file_id, elec)
   implicit none
   integer(HID_T), intent(in) :: file_id
   type(electron_wann), intent(inout) :: elec
   !local
   integer :: m, nelem, nr
   integer(HID_T) :: group_id
   real(dp), allocatable :: r_val(:), i_val(:)
   character(len=120) :: dset_name
   
   nelem = elec%nb * (elec%nb + 1) / 2
   !load H(R) data
   allocate( r_val(elec%max_nr),  i_val(elec%max_nr) )
   !
   call hdf_open_group(file_id, 'electron_wannier', group_id)
   do m = 1, nelem
      ! allocate space
      nr = elec%ham_r(m)%ws_el%nr
      !read data
      dset_name = 'hopping_r'// trim( int_to_char(m) )
      call hdf_read_dataset(group_id, trim(dset_name), r_val(1:nr))
      dset_name = 'hopping_i'// trim( int_to_char(m) )
      call hdf_read_dataset(group_id, trim(dset_name), i_val(1:nr))
      
      elec%ham_r(m)%hop(1:nr) = cmplx(r_val(1:nr), i_val(1:nr), kind=dp)
   enddo
   call hdf_close_group(group_id)
   
   deallocate( r_val, i_val )
end subroutine read_electron_wannier


subroutine write_ephmat_part(file_id, iatom, icart, nb, nk, nq, epmat)
   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: iatom, icart, nb, nk, nq
   complex(dp), intent(in) :: epmat(nb, nb, nk, nq)
   !
   character(len=120) :: dset_name

   dset_name = "elph_" // trim(int_to_char(iatom)) &
               // "_" // trim(int_to_char(icart)) // "_r"
   call hdf_write_dataset(file_id, trim(dset_name), real(epmat))

   dset_name = "elph_" // trim(int_to_char(iatom)) &
               // "_" // trim(int_to_char(icart)) // "_i"
   call hdf_write_dataset(file_id, trim(dset_name), aimag(epmat))

end subroutine


subroutine read_ephmat_part(file_id, iatom, icart, nb, nk, nq, epmat)
   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: iatom, icart, nb, nk, nq
   complex(dp), intent(out) :: epmat(nb, nb, nk, nq)
   !
   character(len=120) :: dset_name
   real(dp), allocatable :: r_val(:,:,:,:), i_val(:,:,:,:)
   
   epmat = (0.0_dp, 0.0_dp)
   allocate( r_val(nb, nb, nk, nq), i_val(nb, nb, nk, nq) )
   
   dset_name = "elph_" // trim(int_to_char(iatom)) &
               // "_" // trim(int_to_char(icart)) // "_r"
   call hdf_read_dataset(file_id, trim(dset_name), r_val)
   
   dset_name = "elph_" // trim(int_to_char(iatom)) &
               // "_" // trim(int_to_char(icart)) // "_i"
   call hdf_read_dataset(file_id, trim(dset_name), i_val)
   
   epmat = cmplx(r_val, i_val, kind=dp)
   deallocate(r_val, i_val)
end subroutine


subroutine write_elph_mat_wann_part(group_id, iwan, jwan, iatom, pep)
   implicit none
   integer(HID_T), intent(in) :: group_id
   integer :: iwan, jwan, iatom
   type(eph_wannier), intent(in) :: pep
   !local
   character(len=120) :: dset_name

   dset_name = "ep_hop_r_"// trim(int_to_char(iatom)) &
      //'_'// trim(int_to_char(jwan)) //'_'// trim(int_to_char(iwan))
   call hdf_write_dataset(group_id, trim(dset_name), real(pep%ep_hop))
   
   dset_name = "ep_hop_i_"// trim(int_to_char(iatom)) &
      //'_'// trim(int_to_char(jwan)) //'_'// trim(int_to_char(iwan))
   call hdf_write_dataset(group_id, trim(dset_name), aimag(pep%ep_hop))

end subroutine


subroutine read_elph_mat_wann_part(group_id, iwan, jwan, iatom, pep)
   implicit none
   integer(HID_T), intent(in) :: group_id
   integer :: iwan, jwan, iatom
   type(eph_wannier), intent(inout) :: pep
   !local
   integer :: nre, nrp
   character(len=120) :: dset_name
   real(dp), allocatable :: r_val(:,:,:), i_val(:,:,:)

   nrp = pep%ws_ph%nr
   nre = pep%ws_el%nr
   allocate( r_val(3, nre, nrp), i_val(3, nre, nrp) )
      
   dset_name = "ep_hop_r_"// trim(int_to_char(iatom)) &
      //'_'// trim(int_to_char(jwan)) //'_'// trim(int_to_char(iwan))
   call hdf_read_dataset(group_id, trim(dset_name), r_val )

   dset_name = "ep_hop_i_"// trim(int_to_char(iatom)) &
      //'_'// trim(int_to_char(jwan)) //'_'// trim(int_to_char(iwan))
   call hdf_read_dataset(group_id, trim(dset_name), i_val )

   pep%ep_hop = cmplx(r_val, i_val, kind=dp)
   !
   deallocate(r_val, i_val)
end subroutine


subroutine write_elph_mat_wann(file_id, elph)
   implicit none
   integer(HID_T), intent(in) :: file_id
   type(elph_mat_wann), intent(in) :: elph
   !local
   integer :: ia, jw, iw
   integer(HID_T) :: group_id
   character(len=120) :: dset_name
   type(eph_wannier), pointer :: pep

   call hdf_create_group(file_id, 'eph_matrix_wannier')
   call hdf_open_group(file_id, "eph_matrix_wannier", group_id)

   do ia = 1, elph%na
   do jw = 1, elph%nb
   do iw = 1, elph%nb
      pep => elph%epr(iw, jw, ia)

      dset_name = "ep_hop_r_"// trim(int_to_char(ia)) &
         //'_'// trim(int_to_char(jw)) //'_'// trim(int_to_char(iw))
      call hdf_write_dataset(group_id, trim(dset_name), real(pep%ep_hop))
      
      dset_name = "ep_hop_i_"// trim(int_to_char(ia)) &
         //'_'// trim(int_to_char(jw)) //'_'// trim(int_to_char(iw))
      call hdf_write_dataset(group_id, trim(dset_name), aimag(pep%ep_hop))
      !
   enddo; enddo; enddo
   call hdf_close_group( group_id )  

end subroutine write_elph_mat_wann


subroutine read_elph_mat_wann(file_id, elph)
   implicit none
   integer(HID_T), intent(in) :: file_id
   type(elph_mat_wann), intent(inout) :: elph
   !local
   integer(HID_T) :: group_id
   integer :: ia, jw, iw, nre, nrp
   real(dp), allocatable :: r_val(:,:,:), i_val(:,:,:)
   character(len=120) :: dset_name

   !read eph matrix elements 
   call hdf_open_group(file_id, "eph_matrix_wannier", group_id)
   do ia = 1, elph%na
   do jw = 1, elph%nb
   do iw = 1, elph%nb
      nrp = elph%epr(iw, jw, ia)%ws_ph%nr
      nre = elph%epr(iw, jw, ia)%ws_el%nr
      !
      allocate( r_val(3, nre, nrp), i_val(3, nre, nrp) )
      
      dset_name = "ep_hop_r_"// trim(int_to_char(ia)) &
         //'_'// trim(int_to_char(jw)) //'_'// trim(int_to_char(iw))
      call hdf_read_dataset( group_id, trim(dset_name), r_val )

      dset_name = "ep_hop_i_"// trim(int_to_char(ia)) &
         //'_'// trim(int_to_char(jw)) //'_'// trim(int_to_char(iw))
      call hdf_read_dataset( group_id, trim(dset_name), i_val )

      elph%epr(iw, jw, ia)%ep_hop = cmplx(r_val, i_val, kind=dp)
      !
      deallocate(r_val, i_val)
   enddo; enddo; enddo
   call hdf_close_group( group_id )

end subroutine read_elph_mat_wann

end module epr_hdf5_io
