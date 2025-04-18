!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
! Copyright (C) 2021-2023 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Ivan Maliyov
!                         Dhruv Desai, Sergio Pineda Flores, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>; 
! Comment:
!  interface to the TDEP code.
!  mainly interact with the following code within TDEP/libolle
!  gottochblandat.f90
!  type_forceconstant_secondorder.f90
!  type_forceconstant_secondorder_io.f90
!  type_forceconstant_secondorder_loto.f90
!  type_forceconstant_secondorder_dynamicalmatrix.f90
!
! Maintenance:
!===============================================================================
module lattice_tdep
   ! modules from QE
   use kinds, only: dp
   use qe_mpi_mod, only: stdout
   use cell_base, only: alat, tpiba
   use force_constant, only: lattice_ifc, force_const
   use force_constant_thirdorder, only:lattice_ifct, force_const_third

#if defined(__TDEP)
   ! modules from TDEP
   !use konstanter, only: flyt
   use gottochblandat
   use type_crystalstructure
   use type_forceconstant_secondorder
   use type_forceconstant_thirdorder
   implicit none
   public

   !> crystal structure of the unit cell
   type(lo_crystalstructure), save :: uc_tdep
   !> second order force constant, 
   !> fc also contains born effective charges and dielectric tensor.
   type(lo_forceconstant_secondorder), save :: fc_tdep
   !> third order force constant, 
   !> fct only contains force constants
   type(lo_forceconstant_thirdorder), save :: fct_tdep
   
   public :: init_ifc_tdep
contains
   !initialize lattice_ifc with data from tdep/forceconstant fc.
   subroutine init_ifc_tdep(uc, fc, fct, ph, pht, at, cryst_tau)
      implicit none
      type(lo_crystalstructure), intent(in) :: uc
      type(lo_forceconstant_secondorder), intent(in) :: fc
      type(lo_forceconstant_thirdorder), intent(in), optional :: fct
      !
      type(lattice_ifc), intent(inout) :: ph
      type(lattice_ifct), intent(inout), optional :: pht
      !lattice vector in the unit of alat
      real(dp), intent(in) :: at(3,3)
      ! atomic position in crystal coordinate
      real(dp), intent(in) :: cryst_tau(:,:) !cryst_tau(3, ph%na)
      !
      integer :: idx, a1, a2, ir, i
      integer, allocatable :: ph2fc(:,:)
      type(force_const), pointer :: ptr
      !
      integer :: a3, irt
      integer, allocatable :: pht2fct(:,:)
      type(force_const_third), pointer :: ptrt     
      !unit of ifc should be rydberg/bohr^2, to be consistent with QE.
      real(dp) :: fc_m(3, 3), fct_m(3, 3 ,3)
      ! kelly
      real(dp) :: lambda
      !
      !the latest version of TDEP using Hatree atomic unit.
      !check "type_forceconstant_secondorder_io.f90/readfromfile"
      ! > fc%atom(a1)%pair(i)%m(j,:)=v1*lo_forceconstant_2nd_eVA_to_HartreeBohr
      ! 
      !conversion factor from Hartree/Borh^2 to Ryd/Bohr^2
      real(dp), parameter :: conv = 2.0_dp
      !debug
      real(dp) :: rvec(3,2), diff(3,2)
      integer :: i1, j1, k1

      !make sure the crystal structure in TDEP and QE are consistent.
      call check_consistency_tdep(uc, fc, ph, at, cryst_tau)
      !make sure any pairs in infile.forceconstant cann be mapped to wscell
      allocate( ph2fc(ph%max_nr, ph%na*(ph%na+1)/2) )
      call check_wscell_coverage(uc, fc, ph, ph2fc)

      idx = 0
      do a2 = 1, ph%na
      do a1 = 1, a2
         ! idx = (a2 * (a2 - 1)) / 2 + a1;  a2 >= a1
         idx = idx + 1
         ptr => ph%phi(idx)
         do ir = 1, ptr%ws_ph%nr
            fc_m(:,:) = 0.0_dp
            ! map (ir, idx) in ph to (a1, i) in fc
            i = ph2fc(ir, idx)
            !the unit of ifc in tdep is Hartree/Borh^2. Converted to rydberg/bohr^2.
            if(i > 0) fc_m(:,:) = real(fc%atom(a1)%pair(i)%m(:,:), dp) * conv
            !assign fc_m to perturbo/lattice_ifc
            ptr%ifc(:,:,ir) = cmplx(fc_m, 0.0_dp, kind=dp)
         enddo
      enddo; enddo


      !third order force constant
      !any triplets in infile.forceconstant_thirdorder should map to tripcell
      !allocate( pht2fct(pht%max_nrt, (pht%na**3+3*pht%na**2+2*pht%na)/6) )
      if (present(fct)) then
         allocate( pht2fct(pht%max_nrt, pht%na**3) )
         call check_tripcell_coverage(uc, fct, pht, pht2fct)

         idx = 0
         do a3 = 1, pht%na
         do a2 = 1, pht%na!a3
         do a1 = 1, pht%na!a2
            ! idx = ((a3-1)**3+3*(a3-1)**2+2*(a3-1))/6+(a2*(a2-1))/2+a1; 
            ! with a3 >= a2,  a2 >= a1
            idx = idx + 1
            ptrt => pht%phit(idx)
            do ir = 1, ptrt%trip_ph%nrt
               fct_m(:,:,:) = 0.0_dp
               ! map (ir, idx) in ph to (a1, i) in fc
               i = pht2fct(ir, idx)

               !the unit of ifct in tdep is Hartree/Borh^3. Converted to
               !rydberg/bohr^3.

               if(i > 0) then
                fct_m(:,:,:) = real(fct%atom(a1)%triplet(i)%m(:,:,:), dp) * conv

               !sanity check for cell distances
               rvec(:,1) = real(uc%cartesian_to_fractional(fct%atom(a1)%triplet(i)%lv2, .false., .false.), dp)
               rvec(:,2) = real(uc%cartesian_to_fractional(fct%atom(a1)%triplet(i)%lv3, .false., .false.), dp)

               diff(:,1) = pht%rvect_set(:, 1, ptrt%trip_ph%rvect(ir)) - rvec(:,1)
               diff(:,2) = pht%rvect_set(:, 2, ptrt%trip_ph%rvect(ir)) - rvec(:,2)
               if( norm2(diff(:,1))+norm2(diff(:,2)) > 2*1.0E-6_dp ) then
               call errore('init_ifc_tdep','mismatch of cell distances in TDEP',1)
               endif

               
               !sanity check for atomic distances
               rvec(:,1) = real(uc%cartesian_to_fractional(fct%atom(a1)%triplet(i)%rv2, .false., .false.), dp)
               rvec(:,2) = real(uc%cartesian_to_fractional(fct%atom(a1)%triplet(i)%rv3, .false., .false.), dp)
               
               diff(:,1) = pht%rvect_set(:, 1, ptrt%trip_ph%rvect(ir)) + cryst_tau(:,a2) - cryst_tau(:,a1)
               diff(:,1) = diff(:,1) - rvec(:,1)
               diff(:,2) = pht%rvect_set(:, 2, ptrt%trip_ph%rvect(ir)) + cryst_tau(:,a3) - cryst_tau(:,a1)
               diff(:,2) = diff(:,2) - rvec(:,2)

               if( norm2(diff(:,1))+norm2(diff(:,2)) > 2*1.0E-6_dp ) then
               call errore('init_ifc_tdep','mismatch of atomic distances in TDEP',1)
               endif

               ! sanity check for atomic indices
               if (fct%atom(a1)%triplet(i)%i2 .ne. a2 &
                   .or.fct%atom(a1)%triplet(i)%i3 .ne. a3) then 
                   call errore('init_ifc_tdep','mismatch of atomic indices in TDEP',1)
               endif

               endif
               !assign fc_m to perturbo/lattice_ifc
               ptrt%ifct(:,:,:,ir) = cmplx(fct_m, 0.0_dp, kind=dp)
            enddo
         enddo; enddo; enddo
      endif

      !now some polar parameter
      !check type_forceconstant_secondorder_io.f90/readfromfile
      !and type definition in type_forceconstant_secondorder.f90
      ph%lpol = fc%polar
      allocate( ph%pol%bcharge(3,3,ph%na) )
      !only polarcorrectiontype=3 is supported
      if(fc%polar .and. (fc%loto%correctiontype .ne. 3)) &
         call errore('init_ifc_tdep','Only support polarcorrectiontype=3 in TDEP',1)
      
      if( ph%lpol ) then
         ph%pol%epsil(:,:) = real(fc%loto%eps(:,:), dp)
         ! lambda^2 = alpha * tpiba^2
         ! kelly
         call read_tdep_lambda(lambda, fc, 'infile.forceconstant')
         ph%pol%alpha = ( lambda / tpiba )**2
         !ph%pol%alpha = ( real(fc%loto%lambda, dp) / tpiba )**2
         !
         do a1 = 1, ph%na
            !ph%pol%bcharge(:,:,a1) = real(fc%atom(a1)%Z(:,:), dp)
            ! kelly
            ph%pol%bcharge(:, :, a1) = real(fc%loto%born_effective_charges(:, :, a1), dp)
         enddo
      else
         ph%pol%epsil(:,:) = 1.0_dp
         ph%pol%alpha = 1.0_dp
         !
         ph%pol%bcharge(:,:,:) = 0.0_dp
      endif
      
      deallocate( ph2fc )
      if ( allocated(pht2fct) )  deallocate( pht2fct )
   end subroutine 


   subroutine check_consistency_tdep(uc, fc, ph, at, cryst_tau)
      implicit none
      type(lo_crystalstructure), intent(in) :: uc
      type(lo_forceconstant_secondorder), intent(in) :: fc
      !
      type(lattice_ifc), intent(inout) :: ph
      !lattice vector in the unit of alat
      real(dp), intent(in) :: at(3,3)
      !atomic position in crystal coordinate
      real(dp), intent(in) :: cryst_tau(:,:) !cryst_tau(3, ph%na)
      !
      integer :: i, a1
      real(dp) :: diff(3)
      character(len=120) :: trams, ctmp

      !print out basic structure info
      !Adapted from TDEP/type_crystalstructure.f90
      write(stdout,'(4x,a)') ' lattice vectors in TDEP calculations (in Bohr):'
      write(stdout,"(4X,A,3(2X,F18.12))") ' a1:',uc%latticevectors(:,1)
      write(stdout,"(4X,A,3(2X,F18.12))") ' a2:',uc%latticevectors(:,2)
      write(stdout,"(4X,A,3(2X,F18.12))") ' a3:',uc%latticevectors(:,3)
      write(stdout,'(4x,a)') ' ... reciprocal lattice vectors:'
      write(stdout,"(4X,A,3(2X,F18.12))") ' b1:',uc%reciprocal_latticevectors(:,1)
      write(stdout,"(4X,A,3(2X,F18.12))") ' b2:',uc%reciprocal_latticevectors(:,2)
      write(stdout,"(4X,A,3(2X,F18.12))") ' b3:',uc%reciprocal_latticevectors(:,3)
      trams=''
      do i = 1, uc%nelements
         write(ctmp,'(a,i2)') trim(uc%atomic_symbol(i))//':',uc%element_counter(i)
         trams=trim(trams)//'   '//trim(ctmp)
      enddo
      write(stdout,'(5x,a )') '... with composition   '//trim(adjustl(trams))
      write(stdout,'(5x,a/)') "----------------------------------------------"
      !
      if(fc%na.ne.ph%na) call errore("lattice_tdep","different number of atoms in TDEP and QE",1)

      !atomic position
      do a1 = 1, fc%na
         ! uc%r(:,a1): atomic position in crystal coordinate
         diff = real(uc%r(:,a1), dp) - cryst_tau(:,a1)
         if(norm2(diff) > 1.0E-6_dp) &
            call errore('lattice_tdep','different atomic position in TDEP and QE',1)
      enddo

      !lattice vector
      do i = 1, 3
         diff = real(uc%latticevectors(:,i), dp) - at(:,i)*alat
         if(norm2(diff) > 1.0E-6_dp*alat) &
            call errore('lattice_tdep','different lattice vectors in TDEP and QE',1)
      enddo
   end subroutine


   subroutine check_wscell_coverage(uc, fc, ph, map_ph2fc)
      implicit none
      type(lo_crystalstructure), intent(in) :: uc
      type(lo_forceconstant_secondorder), intent(in) :: fc
      type(lattice_ifc), intent(in) :: ph
      !map_ph2fc( ph%max_nr, na((na+1)/2) )
      integer, intent(out) :: map_ph2fc(:,:)
      !
      logical :: found
      integer :: a1, a2, i, idx, ir
      real(dp) :: rvec(3), diff(3)
      type(force_const), pointer :: ptr
      !
      real(dp), parameter :: r_tol = 1.0E-6_dp
      !init

      map_ph2fc(:,:) = 0

      do a1 = 1, fc%na
      do i = 1, fc%atom(a1)%n
         !this pair: a1, a2 + Rvec
         a2 = fc%atom(a1)%pair(i)%i2
         !Note that lv2 is in fractional in infile.forceconstant, but is converted to 
         ! Cartesian after readin via readfromfile 
         ! (see type_forceconstant_secondorder_io.f90/readfromfile)
         !here we convert it back to fractional for convenienice.
         rvec = real(uc%cartesian_to_fractional(fc%atom(a1)%pair(i)%lv2, .false., .false.), dp)
         !since (a1, a2+Rvec) and (a2, a1-Rvec) are equivalent, and only one of them
         ! are stored in Perturbo/lattice_ifc (a2 >= a1), we only check one case.
         if ( a1 > a2 ) cycle
         ! map (a1, a2) to the index used in lattice_ifc (see set_ws_cell_ph)
         ! e.g. (1, 1)->1, (1,2)->2, (2,2)->3, (1,3)->4, (2,3)->5, ...
         idx = a2 * (a2 - 1) / 2 + a1
         ptr => ph%phi(idx)
         
         found = .false.
         do ir = 1, ptr%ws_ph%nr
            !skip these already mapped (ir, idx)
            if( map_ph2fc(ir, idx) > 0 ) cycle
            diff(:) = ph%rvec_set(:, ptr%ws_ph%rvec(ir)) - rvec(:)
            ! if matches, label the matched (ir, idx), and exit the loop
            if( norm2(diff) < r_tol ) then
               found = .true.
               ! map (ir, idx) to (a1, i), and visa versa
               map_ph2fc(ir, idx) = i
               exit
            endif
         enddo
         !if not found, throw an error
         if(.not. found) call errore('check_wscell_coverage', &
            'Map tdep-ifc to perturbo failed. Try smaller rc2 in TDEP or larger q-grid in Perturbo',1)
      enddo
      enddo

   end subroutine

   subroutine check_tripcell_coverage(uc, fct, pht, map_pht2fct)
      implicit none
      type(lo_crystalstructure), intent(in) :: uc
      type(lo_forceconstant_thirdorder), intent(in) :: fct
      type(lattice_ifct), intent(in) :: pht
      !map_ph2fc( pht%max_nrt, (na**3+3*na**2+2*na)/6 )
      integer, intent(out) :: map_pht2fct(:,:)
      !
      logical :: found
      integer :: a1, a2, a3, t, idx, ir
      real(dp) :: rvec(3,2), diff(3,2)
      type(force_const_third), pointer :: ptrt
      !
      real(dp), parameter :: r_tol = 1.0E-6_dp
      
      !init
      map_pht2fct(:,:) = 0

      do a1 = 1, fct%na
      do t = 1, fct%atom(a1)%n
         !this triplet: a1, a2, a3 + Rvec
         a2 = fct%atom(a1)%triplet(t)%i2
         a3 = fct%atom(a1)%triplet(t)%i3
         !Note that lv2, lv3 is in fractional in infile.forceconstant_thirdorder
         !but is converted to Cartesian after readin via readfromfile 
         ! (see type_forceconstant_thirdorder_io.f90/readfromfile)
         !here we convert it back to fractional for convenienice.
         rvec(:,1) = real(uc%cartesian_to_fractional(fct%atom(a1)%triplet(t)%lv2, .false., .false.), dp)
         rvec(:,2) = real(uc%cartesian_to_fractional(fct%atom(a1)%triplet(t)%lv3, .false., .false.), dp)
         !since (a1, a2+r1, a3+r2) and (a2, a1-r1, a3+r2-r1) are equivalent, and only one of
         !them
         ! are stored in Perturbo/lattice_ifc (a3 >= a2 >= a1), we only check one case.
!!========debug========
!         if ( a1 > a2 .or. a2 > a3 ) cycle
!         ! map (a1, a2) to the index used in lattice_ifc (see set_ws_cell_ph)
!         ! e.g. (1, 1, 1)->1, (1, 1, 2)->2, (1, 2, 2)->3, (2, 2, 2)->4, (1, 1, 3)->5, ...
!
!!  debug -symmetry
!!         idx = ((a3-1)**3+3*(a3-1)**2+2*(a3-1))/6
!!         idx = idx + a2 * (a2 - 1) / 2 + a1
!!========debug========
         idx = (a3-1)*(fct%na**2)+(a2-1)*fct%na+a1
         ptrt => pht%phit(idx)

         found = .false.
         do ir = 1, ptrt%trip_ph%nrt
            !skip these already mapped (ir, idx)
            if( map_pht2fct(ir, idx) > 0 ) cycle
            diff(:,1) = pht%rvect_set(:, 1, ptrt%trip_ph%rvect(ir)) - rvec(:,1)
            diff(:,2) = pht%rvect_set(:, 2, ptrt%trip_ph%rvect(ir)) - rvec(:,2)
            ! if matches, label the matched (ir, idx), and exit the loop
            if( norm2(diff(:,1))+norm2(diff(:,2)) < 2*r_tol ) then
               found = .true.
               ! map (ir, idx) to (a1, i), and visa versa
               map_pht2fct(ir, idx) = t
               exit
            endif
         enddo
         !if not found, throw an error
         if(.not. found) call errore('check_tripcell_coverage', &
            'Map tdep-ifc to perturbo failed. Try smaller rc3 in TDEP or larger q-grid in Perturbo',1)
      enddo
      enddo

   end subroutine

   subroutine read_tdep_lambda(lambda, fc, fn)
      use gottochblandat, only: open_file
      implicit none

      real(dp), intent(out) :: lambda
      class(lo_forceconstant_secondorder), intent(in) :: fc
      character(len=*), intent(in) :: fn
      integer :: u, a1, i, j

      u = open_file('in', trim(fn)) 
      read (u, *)
      read (u, *)
      do a1 = 1, fc%na
         read (u, *)
         do i = 1, fc%atom(a1)%n
            read (u, *)
            read (u, *)
            do j = 1, 3
               read(u, *)
            enddo
         enddo
      enddo   
      read (u, *) i
      if (i .eq. 0) then
         call errore('lattice_tdep',"material is not polar",1)
      else
         do j = 1, 3
            read (u, *)
         enddo
         read (u, *) lambda
      endif
      close(u)
   end subroutine

#endif

end module lattice_tdep
