module force_constant_thirdorder
   use kinds, only: dp
   use polar_correction, only: polar_parameter
   use triplet_cell, only: trip_cell, trip_vector_set, init_rvect_images, &
      reset_rvect_images, set_trip_cell
   implicit none
   private

   type, public :: lattice_ifct
      integer :: na  ! number of atoms per unit cell
      integer :: nm  ! number of modes = 3*na
      ! phit() -xt
      type(force_const_third), pointer :: phit(:) => null()
      !
      integer :: max_nrt !  max value of phit(:)%wst_ph%nrt
      !
      integer :: nrvect
      real(dp), pointer :: rvect_set(:,:,:) => null() ! rvec_set(3,2,nrvect)

   end type

   !-xt
   type, public :: force_const_third
      type(trip_cell), pointer :: trip_ph => null()
      !real(dp), pointer :: ifct(:,:,:,:) => null() ! ifc(3, 3, 3, ws_ph%nr)
      complex(dp), pointer :: ifct(:,:,:,:) => null() ! ifc(3, 3, 3, ws_ph%nr)
   end type
   
   type, public :: expikr_set
      type(expikr), pointer :: exps_storage(:,:) => null() ! expstorage(nq,pht%na^3) 
   end type

   type, public :: expikr
      complex(dp), pointer :: exps(:) => null() ! exps(tc%nrt)
   end type
   
   public :: set_tripcell_ph
contains

subroutine set_tripcell_ph(pht, qdim, nat, at, cryst_tau)
   implicit none
   type(lattice_ifct), intent(inout) :: pht
   integer, intent(in) :: qdim(3), nat
   ! atomic position in crystal coordinate
   real(dp), intent(in) :: at(3,3), cryst_tau(3,nat) 
   !
   integer :: nelem, m, ia, ja, ka, nvect, max_nvect
   type(trip_vector_set) :: rvect_images

   ! make sure ph is an empty object
   if(associated(pht%phit) .or. associated(pht%rvect_set)) &
      call errore('set_ws_tripcell_ph','ws_cell object is already initialized',1)
   
   pht%na = nat
   pht%nm = 3*nat
!   nelem = (nat**3+3*nat**2+2*nat)/6
   nelem = nat**3
   !setup all the the vector set for wigner seitz cell
   allocate( pht%phit(nelem) )

   call init_rvect_images( qdim(1), qdim(2), qdim(3), at, rvect_images )

   max_nvect = 0
   m = 0
! debug - symmetry turned off
!   do ka = 1, nat
!   do ja = 1, ka
!   do ia = 1, ja
      !m = ((ka-1)**3+3*(ka-1)**2+2*(ka-1))/6+(ja*(ja-1))/2+ia

    do ka = 1, nat
    do ja = 1, nat
    do ia = 1, nat
! debug
      !m = (ka-1)*(nat**2)+(ja-1)*nat+ia
      m = m + 1
      allocate( pht%phit(m)%trip_ph )
      !
      call set_trip_cell(m, qdim(1), qdim(2), qdim(3), at, rvect_images, &
      pht%phit(m)%trip_ph, cryst_tau(:,ia), cryst_tau(:,ja), cryst_tau(:,ka))

      !allocate space and init ham_r
      nvect = pht%phit(m)%trip_ph%nrt
      allocate( pht%phit(m)%ifct(3, 3, 3, nvect) )
      !pht%phit(m)%ifct(:,:,:) = cmplx(0.0_dp, 0.0_dp, 0.0_dp, kind=dp)

      if( max_nvect < nvect ) max_nvect = nvect
   enddo; enddo; enddo
   !
   pht%max_nrt = max_nvect
   
   call setup_rvect_set_ph(pht, rvect_images)
   call reset_rvect_images( rvect_images )
end subroutine set_tripcell_ph


subroutine setup_rvect_set_ph(pht, r_images)
   implicit none
   type(lattice_ifct), intent(inout) :: pht
   type(trip_vector_set), intent(in) :: r_images
   !
   integer :: m, ir, irp, nelem
   integer, allocatable :: rvec_label(:)
   type(force_const_third), pointer :: ptrt
   
   !nelem = (pht%na**3+3*pht%na**2+2*pht%na)/6
   nelem = pht%na**3

   !set up rvec_set
   allocate( rvec_label( r_images%nvect ) )

   rvec_label = 0
   do m = 1, nelem
      ptrt => pht%phit(m)
      
      do ir = 1, ptrt%trip_ph%nrt
         irp = ptrt%trip_ph%rvect(ir)
         rvec_label( irp ) = rvec_label( irp ) + 1
      enddo
   enddo
   
   pht%nrvect = count(rvec_label > 0)
   allocate( pht%rvect_set(3, 2, pht%nrvect) )
   pht%rvect_set = 0.0E0_dp
   
   m = 0
   do irp = 1, r_images%nvect
      !this R is used in the wigner seitz cell
      if( rvec_label(irp) > 0 ) then
         m = m + 1
         pht%rvect_set(1:3, 1:2, m) = real(r_images%vect(1:3, 1:2, irp), dp)
         rvec_label(irp) = m
      endif
   enddo
   
   !update trip_ph%rvect
   do m= 1, nelem
      ptrt => pht%phit(m)

      do ir = 1, ptrt%trip_ph%nrt
         irp = ptrt%trip_ph%rvect(ir)
         ptrt%trip_ph%rvect(ir) = rvec_label(irp)
      enddo
   enddo

   deallocate( rvec_label )
end subroutine setup_rvect_set_ph

end module force_constant_thirdorder
