!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  setup the Wigner Seitz supercells for wannier interpolation.
!
! Maintenance:
!===============================================================================

module triplet_cell
   use kinds, only: dp
   use wigner_seitz_cell, only: get_length, set_cutoff_large, set_cutoff_small
   implicit none
   public

   type :: trip_cell
      integer :: nrt
      integer :: reserved_ = 0 !not used, for data alignment.
      integer, pointer :: ndegt(:) => null() !ndegt(nrt)
      integer, pointer :: rvect(:) => null() !rvect(nrt)
   end type

   type :: trip_vector_set
      integer :: nvect
      integer :: reserved_ = 0 !not used, for data alignment.
      !syp - this is the most different point with e-ph
      !    - there are three atoms involved in the unharmoic term
      !    - assume the first aton in the home unit cell
      !    - so there are two other atoms in other unit cells in the Bvk supercell
      !    - vect(:,1,nvect) store one of the two atoms
      !    - vect(:,2,nvect) store the other atom
      integer, pointer :: vect(:,:,:) => null()  !vect(3,2,nvect)
      integer, pointer :: idx(:) => null()  ! idx(nr1*nr2*nr3,2)
      integer, pointer :: nim(:) => null()  ! nim(nr1*nr2*nr3,2)
      !NOTE: idx and nim are used to access image vectors of (i, j, k) pairs
      ! with ir = (k-1)*nr1*nr2 + (j-1)*nr1 + (i-1) + 1
      !
      ! the images of R=(i,j,k) are stored in vec: 
      !    from index idx(ir,1/2) to [idx(ir,1/2)+nim(ir,1/2)-1]
      ! where nim(ir,1/2) stores the number of images of R
      ! and nvec = sum(nim)
   end type

   integer, parameter, private :: trip_search_range = 2
   real(dp), parameter, private :: eps6 = 1.0E-6_dp

   public :: init_rvect_images, reset_rvect_images, set_trip_cell
contains

subroutine init_rvect_images(nr1, nr2, nr3, at, rvect_images, large_cutoff)
   implicit none
   integer, intent(in) :: nr1, nr2, nr3
   real(dp), intent(in) :: at(3,3)
   type(trip_vector_set), intent(inout) :: rvect_images
   logical, intent(in), optional :: large_cutoff
   !local
   logical :: lcut
   real(dp) :: cutoff
   integer :: ws_dim(3), vec(3), rvec(3), nim, i, j, k, mi, mj, mk, tot, ir
   integer :: vec2(3), rvec2(3), nim2, i2, j2, k2, mi2, mj2, mk2, tot2, ir2

   ws_dim(1) = nr1
   ws_dim(2) = nr2
   ws_dim(3) = nr3

   lcut = .false.
   if(present(large_cutoff)) lcut = large_cutoff
   if( lcut ) then
      cutoff = set_cutoff_large(ws_dim, at)
   else
      cutoff = set_cutoff_small(ws_dim, at)
   endif

   if( associated(rvect_images%idx) ) deallocate( rvect_images%idx )
   if( associated(rvect_images%nim) ) deallocate( rvect_images%nim )
   allocate( rvect_images%idx( (nr1*nr2*nr3)**2), rvect_images%nim( (nr1*nr2*nr3)**2) )
   
   tot = 0
   ir  = 0
   do k = 1, nr3; do j = 1, nr2; do i = 1, nr1
      vec(1:3) = (/i-1, j-1, k-1/)
      do k2 = 1, nr3; do j2 = 1, nr2; do i2 = 1, nr1
            vec2(1:3) = (/i2-1, j2-1, k2-1/)

            ir = ir + 1
            nim = 0

            do mk = -trip_search_range, trip_search_range
            do mj = -trip_search_range, trip_search_range
            do mi = -trip_search_range, trip_search_range
               rvec(1:3) = (/mi, mj, mk/) * ws_dim(1:3) + vec(1:3)

              do mk2 = -trip_search_range, trip_search_range
              do mj2 = -trip_search_range, trip_search_range
              do mi2 = -trip_search_range, trip_search_range
                 rvec2(1:3) = (/mi2, mj2, mk2/) * ws_dim(1:3) + vec2(1:3)
         !
                 if( get_length(real(rvec, kind=dp), at) < cutoff .and. &
                 get_length(real(rvec2, kind=dp), at) < cutoff) nim = nim + 1

              enddo; enddo; enddo
            enddo; enddo; enddo
            if(nim < 1) call errore('init_trip_cell','nim < 1',1)

            ! init idx, nim
            rvect_images%idx(ir) = tot + 1
            rvect_images%nim(ir) = nim
            tot = tot + nim
      enddo; enddo; enddo
   enddo; enddo; enddo

   rvect_images%nvect = tot
   if( associated(rvect_images%vect) ) deallocate( rvect_images%vect )
   allocate( rvect_images%vect(3, 2, tot) )

   ! collect all the vectors 
   nim = 0
   do k = 1, nr3; do j = 1, nr2; do i = 1, nr1
      vec(1:3) = (/i-1, j-1, k-1/)
      do k2 = 1, nr3; do j2 = 1, nr2; do i2 = 1, nr1
         vec2(1:3) = (/i2-1, j2-1, k2-1/)

         do mk = -trip_search_range, trip_search_range
         do mj = -trip_search_range, trip_search_range
         do mi = -trip_search_range, trip_search_range
            rvec(1:3) = (/mi, mj, mk/) * ws_dim(1:3) + vec(1:3)

            do mk2 = -trip_search_range, trip_search_range
            do mj2 = -trip_search_range, trip_search_range
            do mi2 = -trip_search_range, trip_search_range
               rvec2(1:3) = (/mi2, mj2, mk2/) * ws_dim(1:3) + vec2(1:3)

               if( get_length(real(rvec, kind=dp), at) < cutoff .and. &
                   get_length(real(rvec2, kind=dp), at) < cutoff) then
                   nim = nim + 1
                   rvect_images%vect(:,1,nim) = rvec
                   rvect_images%vect(:,2,nim) = rvec2
               endif
            enddo; enddo; enddo
         enddo; enddo; enddo
      enddo; enddo; enddo
   enddo; enddo; enddo
end subroutine init_rvect_images


! <0,a| H_[c] | R,b>
subroutine set_trip_cell &
      (matom, nr1, nr2, nr3, at, rvect_images, trip, tau_a, tau_b, tau_c, ref_c)
   implicit none
   integer,  intent(in) :: matom, nr1, nr2, nr3
   real(dp), intent(in) :: at(3,3)
   type(trip_vector_set), intent(in) :: rvect_images
   type(trip_cell), intent(inout) :: trip
   ! tau_a and tau_b are within unit cell and in crystal coordinate.
   real(dp), intent(in) :: tau_a(3), tau_b(3), tau_c(3)
   real(dp), intent(in), optional :: ref_c(3) ! in crystal coordinate
   !local
   integer :: nsize, idx, nim, n, m, tot, ir, tot_r
   real(dp) :: vec_b(3), vec_b2(3)
   integer, allocatable :: itmp(:), ndegt(:)
   real(dp), allocatable :: dist(:), dist2(:)
   
   logical :: output_rvec
   integer :: uout
   integer :: irim, rv1(3), rv2(3)
   character(len=80) :: fname
   ! output rvec
   output_rvec = .false.

   nsize = maxval( rvect_images%nim(:) )
   tot_r = (nr1 * nr2 * nr3)**2
   !work array
   allocate( itmp( rvect_images%nvect ), ndegt(tot_r) )
   itmp = 0;  ndegt = 0

!$omp parallel default(shared) private(ir, idx, nim, n, vec_b, vec_b2, dist, dist2)
   allocate( dist(nsize), dist2(nsize) )
!$omp do schedule(guided)
   do ir = 1, tot_r
      ! R = (i, j, k) => 
      !    ir = (k-1)*nr1*nr2 + (j-1)*nr1 + (i-1) + 1
      !
      idx = rvect_images%idx(ir)
      nim = rvect_images%nim(ir)
      do n = 1, nim
         vec_b = real(rvect_images%vect(:, 1, n+idx-1), kind=dp) + tau_b
         vec_b2 = real(rvect_images%vect(:, 2, n+idx-1), kind=dp) + tau_c
         dist(n) = get_length(vec_b-tau_a, at)
         dist2(n) = get_length(vec_b2-tau_a, at)
         ! 
         !if(present(ref_c)) dist(n) = dist(n) + get_length(vec_b-ref_c, at)
      enddo
      ! find the shortest one(s) or its equivalences
      dist(1:nim) = dist(1:nim) - minval(dist(1:nim))
      dist2(1:nim) = dist2(1:nim) - minval(dist2(1:nim))
      ndegt(ir) = count( dist(1:nim)+dist2(1:nim) < eps6 )
      
      do n = 1, nim
         if( dist(n)+dist2(n) < eps6 ) itmp( n+idx-1 ) = ndegt(ir)
      enddo
   enddo
!$omp end do
   deallocate( dist, dist2 )
!$omp end parallel
   !sanity check
   if( any( ndegt < 1 ) ) call errore('set_wigner_seitz_cell','ndegt < 1',1)
   
   tot = sum( ndegt(:) )
   trip%nrt = tot
   if( associated( trip%ndegt ) ) deallocate( trip%ndegt )
   if( associated( trip%rvect ) ) deallocate( trip%rvect )
   allocate( trip%rvect(tot), trip%ndegt(tot) )

   !collect vectors into trip%rvect, only store their index in rvect_images%vec
   m = 0
   do n = 1, rvect_images%nvect
      if( itmp(n) > 0 ) then
         m = m + 1
         trip%rvect(m) = n
         trip%ndegt(m) = itmp(n)
      endif
   enddo
   
   deallocate(itmp, ndegt)

   if (output_rvec) then
      uout = 999 ! hard coded
      fname = 'ifc_index.out'
      open(unit=uout, file=trim(fname), status='unknown', form='formatted')
      write(uout, 1030)
      write(uout, 1031)
      do ir = 1, trip%nrt
         ! irim index the image in rvect_images%vect
         irim = trip%rvect(ir)
         rv1 = rvect_images%vect(:,1,irim)
         rv2 = rvect_images%vect(:,2,irim)
         write(uout, 1032) matom, ir, irim, rv1, rv2
      enddo
   endif

1030 format('# IFC index mapping #')
1031 format('# atomic index        ir          iR          rvec1                rvec2 #')
1032 format(i10, 2x, i10, 2x, i10, 2x, 3(i10,",",2x), 3(i10, ",", 2x))
end subroutine set_trip_cell


subroutine reset_rvect_images(rvect_images)
   implicit none
   type(trip_vector_set), intent(inout) :: rvect_images
   
   rvect_images%nvect = 0
   if( associated(rvect_images%vect) ) deallocate( rvect_images%vect )
   if( associated(rvect_images%idx) ) deallocate( rvect_images%idx )
   if( associated(rvect_images%nim) ) deallocate( rvect_images%nim )
end subroutine reset_rvect_images

end module triplet_cell
