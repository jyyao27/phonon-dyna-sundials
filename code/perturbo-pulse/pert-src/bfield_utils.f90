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
!
! Maintenance:
!===============================================================================
!> This module defines functions to be used in
!! the magnetotransport version of the code
module bfield_utils
   use pert_const, only: dp, pi
   use pert_data, only: tpiba,at,system_2d
   use boltz_grid, only: grid
   use pert_utils, only: fermi
   implicit none
   public :: deriv, crossprod
contains

!> Compute cross product of two vectors
function crossprod(vec1, vec2)
   implicit none
   real(dp), dimension(3), intent(in) :: vec1           !< vector 1
   real(dp), dimension(3), intent(in) :: vec2           !< vector 2
   real(dp), dimension(3) :: crossprod                  !< crossproduct(vec1,vec2)
   crossprod(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
   crossprod(2) = -vec1(1)*vec2(3) + vec1(3)*vec2(1)
   crossprod(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
end function crossprod

!> Compute the B field term in BTE
function deriv(kgrid, force, mfd, ib, ik)
   implicit none
   type(grid), intent(in) :: kgrid                      !< k-point grid
   integer, intent(in) :: ik                            !< k-point index
   integer, intent(in) :: ib                            !< band index
   real(dp), intent(in) :: mfd(3, kgrid%numb, kgrid%nk) !< vector used to compute Lorentz
   real(dp), intent(in) :: force(3)                     !< Force term - vxB for Lorentz
   real(dp) :: deriv(3)                                 !< \f$ deriv(mfd) = (force_{i}\frac{d}{dk_{i}}) mfd \f$
   real(dp) :: tmp(3)
   integer :: i
   deriv = 0.0_dp
   tmp = 0.0_dp

   !Use Marzari's central finite difference scheme
   do i = 1, kgrid%num_neighbors
     if( kgrid%neighbors(i,ik) .eq. 0 ) cycle  !no neighbors then skip this step
      tmp = kgrid%bweight(i) * dot_product( force(:), kgrid%bvec(:,i)) * &
            (mfd(:,ib,kgrid%neighbors(i,ik)) - mfd(:,ib,ik))/tpiba
      deriv = deriv + tmp 
   enddo
end function deriv

end module bfield_utils

