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
!   basic functions used in boltzmann subroutines
!
! Maintenance:
!===============================================================================

module boltz_utils
   use pert_const, only: dp
   private
   public :: num2kpt
   public :: kpt2num
   public :: idx2band
   public :: band2idx
   public :: kpts_plus
   public :: kpts_minus
   public :: inside_win
   public :: velprod
   public :: kpts_plus_invert
   public :: calc_adaptive_smear
   public :: find_q_neighbor
   public :: fold_k
   public :: num2ipt
contains

pure function num2ipt(num, nk) result(ipt)
   implicit none
   integer, intent(in) :: num   !< non-negative integer used to determine `kpt`
   integer, intent(in) :: nk(3) !< number of k points in each direction 
   integer :: ipt(3)            !< used to store [i, j, k]
   ! k = mod(num, nk3)
   ipt(3) = mod(num, nk(3))
   ! j = mod(num / nk3, nk2)
   ipt(2) = mod(num / nk(3), nk(2))
   ! i = num / (nk2*nk3)
   ipt(1) = num / (nk(2) * nk(3))
end function num2ipt

!> Take a non-negative integer; `num`, and the number of k points
!! in 3D; `nk`= (nk1, nk2, nk3). `num2kpt` finds the corresponding
!! k point; `kpt`, based on splitting up `num` as follows
!! `num` = k + j * nk3 + i * nk2 * nk3
!! where \f$ i \in [0, nk1 - 1] \f$, \f$ j \in [0, nk2 - 1] \f$
!! and \f$ k \in [0, nk3 - 1] \f$. 
!!
!! The k point; `kpt`, is given as a fraction [i/nk1, j/nk2, k/nk3]
!! and folded to the Gamma-centered FBZ in crystal coordinates such
!! that each `kpt` value is between -0.5 and 0.5
pure function num2kpt(num, nk) result(kpt)
   implicit none
   integer, intent(in) :: num   !< non-negative integer used to determine `kpt`
   integer, intent(in) :: nk(3) !< number of k points in each direction 
   integer :: ipt(3)            !< used to store [i, j, k]
   integer :: i                 !< used in loop
   real(dp) :: kpt(3)           !< k point in Gamma-centered FBZ crystal coordinates
   ! k = mod(num, nk3)
   ipt(3) = mod(num, nk(3))
   ! j = mod(num / nk3, nk2)
   ipt(2) = mod(num / nk(3), nk(2))
   ! i = num / (nk2*nk3)
   ipt(1) = num / (nk(2) * nk(3))

   ! kpt as a fraction of total k points
   kpt(:) = real(ipt(:), dp) / real(nk(:), dp)
   ! fold to the Gamma-centered FBZ, in crystal coordinate
   ! if kpt(i) rounds to 1, then subtract nk(i) so it is
   ! between -nk(i)/2 and 0
   do i = 1, 3
      ipt(i) = merge(ipt(i)-nk(i), ipt(i), nint(kpt(i)) > 0)
   enddo
   ! get the folded kpt as a fraction
   kpt(:) = real(ipt(:), dp) / real(nk(:), dp)
end function num2kpt

! map (kx, ky, kz) to num, if num < 0, then (kx, ky, kz) is not in the list
pure function kpt2num(kpt, nk) result(num)
   implicit none
   real(dp), intent(in) :: kpt(3)
   integer, intent(in) :: nk(3)
   integer :: num
   real(dp), parameter :: eps = 1.0E-5_dp
   ! local variables
   integer :: r(3)  !, i 
   real(dp) :: xkr(3), dis(3)
   ! init num to the default value 
   num = -1
   ! fold to the Gamma-centered FBZ, in crystal coordinate
   ! and check if kpt(i)*nk(i) is a integer or not.
   xkr(:) = (kpt(:) - nint(kpt(:))) * nk(:)
   dis(:) =  xkr(:) - nint(xkr(:))
   ! return -1 if (kx, ky, kz) is not in the k-mesh.
   if( sqrt(dot_product(dis, dis)) > eps ) return
   ! ri = 0...nki-1; r(1)->i; r(2)->j; r(3)->k
   r(:) = mod( nint(xkr(:)+2*nk(:)), nk(:) )
   ! num = k + j*nk3 + i*nk2*nk3 + 1
   num = r(3) + r(2)*nk(3) + r(1)*nk(2)*nk(3)
end function kpt2num

! given index of ikq and ik, compute the index of xq = ikq-ik on q-grid.
pure function kpts_minus(ikq, ik, nk, nq) result(iq)
   implicit none
   integer, intent(in) :: ikq, ik, nk(3), nq(3)
   integer :: iq, q(3) !, i
   real(dp) :: xq(3)
   ! compute difference in i, j, k
   q(1) =     ikq/(nk(2)*nk(3)) -     ik/(nk(2)*nk(3))
   q(2) = mod(ikq/nk(3), nk(2)) - mod(ik/nk(3), nk(2))
   q(3) = mod(ikq, nk(3))       - mod(ik, nk(3))
   !
   !coordinates of q-points
   xq(:) = real(q(:), dp) / real(nk(:), dp)
   !get the index, if xq is not on the q-grid, then iq = -1
   iq = kpt2num(xq, nq)
end function kpts_minus

! given index of ikq and ik, compute the index of xq = ikq+ik on q-grid.
pure function kpts_plus(ikq, ik, nk, nq) result(iq)
   implicit none
   integer, intent(in) :: ikq, ik, nk(3), nq(3)
   integer :: iq, q(3) !, i
   real(dp) :: xq(3)
   ! compute difference in i, j, k
   q(1) =     ikq/(nk(2)*nk(3)) +     ik/(nk(2)*nk(3))
   q(2) = mod(ikq/nk(3), nk(2)) + mod(ik/nk(3), nk(2))
   q(3) = mod(ikq, nk(3))       + mod(ik, nk(3))
   ! 
   !coordinates of q-points
   xq(:) = real(q(:), dp) / real(nk(:), dp)
   !get the index, if xq is not on the q-grid, then iq = -1
   iq = kpt2num(xq, nq)
end function kpts_plus

! given index of ikq and ik, compute the index of xq = -(ikq+ik) on q-grid.
pure function kpts_plus_invert(ikq, ik, nq) result(iq)
   implicit none
   !integer, intent(in) :: ikq, ik, nk(3), nq(3)
   !integer :: iq, q(3) !, i
   !real(dp) :: xq(3)
   integer, intent(in) :: ikq, ik, nq(3) ! ikq, ik starts with zero
   integer :: ipt0(3), ipt1(3)
   integer :: iq
   !! compute difference in i, j, k
   !q(1) =     ikq/(nk(2)*nk(3)) +     ik/(nk(2)*nk(3))
   !q(2) = mod(ikq/nk(3), nk(2)) + mod(ik/nk(3), nk(2))
   !q(3) = mod(ikq, nk(3))       + mod(ik, nk(3))
   ! 
   !!coordinates of q-points
   !xq(:) = - real(q(:), dp) / real(nk(:), dp)
   !!get the index, if xq is not on the q-grid, then iq = -1
   !iq = kpt2num(xq, nq)
   
   ipt0 = num2ipt(ikq, nq)
   ipt1 = num2ipt(ik, nq)
   iq = dot_product(mod(2*nq-ipt0-ipt1, nq), (/nq(3)*nq(2), nq(3), 1/))
end function kpts_plus_invert

! map index to (mkq, nk, mu)
pure function idx2band(idx, nbnd) result(band)
   implicit none
   integer, intent(in) :: idx
   integer, intent(in) :: nbnd
   integer :: band(3)
   ! mu <= nmodes, nk <= nbnd, mkq <= nbnd
   ! index = (mu-1)*nbnd*nbnd + (nk-1)*nbnd + (mkq-1)
   band(1) = mod( idx, nbnd ) + 1
   band(2) = mod( idx/nbnd, nbnd ) + 1
   band(3) = idx / (nbnd*nbnd) + 1
end function idx2band

! map (mkq, nk, mu) to index
pure function band2idx(band, nbnd) result(idx)
   implicit none
   integer, intent(in) :: band(3) ! mkq, nk, mu
   integer, intent(in) :: nbnd ! numb, numb, nmodes
   integer :: idx
   ! index = (mu-1)*nbnd*nbnd + (nk-1)*nbnd + (mkq-1)
   idx = ( band(3)-1 )*nbnd*nbnd + ( band(2)-1 )*nbnd + ( band(1)-1 )
end function band2idx

! check if xkg has energy levels inside [emin, emax]
! only bands in [bmin, bmax] are considered.
function inside_win(el, xkg, emin, emax, bmin, bmax)
   use band_structure, only: electron_wann, solve_eigenvalue_vector
   implicit none
   type(electron_wann), intent(in) :: el
   integer,  intent(in) :: bmin, bmax
   real(dp), intent(in) :: emin, emax, xkg(3)
   integer :: inside_win
   ! local variables
   real(dp) ::  eval(el%nb)
   integer :: lower, upper, ib 

   call solve_eigenvalue_vector(el, xkg, eval)
   ! only bands in [bmin, bmax] are considered.
   lower = bmin - 1
   upper = bmax + 1
   do ib = bmin, bmax
      if( eval(ib) < emin ) lower = lower + 1
      if( eval(ib) > emax ) upper = upper - 1
   enddo
   ! no level inside [emin, emax]
   ! e.g. all levels are lower than emin .or larger than emax
   ! or have levels lower than emin and levels larger than emax.
   if( upper .eq. (lower+1) ) then
      inside_win = lower
   ! have levels inside [emin, emax]
   elseif(upper > (lower+1) ) then
      inside_win = -2
   else
      ! error appears, not a logical result
      inside_win = -4
   endif
   return
end function inside_win

pure function velprod(v1, v2) result(vv)
   implicit none
   real(dp), intent(in) :: v1(3)
   real(dp), intent(in) :: v2(3)
   !
   real(dp) :: vv(6)

   vv(1) = v1(1) * v2(1)  ! XX: v_x * v_x
   vv(2) = v1(1) * v2(2)  ! XY: v_x * v_y
   vv(3) = v1(2) * v2(2)  ! YY: v_y * v_y
   vv(4) = v1(1) * v2(3)  ! XZ: v_x * v_z
   vv(5) = v1(2) * v2(3)  ! YZ: v_y * v_z
   vv(6) = v1(3) * v2(3)  ! ZZ: v_z * v_z
end function velprod

! computes adaptive smearing for phonon grids 
!(or electron-phonon) given group velocities and lattice structure 
! max{ |(vk - vq) * G1/N1, G2/N2, G3/N3| }
! using v in cartesian coordinates
pure function calc_adaptive_smear(vk,vkq,ndim) result(smear)
   use pert_data, only: nat,bg
   implicit none
   ! group velocities at vk at k and vkq at k+q, in cartesian coords 
   real(dp), intent(in) :: vk(3), vkq(3)
   ! largest dimension
   integer, intent(in) :: ndim
   ! adaptive smearing value
   real(dp) :: smear
!   do idir = 1, nat
!      del_e(idir) = abs(dot_product( (vk(:) - vkq(:)), cryst_tau(:,idir)) ) / ndim(idir)
!   enddo
   !smear = (del_e(1)+del_e(2))/2.0_dp
   !smear = maxval(del_e)
   
   ! just use dv  
   smear = norm2(vk(:) - vkq(:)) * norm2(bg(1,:))/ ndim
end function calc_adaptive_smear

subroutine random_seed_gen()
   integer :: i, n, clock
   integer, allocatable :: seed(:)

   call random_seed(size=n)
   allocate(seed(n))

   call system_clock(count=clock)
   seed = clock + 37 * (/(i-1, i=1, n)/) + my_pool_id
   call random_seed(put=seed)
   
   deallocate(seed)
end subroutine 

subroutine random_choice(list_len, ind)
   implicit none
   integer, intent(in) :: list_len
   integer, intent(out) :: ind 
   real(dp) :: rand
   call random_seed_gen()
   
   call random_number(rand)
   ! determine index based on random number
   ! interval is 
   ind = int(rand*list_len)+1

end subroutine 

pure function fold_k(k) result(kf)
   implicit none
   real(dp), intent(in) :: k
   real(dp) :: kf
   ! k is from -1 to 1
   ! if k is larger than or equal to 0.5, k - 1.0
   ! if k is less than -0.5, k + 1.0
   ! k'=k+0.5 -> if k'>=1, int(k')=1, k-floor(k')
   ! if k'<0, floor(k')=-1,
   kf = k - floor(k+0.5_dp)
   !kf = merge(k-1.0_dp, k, k > 0.5_dp - 1.0E-8_dp)
   !kf = merge(k+1.0_dp, k, k < -0.5_dp - 1.0E-8_dp)

end function fold_k

!subroutine find_q_neighbor(qpt, coarse_dim, fine_weights, iq)
subroutine find_q_neighbor(qpt, coarse_dim, fine_weights, iqlist_out)
   implicit none
   real(dp), intent(in) :: qpt(3)
   integer, intent(in) :: coarse_dim(3)
   real(dp), intent(in) :: fine_weights(3)
   !integer, intent(out) :: iq
   integer :: iq

   real(dp), allocatable :: q(:,:)
   integer :: on_grid(3)
   integer :: i, j, i2, i3, num_nb, ind, temp_index, temp
   real(dp), allocatable :: rnd(:)
   integer, allocatable :: indices(:), iqlist(:)
   integer, allocatable, intent(out) :: iqlist_out(:)
   !integer, allocatable :: iqlist(:)
   real(dp) :: modqp, modqm
   ! if iq_fine falls on qg
   on_grid(:) = 1
   do i = 1,3
      if (mod(abs(nint(qpt(i)/fine_weights(i))),2) == 1) then 
         on_grid(i) = 2
      endif
      !write(*,*) mod(nint(qpt(i)/fine_weights(i)),2), on_grid(i)
   enddo
   num_nb = on_grid(1)*on_grid(2)*on_grid(3)
   if (num_nb == 1) then
      iq = kpt2num(qpt, coarse_dim) + 1
      allocate(iqlist(1))
      iqlist(1) = iq
   else
      allocate(q(num_nb, 3))  
      allocate(iqlist(num_nb))  
      iqlist(:) = 1
      if (num_nb == 8) then
         ! set x negative at 1234
         do j = 1, 4
            modqm = fold_k(qpt(1) - fine_weights(1))
            modqp = fold_k(qpt(1) + fine_weights(1))
            q(j,1) = modqm 
            q(j+4,1) = modqp 
         enddo
         ! set y negative at 1357
         do j = 1, 4
            modqm = fold_k(qpt(2) - fine_weights(2))
            modqp = fold_k(qpt(2) + fine_weights(2))
            q(2*j-1,2) = modqm 
            q(2*j,2) = modqp
         enddo
         ! set z (negative at 1256)
         do j = 1, 2
            modqm = fold_k(qpt(3) - fine_weights(3))
            modqp = fold_k(qpt(3) + fine_weights(3))
            q(j,3) = modqm 
            q(j+4,3) = modqm 
            q(j+2,3) = modqp 
            q(j+6,3) = modqp
         enddo
      else if (num_nb == 2) then
         do i = 1, 3
            if (on_grid(i) > 1) then ! not ongrid
               modqm = fold_k(qpt(i) - fine_weights(i))
               modqp = fold_k(qpt(i) + fine_weights(i))
               q(1,i) = modqm 
               q(2,i) = modqp 
               i2 = mod(i,3)+1
               i3 = mod(i+1,3)+1
               q(:,i2) = qpt(i2) 
               q(:,i3) = qpt(i3)
            endif
         enddo
      else ! (num_nb == 4)
         do i = 1, 3
            if (on_grid(i) < 2) then ! ongrid
               q(:,i) = qpt(i)
               i2 = mod(i,3)+1
               i3 = mod(i+1,3)+1
               modqm = fold_k(qpt(i2) - fine_weights(i2))
               modqp = fold_k(qpt(i2) + fine_weights(i2))
               q(1,i2) = modqm 
               q(2,i2) = modqm
               q(3,i2) = modqp 
               q(4,i2) = modqp
               modqm = fold_k(qpt(i3) - fine_weights(i3))
               modqp = fold_k(qpt(i3) + fine_weights(i3))
               q(1,i3) = modqm 
               q(2,i3) = modqp 
               q(3,i3) = modqm 
               q(4,i3) = modqp
            endif
         enddo
      endif
      do j = 1, num_nb 
         iqlist(j) = kpt2num(q(j,:), coarse_dim) + 1
      enddo
      !allocate(rnd(num_nb))
      !allocate(indices(num_nb))
      !do i = 1, num_nb
      !   call random_number(rnd(i))
      !   indices(i) = i
      !enddo
      !do i = 1, num_nb-1
      !   do j = 1, num_nb-i
      !      if (rnd(j) > rnd(j+1)) then
      !         ! Swap numbers
      !         temp = rnd(j)
      !         rnd(j) = rnd(j+1)
      !         rnd(j+1) = temp
      !         ! Swap indices
      !         temp_index = indices(j)
      !         indices(j) = indices(j+1)
      !         indices(j+1) = temp_index
      !      endif
      !   enddo
      !enddo
      !write(*,*) indices
      !do i = 1, num_nb
      !   iqlist_out(i) = iqlist(indices(i))
      !enddo
      !deallocate(rnd, indices)
   endif
   allocate(iqlist_out(num_nb))
   iqlist_out = iqlist
  
end subroutine find_q_neighbor


end module boltz_utils
