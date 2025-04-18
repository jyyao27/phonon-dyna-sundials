!!===============================================================================
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

module phonon_dispersion
   use pert_const, only: dp, czero, twopi, ci, cone, e2, pi
   use pert_utils, only: hmat_upper_diag, get_exp_ikr, abs2
   use wigner_seitz_cell, only: ws_cell
   use triplet_cell, only:trip_cell
   use polar_correction, only: set_polar_parameter, dyn_mat_longrange
   use force_constant, only: lattice_ifc, force_const, set_ws_cell_ph
   use force_constant_thirdorder, only: lattice_ifct, force_const_third, set_tripcell_ph, expikr_set, expikr
   use qe_mpi_mod, only: ionode, mp_bcast, inter_pool_comm, ionode_id, ionode, stdout, mp_barrier
   implicit none
   public

   type(lattice_ifc), save :: phon
   type(lattice_ifct), save :: phont
   
   public :: init_lattice_ifc
   public :: init_lattice_ifct
   public :: solve_phonon_modes
   public :: solve_phonon_velocity
   public :: solve_phonon_velocity_findiff
   public :: solve_phi3mat
   public :: solve_phi3mat_fast
   public :: trim_lattice_ifct
   public :: solve_gruneisen
contains

subroutine init_lattice_ifc(file_id, qdim, ph)
   use hdf5_utils
   use epr_hdf5_io, only: read_force_constant
   use pert_data,  only: zstar, epsil, nat, at, bg, &
                           tau, volume, tpiba, thickness_2d, loto_alpha, lpolar

   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: qdim(3)
   type(lattice_ifc), intent(out) :: ph
   !
   integer :: i, nelem
   real(dp), allocatable :: cryst_tau(:,:)

   allocate( cryst_tau(3,nat) )
   !transform atomic position to crystal coordinate
   cryst_tau(:,:) = tau(:,:)
   call cryst_to_cart(nat, cryst_tau, bg, -1)

   call set_ws_cell_ph(ph, qdim, nat, at, cryst_tau)
   if(ionode) call read_force_constant(file_id, ph)

   !bcast
   nelem = ph%na * (ph%na + 1) / 2
   do i = 1, nelem
      call mp_bcast(ph%phi(i)%ifc(:,:,:), ionode_id, inter_pool_comm)
   enddo

   ph%lpol = lpolar
   ph%l_2d = merge(.true., .false., thickness_2d > 0.0_dp)
   !setup for polar correction.
   !N.B. to be consistent with rigid_bulk, alpha=1.0 is enforced
   if(ph%lpol) call set_polar_parameter(qdim, nat, volume, tpiba, bg, &
                        epsil, zstar, tau, ph%pol, thickness_2d, loto_alpha)
   
   deallocate( cryst_tau )
end subroutine init_lattice_ifc

subroutine init_lattice_ifct(file_id, qdim, pht, trim_pht)
   use hdf5_utils
   use epr_hdf5_io, only: read_force_constant_thirdorder
   use pert_data,  only: epsil, nat, at, bg, tau

   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: qdim(3)
   type(lattice_ifct), intent(out) :: pht
   logical, intent(in), optional :: trim_pht
   !
   integer :: i, nelem
   real(dp), allocatable :: cryst_tau(:,:)
   logical :: trim_pht_local

   if (.not. present(trim_pht)) then
      trim_pht_local = .false. 
   endif

   allocate( cryst_tau(3,nat) )
   !transform atomic position to crystal coordinate
   cryst_tau(:,:) = tau(:,:)
   call cryst_to_cart(nat, cryst_tau, bg, -1)

   call set_tripcell_ph(pht, qdim, nat, at, cryst_tau)
   if(ionode) call read_force_constant_thirdorder(file_id, pht)
!   nelem = (pht%na**3+3*pht%na**2+2*pht%na)/6
   nelem = pht%na**3
   do i = 1, nelem
      call mp_bcast(pht%phit(i)%ifct(:,:,:,:), ionode_id, inter_pool_comm)
      call mp_bcast(pht%phit(i)%trip_ph%rvect(:), ionode_id, inter_pool_comm)
      call mp_bcast(pht%phit(i)%trip_ph%ndegt(:), ionode_id, inter_pool_comm)
      call mp_bcast(pht%phit(i)%trip_ph%nrt, ionode_id, inter_pool_comm)
   enddo

   ! kelly yao - trim zero pht values to save time in solve phi3mat
   if (trim_pht_local) then
      write(stdout, '(5x, a)') 'Trimming ifcts'
      call trim_lattice_ifct(pht)
      do i = 1, nelem
         call mp_bcast(pht%phit(i)%ifct(:,:,:,:), ionode_id, inter_pool_comm)
         call mp_bcast(pht%phit(i)%trip_ph%rvect(:), ionode_id, inter_pool_comm)
         call mp_bcast(pht%phit(i)%trip_ph%ndegt(:), ionode_id, inter_pool_comm)
         call mp_bcast(pht%phit(i)%trip_ph%nrt, ionode_id, inter_pool_comm)
      enddo
   endif

   call mp_barrier(inter_pool_comm)
   deallocate( cryst_tau )
end subroutine init_lattice_ifct

subroutine trim_lattice_ifct(pht)

   implicit none
   type(lattice_ifct), intent(inout) :: pht

   ! local variables
   integer :: idx, ir, irt, a1, a2, a3, max_nvect
   integer :: idx_count, idx_recount
   logical, allocatable :: large_irs(:)
   type(trip_cell), pointer :: tc
   integer, allocatable :: rvect_copy(:), ndegt_copy(:)
   integer :: nrt_copy
   
   complex(dp) :: ifct(3,3,3)
   !real(dp) :: ifct(3,3,3)
   complex(dp), allocatable :: ifct_zeros(:,:,:,:)
   !real(dp), allocatable :: ifct_zeros(:,:,:,:)
   idx = 0
   max_nvect = pht%max_nrt
   do a3 = 1, pht%na
   do a2 = 1, pht%na
   do a1 = 1, pht%na
      idx = idx + 1
   !   write(*,*) 'idx', idx
      tc => pht%phit(idx)%trip_ph

      nrt_copy = tc%nrt
      allocate( rvect_copy(nrt_copy) ) 
      allocate( ndegt_copy(nrt_copy) ) 
      rvect_copy(:) = tc%rvect(:)
      ndegt_copy(:) = tc%ndegt(:)

      allocate( large_irs(nrt_copy) )
      large_irs = .false.

      allocate( ifct_zeros(3,3,3,nrt_copy) ) 
      ifct_zeros(:,:,:,:) = pht%phit(idx)%ifct(:,:,:,:)

      idx_count = 0
      do ir = 1, nrt_copy 
         ! for each idx, find the large ir indices
         ifct = pht%phit(idx)%ifct(:,:,:,ir)
         if (any(abs(real(ifct)) > 0.0_dp)) then
            large_irs(ir) = .true.
            idx_count = idx_count + 1
         endif
      enddo
      deallocate( pht%phit(idx)%ifct )
      deallocate( pht%phit(idx)%trip_ph%rvect )
      deallocate( pht%phit(idx)%trip_ph%ndegt )

      ! allocate new arrays for ifct(:,:,:,ir) and tc%rvect(ir)
      allocate( pht%phit(idx)%ifct(3,3,3,idx_count) )
      allocate( pht%phit(idx)%trip_ph%rvect(idx_count) )
      allocate( pht%phit(idx)%trip_ph%ndegt(idx_count) )
      pht%phit(idx)%trip_ph%nrt = idx_count
      idx_recount = 1
      do ir = 1, nrt_copy
         ! for each idx, find the large ir indices
         if ( large_irs(ir) ) then
            pht%phit(idx)%ifct(:,:,:,idx_recount) = ifct_zeros(:,:,:,ir)
            pht%phit(idx)%trip_ph%rvect(idx_recount) = rvect_copy(ir)
            pht%phit(idx)%trip_ph%ndegt(idx_recount) = ndegt_copy(ir)

            idx_recount = idx_recount + 1
         endif            
      enddo
      deallocate( large_irs )
      deallocate( rvect_copy, ndegt_copy, ifct_zeros)

   enddo; enddo; enddo
   

end subroutine trim_lattice_ifct



subroutine solve_phonon_modes(ph, xqt, pheig, phmode)
   use pert_data,  only: mass
   implicit none
   type(lattice_ifc), intent(in) :: ph
   real(dp), intent(in) :: xqt(3)
   real(dp), intent(out) :: pheig(:) !pheig(ph%nm)
   complex(dp), intent(out), optional :: phmode(:,:) !phmode(ph%nm, ph%nm)
   !local variables
   integer :: i, j , ia, ja, ii, jj, m, n, ir, nmodes, nelem, irp
   real(dp) :: mfactor, w2( 3*ph%na ), tmp( 3*ph%na )
   complex(dp), allocatable :: dmat_lr(:,:,:), dmat(:,:,:), ev(:,:), dyn_upper(:), exp_ikr(:)
   type(ws_cell), pointer :: ws
   
   nmodes = 3 * ph%na
   nelem = ph%na * (ph%na + 1) / 2
   !
   allocate( dmat(3,3,nelem), dyn_upper( (nmodes*(nmodes+1))/2 ), exp_ikr(ph%nrvec) )
   !
   call get_exp_ikr(xqt, ph%rvec_set, exp_ikr)
   dmat = czero
   do m = 1, nelem
      ws => ph%phi(m)%ws_ph

      do ir = 1, ws%nr
         irp = ws%rvec(ir)
         ! the 1/ws%ndeg(ir) is already included in ifc(:) !/ real(ws%ndeg(ir), dp)
         !dmat(:,:, m) = dmat(:,:, m) + exp(iqR) * ph%phi(m)%ifc(:,:, ir)
         ! kelly yao - same as
         !dmat(:,:, m) = dmat(:,:, m) + exp_ikr(irp) * ph%phi(m)%ifc(:,:, ir)
         call zaxpy(9, exp_ikr(irp), ph%phi(m)%ifc(1,1,ir), 1, dmat(1,1,m), 1)
      enddo
   enddo
   
   !apply polar correction if needed.
   if( ph%lpol ) then
      allocate( dmat_lr(3, 3, nelem) )
      call dyn_mat_longrange(ph%pol, xqt, dmat_lr)
      dmat = dmat + dmat_lr
      deallocate(dmat_lr)
   endif
   
   dyn_upper = czero
   !create upper trianglar of the dynamical matrix
   do jj = 1, nmodes
   do ii = 1, jj
      !index of the upper triangluatr dynamical matrix
      n = (jj * (jj-1)) / 2 + ii
      
      !atomic index
      ia = (ii - 1) / 3 + 1
      ja = (jj - 1) / 3 + 1
      !cartesian directions
      i = mod(ii-1, 3) + 1
      j = mod(jj-1, 3) + 1
      
      m = ( ja * (ja-1) ) / 2 + ia
      mfactor = 1.0_dp / sqrt( mass(ia) * mass(ja) )

      if(ia .ne. ja) then
         dyn_upper(n) = dmat(i,j,m) * mfactor
      else
         !enforce herminicity
         dyn_upper(n) = (dmat(i,j,m) + conjg(dmat(j,i,m))) * 0.5_dp * mfactor
      endif
   enddo; enddo
   
   if( present(phmode) ) then
      allocate( ev(nmodes, nmodes) )
      call hmat_upper_diag(dyn_upper, nmodes, w2, ev)
      !return eigen-displacement if requested.
      do ii = 1, nmodes
         ia = (ii - 1) / 3 + 1
         phmode(ii, 1:nmodes) = ev(ii, 1:nmodes) / sqrt( mass(ia) )
      enddo
      deallocate( ev )
   else
      call hmat_upper_diag(dyn_upper, nmodes, w2)
   endif
   !compute phonon frequencies: omega = sqrt(w2)
   tmp = sqrt( abs(w2) )
   ! if w2(i) < 0, return negative frequency.
   pheig(1:nmodes) = merge(tmp, -tmp, w2 > 0.0_dp)
   
   !write(*,*) xqt
   !write(*,*) pheig
   !do ii = 1, nmodes
   !   write(*,*) phmode(ii,:)
   !enddo
   deallocate(dmat, dyn_upper)
end subroutine solve_phonon_modes

! this is for non-polar only - Kelly Yao
! for polar materials, contribution from the polar correction is needed either analytically or numerically
subroutine solve_phonon_velocity(ph, xqt, eigvalue, eigvector, velocity)
   use band_structure, only: calc_deleig_a
   use pert_data, only: at
   implicit none
   type(lattice_ifc), intent(in) :: ph
   real(dp), intent(in) :: xqt(3)
   real(dp), intent(in) :: eigvalue(ph%nm) 
   complex(dp), intent(in) :: eigvector(ph%nm,ph%nm)
   real(dp), intent(inout) :: velocity(3, ph%nm)

   ! local variables
   integer :: m, i, j, ii, jj, ia, ja, ir, irp, imode
   real(dp) :: del_eig(ph%nm)
   complex(dp) :: hamq(ph%nm, ph%nm, 3)
   complex(dp), allocatable :: exp_iqr(:)
   real(dp) :: eig_approx(ph%nm)
   real(dp), allocatable :: rset_cart(:)
   type(ws_cell), pointer :: ws

   ! check if this is non polar
   if( ph%lpol ) call errore('solve_phonon_velocity', 'only available for nonpolar materials')
   
   ! take in eigen results from solving phonon modes
   hamq = cmplx(0.0_dp, 0.0_dp, kind=dp)

   allocate( exp_iqr(ph%nrvec), rset_cart(3) )
   exp_iqr = cmplx(0.0_dp, 0.0_dp, kind=dp)

   call get_exp_ikr(xqt, ph%rvec_set, exp_iqr)
    
   m = 0
   do jj = 1, ph%nm
   do ii = 1, jj
      !atomic index
      ia = (ii - 1) / 3 + 1
      ja = (jj - 1) / 3 + 1
      !cartesian directions
      i = mod(ii-1, 3) + 1
      j = mod(jj-1, 3) + 1
      
      m = ( ja * (ja-1) ) / 2 + ia
   
      ws => ph%phi(m)%ws_ph

      do ir = 1, ws%nr
         irp = ws%rvec(ir)
         ! the 1/ws%ndeg(ir) is already included in ifc(:) !/ real(ws%ndeg(ir), dp)
         !dmat(:,:, m) = dmat(:,:, m) + exp(iqR) * ph%phi(m)%ifc(:,:, ir)
         rset_cart(:) = ph%rvec_set(:,irp)
         call cryst_to_cart(1, rset_cart, at, 1)
         hamq(ii,jj,:) = hamq(ii,jj,:) + &
            ph%phi(m)%ifc(i,j,ir) * exp_iqr(irp) * ci * rset_cart 
      enddo

      if (i .ne. j) then
         hamq(jj,ii,:) = conjg( hamq(ii,jj,:))
      else
         ! hamq is real when i = j
         hamq(jj,ii,:) = cmplx(real(hamq(ii,jj,:)), 0.0_dp, kind=dp)
      endif
   enddo; enddo
   
   eig_approx = 0.0_dp
   do ia = 1, 3
      call calc_deleig_a(hamq(:,:,ia), ph%nm, eigvalue, eigvector, del_eig)
      ! appoximate small eigenvalues smaller than 1E-12
      do imode = 1, ph%nm
         eig_approx(imode) = max(1.0E-12_dp, eigvalue(imode))
      enddo
      velocity(ia, :) = del_eig / 2.0_dp / eig_approx(:)
   enddo

   deallocate(exp_iqr)
end subroutine solve_phonon_velocity


subroutine solve_phonon_velocity_findiff(iq0, enk, num_bvec, neighbors, weight, bvec, velocity)
   use pert_const, only : dp
   implicit none
   ! q index, base from 1
   integer, intent(in) :: iq0
   ! energy of different modes
   real(dp), intent(in) :: enk(:,:)
   ! number of b vectors
   integer, intent(in) :: num_bvec
   ! neighbors indices 
   integer, intent(in) :: neighbors(:)
   ! weight
   real(dp), intent(in) :: weight(:)
   ! bvecs
   real(dp), intent(in) :: bvec(:,:)
   !syp - TODO change it from inout to out
   real(dp), intent(inout) :: velocity(:,:)

   ! local variables
   integer :: nmode, nqpts, im, ib
   real(dp) :: omegadiff 
   nmode = size(enk, 1)
   nqpts = size(enk, 2)
   do im = 1, nmode
      do ib = 1, num_bvec
         if ((neighbors(ib) .gt. 0) .and. (neighbors(ib) .le. nqpts)) then
            omegadiff = enk(im, neighbors(ib)) - enk(im,iq0)
         else
            omegadiff = 0.0_dp
         endif
         !syp - 20240514 see boltz_grid_neighbors.f90
         !bvec: cartesian coordinate, in the unit of tpiba/bohr
         !bweight: in the unit of (tpiba/bohr)^-2 (bweight*bvec*bvec -> dimensionless)
         !so bweight*bvec is in the unit of (tpiba/bohr)^-1 (df/dk is in (tpiba/bohr)^-1)
         !So wherever you use phonon velocity, please divided by two pi.
         !syp - 20241012: Now I decide to devide by 2pi here
         velocity(:, im) = velocity(:, im) + weight(ib) * bvec(:,ib) * omegadiff / (2.0_dp * pi)
         !velocity(:, im) = velocity(:, im) + qg%bweight(ib) * qg%bvec(:,ib) * omegadiff
      enddo
   enddo
end subroutine solve_phonon_velocity_findiff

subroutine solve_phi3mat(pht,xkk,xqp,xkq,wf1,wf2,wf3,uf1,uf2,uf3,jb,im,ib,ph3mat)
   implicit none
   integer, intent(in) :: jb, ib, im
   real(dp), intent(in) ::xkk(3),xqp(3),xkq(3)
   real(dp), intent(in) :: wf1(:),wf2(:),wf3(:)
   complex(dp), intent(in) ::uf1(:,:),uf2(:,:),uf3(:,:)
   !
   type(lattice_ifct), intent(in) ::pht
   type(trip_cell), pointer :: tc
   complex(dp), intent(out) :: ph3mat
  !local parameter
   integer :: i, j, k, k1, k2, k3, a1, a2, a3, ir, irt, idx
   real(dp) :: wq(3), invomegaprod
   complex(dp) :: ifct(3,3,3)
   complex(dp) :: expiqr
   complex(dp), dimension(:,:,:,:,:,:), allocatable :: egvprod
   complex(dp), allocatable :: egv(:,:), exp_ikr1(:), exp_ikr2(:)


   ph3mat=0.0_dp
   wq(1) = wf1(jb)
   wq(2) = wf2(im)
   wq(3) = wf3(ib)
   ! Frequency product
   if (wq(1)<1.0E-7_dp .or. wq(2)<1.0E-7_dp .or. wq(3)<1.0E-7_dp) return
   invomegaprod=1/sqrt(wq(1)*wq(2)*wq(3))
   ! The actual matrix element

   ! reduced mass in units of amass already contained
   ! in units of 1/sqrt(mass)
   allocate(egv(pht%na*3, 3))

   egv(:,1) = uf1(:,jb)
   egv(:,2) = uf2(:,im)
   egv(:,3) = uf3(:,ib)

   allocate(egvprod(3,3,3,pht%na,pht%na,pht%na))
   ! Eigenvector product
   do a3=1,pht%na
   do a2=1,pht%na
   do a1=1,pht%na
       do k=1,3
       do j=1,3
       do i=1,3
           k1=(a1-1)*3+i
           k2=(a2-1)*3+j
           k3=(a3-1)*3+k
           egvprod(i,j,k,a1,a2,a3)=egv(k1,1)*egv(k2,2)*egv(k3,3)
       enddo; enddo; enddo
   enddo; enddo; enddo

   allocate(exp_ikr1(pht%nrvect), exp_ikr2(pht%nrvect))
   call get_exp_ikr(xqp, pht%rvect_set(:,1,:), exp_ikr1)   
   call get_exp_ikr(xkq, pht%rvect_set(:,2,:), exp_ikr2)   
   !
   idx = 0
   do a3 = 1, pht%na
   do a2 = 1, pht%na
   do a1 = 1, pht%na
      idx = idx + 1
      tc => pht%phit(idx)%trip_ph
      !write(*,*) tc%nrt
         do ir = 1, tc%nrt
            ifct = pht%phit(idx)%ifct(:,:,:,ir)
            !if (norm2(ifct) == 0) cycle            
            irt = tc%rvect(ir)
            expiqr = exp_ikr1(irt)*exp_ikr2(irt)
            ph3mat=ph3mat+sum(ifct*egvprod(:,:,:,a1,a2,a3))*expiqr
         enddo
   enddo; enddo; enddo

   ph3mat=ph3mat*invomegaprod/sqrt(2.0_dp)**3/6.0_dp

end subroutine solve_phi3mat


subroutine solve_phi3mat_rank3(pht,xkk,xqp,xkq,wf1,wf2,wf3,uf1,uf2,uf3,jb,im,ib,ph3mat)
   implicit none
   integer, intent(in) :: jb, ib, im
   real(dp), intent(in) ::xkk(3),xqp(3),xkq(3)
   real(dp), intent(in) :: wf1(:),wf2(:),wf3(:)
   complex(dp), intent(in) ::uf1(:,:),uf2(:,:),uf3(:,:)
   !
   type(lattice_ifct), intent(in) ::pht
   type(trip_cell), pointer :: tc
   complex(dp), intent(out) :: ph3mat(pht%na*3,pht%na*3,pht%na*3)
  !local parameter
   integer :: i, j, k, k1, k2, k3, a1, a2, a3, ir, irt, idx
   real(dp) :: wq(3), invomegaprod
   complex(dp) :: ifct(3,3,3)
   complex(dp) :: expiqr
   complex(dp), dimension(:,:,:,:,:,:), allocatable :: egvprod
   complex(dp), allocatable :: egv(:,:), exp_ikr1(:), exp_ikr2(:)


   ph3mat=0.0_dp
   wq(1) = wf1(jb)
   wq(2) = wf2(im)
   wq(3) = wf3(ib)
   ! Frequency product
   if (wq(1)<1.0E-7_dp .or. wq(2)<1.0E-7_dp .or. wq(3)<1.0E-7_dp) return
   invomegaprod=1/sqrt(wq(1)*wq(2)*wq(3))
   ! The actual matrix element

   ! reduced mass in units of amass already contained
   ! in units of 1/sqrt(mass)
   allocate(egv(pht%na*3, 3))

   egv(:,1) = uf1(:,jb)
   egv(:,2) = uf2(:,im)
   egv(:,3) = uf3(:,ib)

   allocate(egvprod(3,3,3,pht%na,pht%na,pht%na))
   ! Eigenvector product
   do a3=1,pht%na
   do a2=1,pht%na
   do a1=1,pht%na
       do k=1,3
       do j=1,3
       do i=1,3
           k1=(a1-1)*3+i
           k2=(a2-1)*3+j
           k3=(a3-1)*3+k
           egvprod(i,j,k,a1,a2,a3)=egv(k1,1)*egv(k2,2)*egv(k3,3)
       enddo; enddo; enddo
   enddo; enddo; enddo

   allocate(exp_ikr1(pht%nrvect), exp_ikr2(pht%nrvect))
   call get_exp_ikr(xqp, pht%rvect_set(:,1,:), exp_ikr1)   
   call get_exp_ikr(xkq, pht%rvect_set(:,2,:), exp_ikr2)   
   !
   idx = 0
   do a3 = 1, pht%na
   do a2 = 1, pht%na
   do a1 = 1, pht%na
      idx = idx + 1
      tc => pht%phit(idx)%trip_ph
      !write(*,*) tc%nrt
         do ir = 1, tc%nrt
            ifct = pht%phit(idx)%ifct(:,:,:,ir)
            !if (norm2(ifct) == 0) cycle            
            irt = tc%rvect(ir)
            expiqr = exp_ikr1(irt)*exp_ikr2(irt)
            ph3mat=ph3mat+sum(ifct*egvprod(:,:,:,a1,a2,a3))*expiqr
         enddo
   enddo; enddo; enddo

   ph3mat=ph3mat*invomegaprod/sqrt(2.0_dp)**3/6.0_dp

end subroutine solve_phi3mat_rank3


subroutine solve_phi3mat_fast(pht,expiqr,wq1,wq2,wq3,uf1,uf2,uf3,ph3mat)
   implicit none
   real(dp), intent(in) :: wq1,wq2,wq3
   complex(dp), intent(in) ::uf1(:),uf2(:),uf3(:)
   complex(dp), intent(in) :: expiqr(:)
   !complex(dp), intent(in) :: exp_ikr1(:), exp_ikr2(:)
   !
   type(lattice_ifct), intent(in) ::pht
   type(trip_cell), pointer :: tc
   complex(dp), intent(out) :: ph3mat
  !local parameter
   integer :: i, j, k, k1, k2, k3, a1, a2, a3, ir, idx
   integer :: nrt, irt
   real(dp) :: invomegaprod, omegaprod
   complex(dp) :: ifct(3,3,3)
   complex(dp), dimension(:,:,:,:,:,:), allocatable :: egvprod
   !complex(dp), allocatable :: exps_vec(:), inc_arr(:)
   complex(dp) :: inc


   ph3mat=0.0_dp
   ! Frequency product
   omegaprod = wq1*wq2*wq3
   !if (wq1<1.0E-7_dp .or. wq2<1.0E-7_dp .or. wq3<1.0E-7_dp) return
   if (omegaprod < 1.0E-16_dp) return
   invomegaprod=1/sqrt(wq1*wq2*wq3)
   ! The actual matrix element
    
   ! reduced mass in units of amass already contained
   ! in units of 1/sqrt(mass)

   allocate(egvprod(3,3,3,pht%na,pht%na,pht%na))
   egvprod = 0.0_dp
   ! Eigenvector product
   
   do a3=1,pht%na
   do a2=1,pht%na
   do a1=1,pht%na
      do k=1,3
      do j=1,3
      do i=1,3
         k1=(a1-1)*3+i
         k2=(a2-1)*3+j
         k3=(a3-1)*3+k
         egvprod(i,j,k,a1,a2,a3)=uf1(k1)*uf2(k2)*uf3(k3)
         !egvprod(i,j,k,a1,a2,a3)=egv(k1,1)*egv(k2,2)*egv(k3,3)
      enddo; enddo; enddo
   enddo; enddo; enddo

   !
   idx = 0
   do a3 = 1, pht%na
   do a2 = 1, pht%na
   do a1 = 1, pht%na
      idx = idx + 1
      tc => pht%phit(idx)%trip_ph
      nrt = tc%nrt
      !allocate( inc_arr(nrt))
      
      do ir = 1, nrt
         ifct = pht%phit(idx)%ifct(:,:,:,ir)
         !inc_arr(ir) = sum(ifct*egvprod(:,:,:,a1,a2,a3))
         irt = tc%rvect(ir)
         inc = sum(ifct*egvprod(:,:,:,a1,a2,a3))*expiqr(irt)
         ph3mat = ph3mat + inc
      enddo

      !ph3mat = ph3mat + dot_product(inc_arr,exps_vec)
      !deallocate(inc_arr)

   enddo; enddo; enddo

   ph3mat=ph3mat*invomegaprod/sqrt(2.0_dp)**3/6.0_dp

end subroutine solve_phi3mat_fast

subroutine solve_gruneisen (pht, xk, wf, uf, grun)
   use pert_data,  only: nat, tau, bg
   use pert_const, only: ryd2mev
   implicit none
   real(dp), intent(in) ::xk(3), wf(:)
   complex(dp), intent(in) ::uf(:,:)
   !
   type(lattice_ifct), intent(in) ::pht
   type(trip_cell), pointer :: tc
   real(dp), intent(out) :: grun(:)
   !local parameter
   integer :: nu, i, j, k, a1, a2, a3, ir, irt, idx
!   integer :: k, k1, k2, k3
   real(dp) :: wfnu, omegaprod, rv3(3)
   complex(dp) :: ifct(3,3,3)
   complex(dp) :: grutmp
   complex(dp) :: expiqr, evprod
   complex(dp), allocatable :: egva(:,:,:), egv(:), exp_ikr(:)
   real(dp), allocatable :: cryst_tau(:,:)


   allocate( cryst_tau(3,pht%na) )
   !transform atomic position to crystal coordinate
   cryst_tau(:,:) = tau(:,:)
   call cryst_to_cart(nat, cryst_tau, bg, -1)

   allocate(egva(3,pht%na,pht%na*3), egv(pht%na*3))
   allocate(exp_ikr(pht%nrvect))
   egva = czero; egv=czero; exp_ikr =czero;
   call get_exp_ikr(xk, pht%rvect_set(:,1,:), exp_ikr)

do nu = 1, pht%na*3
   grun(nu) = 0.0_dp

   wfnu = wf(nu)
   if (wfnu < 0.0001_dp) cycle
   omegaprod = -1.0_dp/(6.0_dp*wfnu**2)
   egv = uf(:,nu)
   grutmp = czero  
   
   do a1=1,pht%na; do i=1,3
      egva(i,a1,nu)=egv((a1-1)*3+i)
   enddo; enddo

   idx = 0

   do a3 = 1, pht%na
   do a2 = 1, pht%na
   do a1 = 1, pht%na
      idx = idx + 1
      tc => pht%phit(idx)%trip_ph
      do ir = 1, tc%nrt
         irt = tc%rvect(ir)
         ifct = pht%phit(idx)%ifct(:,:,:,ir)
         !if (norm2(ifct) == 0) cycle
         do i = 1,3; do j = 1,3;
            evprod = conjg(egva(i,a1,nu))*egva(j,a2,nu)
         do k = 1,3
            rv3 = pht%rvect_set(:,2,irt)+cryst_tau(:,a3) - cryst_tau(:,a1)
            ! from crystal to cartesian coordinates
            call cryst_to_cart(nat, rv3, bg, 1)
            expiqr = exp_ikr(irt)
            grutmp=grutmp+ifct(i,j,k)*expiqr*rv3(k)*evprod
 

         enddo ! k
         enddo; enddo ! i,j
       enddo ! ir
   enddo; enddo; enddo

      grutmp = grutmp*omegaprod
      grun(nu) = real(grutmp,dp)
!
enddo
   deallocate(egva)

end subroutine solve_gruneisen

end module phonon_dispersion
