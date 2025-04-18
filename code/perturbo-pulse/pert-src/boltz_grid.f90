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
!   define data used in boltzmann dynamics and transport calculations
!
!  - at a given kmesh defined by kdim, select and output tetrahedra that
!    contributed to BZ integration. 
!  - collect vertex of all selected tetrahedra, but only output the 
!    irreducible into tet.kpt file.
!  - only levels between [bmin, bmax] are considered.  only resluts (dos) 
!    inside [emin, emax] are wanted.
!
! Maintenance:
!===============================================================================

module boltz_grid
   use pert_const, only: dp, ryd2mev, ryd2ev, pi
   implicit none
   private
 
   type, public :: grid
      !k-point
      integer :: ndim(3) ! size of the grid
      integer :: nk  !number of selected kpoints, includes reducible kpoints.
      integer :: nk_irr !number of irreducible k-points
      !tetrahedra
      !number of tetrahedra, and weight of each tetrahedra
      integer  :: num_tet
      real(dp) :: tweight
      ! tet(2,4,num_tet): tet(1,:,:)->irreducible k; tet(2,:,:)->full grid
      integer, pointer :: tetra(:,:,:)=>null()
      !
      real(dp):: kweight !weight of each k-point, = 1/(kdim(1)*kdim(2)*kdim(3))
      !
      !index representation of kpoints in full grid
      integer, pointer :: kpt(:)=>null() ! kpt(nk)
      !irreducible k-list, irrk(nk_irr)
      !the symmetry operator index in symop_inv_cart to transform the 
      !irreducible k-point to the current reducible k-point
      integer, pointer :: symopindex(:)=>null() ! symopindex(nk)
      ! store index representation of this kpoints
      integer, pointer :: irrk(:)=>null()  !irrk(nk_irr)
      !map the i-th kpt in full grid to its irreducible kpt 
      ! irrk( kpt2ir(i) ): i-th kpoints in full grid -> irr-kpoints.
      integer, pointer :: kpt2ir(:)=>null() ! kpt2ir(nk)
      !
      !syp - map the i-th irreducible kpt to its index in total selected reducible kpoints
      !      counting from 1
      !      for example, if we have 
      !          kpt(1:10)    = [0, 5, 6, 7,10,21,22,34,36,51]
      !          kpt2ir(1:10) = [1, 2, 2, 3, 2, 4, 2, 3, 4, 4]
      !          irrk(1:4)    = [0, 5,    7,   21]
      !          irr2kpt(1:4) = [1, 2,    4,    6] 
      integer, pointer :: ir2kpt(:)=>null() ! ir2kpt(nk_irr)
      !
      !band structure and band velocity information
      ! kgrid_enk(numb, nk)
      real(dp), pointer :: enk(:,:)=>null()
      ! E_nk for irreducible k-points in irrk.
      real(dp), pointer :: enk_irr(:,:)=>null()
      ! band velocity, vnk(3,numb,nk)
      real(dp), pointer :: vnk(:,:,:)=>null()
      ! only bands in [bmin, bmax] contribute to transport
      ! 1 < bmin < bmax < num_bands;  numb = bmax - bmin + 1
      integer :: numb, bmin, bmax

      !neighbor k-points, for kspace derivatives
      !using the finite-difference formula in wannier90,
      !Refer to Mostofi et.al comp. phys. commu. 178 (2008) 685-699
      integer :: num_neighbors   ! number of neighbor
      !b vectors in cartesian coordinates, unit tpiba/bohr.
      real(dp), pointer :: bvec(:,:)=>null()
      !weight used to compute df/dk, unit (tpiba/bohr)^-2
      real(dp), pointer :: bweight(:)=>null()
      !neighbors(i,ik): index of i-th neighbor of the ik-points.
      integer,  pointer :: neighbors(:,:)=>null()
   end type grid

   !real(dp), parameter :: eig_tol = 20.0E0_dp/ryd2mev ! 50 meV

   public :: init_boltz_grid
   public :: init_boltz_qgrid
   public :: init_boltz_grid_sym
   public :: boltz_grid_generate
   public :: boltz_grid_load
   public :: sum_equivalent_kpoint
   !syp - added
   public :: init_boltz_qgrid_vel 
   public :: init_boltz_qgrid_vel_sym
   public :: symcheck_absvbs_ongrid
   public :: symmetrize_absvbs_ongrid
contains

subroutine init_boltz_grid(kg, el, bmin, bmax, load)
   use boltz_utils, only: num2kpt
   use qe_mpi_mod, only: mp_split_pools, mp_sum, inter_pool_comm
   use band_structure, only: electron_wann, solve_band_velocity
   implicit none
   type(electron_wann), intent(in) :: el
   integer, intent(in) :: bmin, bmax
   type(grid), intent(inout) :: kg !kgrid
   !if false, then kg is already initialized, no need to load from files.
   logical, intent(in), optional :: load
   !local variable
   logical :: load_grid
   integer :: ik, kst, kend
   real(dp) :: xk(3), eig(el%nb), vel(3,el%nb)

   !initialize numb
   kg%numb = bmax - bmin + 1
   kg%bmin = bmin
   kg%bmax = bmax
   !
   !read input files, init nk, nk_irr, kweight, irrk_k,
   load_grid = .true.
   if( present(load) ) load_grid = load
   if( load_grid ) call boltz_grid_load(kg, 'k')
   !
   !calc. eig and velocity
   allocate(kg%enk(kg%numb, kg%nk), kg%vnk(3, kg%numb, kg%nk), kg%symopindex(kg%nk))

   kg%enk = 0.0_dp;  kg%vnk = 0.0_dp; kg%symopindex = 0
   call mp_split_pools(kg%nk, kst, kend)
!$omp parallel do schedule(guided) default(shared) private(ik, xk, eig, vel)
   do ik = kst, kend
      xk = num2kpt(kg%kpt(ik), kg%ndim)
      call solve_band_velocity(el, xk, vel, eigvalue=eig)
      kg%enk(:,ik) = eig(bmin:bmax)
      kg%vnk(:,:,ik) = vel(:, bmin:bmax)
   enddo
!$omp end parallel do
   call mp_sum(kg%enk, inter_pool_comm)
   call mp_sum(kg%vnk, inter_pool_comm)
end subroutine init_boltz_grid

!> Initialize q-grid for phonons
subroutine init_boltz_qgrid(qg, ph, load, finegrid)
   !velocity is not computed at this point
   use boltz_utils, only: num2kpt
   use qe_mpi_mod, only: mp_split_pools, mp_sum, inter_pool_comm
   use force_constant, only: lattice_ifc
   use phonon_dispersion, only: solve_phonon_modes, solve_phonon_velocity
   implicit none
   type(lattice_ifc), intent(in)   :: ph     !< Wannier basis inter-atomic force constants for interpolation
   type(grid), intent(inout) :: qg           !< q-grid 
   logical, intent(in), optional :: finegrid
   logical, intent(in), optional :: load     !< load flag, if false, grid initialized, no need to load from files.
   logical :: load_grid                      !< local variable for loading q-grid
   logical :: fine                           !< local variable for loading q-grid
   integer :: iq                             !< q-index
   integer :: qst                            !< starting q-index 
   integer :: qend                           !< ending q-index
   real(dp) :: xq(3)                         !< q vector in lattice coordinates
   real(dp) ::  pheig(ph%nm)                 !< phonon frequencies
   complex(dp) :: phmode(ph%nm,ph%nm)        !< phonon mode matrix at a given q

   !initialize numb
   qg%numb = ph%nm
   qg%bmin = 1
   qg%bmax = ph%nm
   !
   !read input files, init nk, nk_irr, kweight, irrk_k,
   load_grid = .true.
   fine = .false.
   if( present(load) ) load_grid = load
   if (present(finegrid)) fine = finegrid
   if( load_grid ) then
      if (fine) then
         call boltz_grid_load(qg, 'f')
      else
         call boltz_grid_load(qg, 'q')
      endif
   endif
   !calc. pheig
   allocate(qg%enk(qg%numb, qg%nk), qg%vnk(3, qg%numb, qg%nk), qg%symopindex(qg%nk))

   qg%enk = 0.0_dp; qg%symopindex = 0
   call mp_split_pools(qg%nk, qst, qend)
!$omp parallel do schedule(guided) default(shared) private(iq, xq, pheig, phmode)
   do iq = qst, qend
      xq = num2kpt(qg%kpt(iq), qg%ndim)
      call solve_phonon_modes(ph, xq, pheig, phmode)
      qg%enk(:,iq) = pheig(qg%bmin:qg%bmax)
   enddo
!$omp end parallel do
   call mp_sum(qg%enk, inter_pool_comm)
end subroutine init_boltz_qgrid


!after execute init_boltz_qgrid, the velocity is not computed yet.
!the velocity will be computed in this subroutine.
!syp - compute phonon velocity: analytical or numerical 
!      1. analytical: no need to use neighboring q-points
!      2. numerical: need to use neighboring q-points, should execute 
!            after "set_grid_neighbors" outside this function, that means 
!            this subrtoutine will be executed twice
subroutine init_boltz_qgrid_vel(qg, ph, ph_vel_numerical)
   use boltz_utils, only: num2kpt
   use qe_mpi_mod, only: mp_split_pools, mp_sum, inter_pool_comm
   use force_constant, only: lattice_ifc
   use phonon_dispersion, only: solve_phonon_modes, solve_phonon_velocity, &
      solve_phonon_velocity_findiff
   implicit none
   type(lattice_ifc), intent(in) :: ph
   type(grid), intent(inout) :: qg !qgrid
   !syp - whether to adaptively smear the phonon-phonon interaction
   logical, intent(in),optional :: ph_vel_numerical
   !local variable
   integer :: iq, qst, qend
   real(dp) :: xq(3), vel(3,ph%nm), pheig(ph%nm)
   complex(dp) :: phmode(ph%nm,ph%nm)
   logical :: ph_vel_numerical_local
 
   !initialize adapt_smear_phph
   ph_vel_numerical_local = .false.
   if( present(ph_vel_numerical) ) ph_vel_numerical_local = ph_vel_numerical

   allocate(qg%vnk(3, qg%numb, qg%nk))
   !phonon velocity is not considered right now
   qg%vnk = 0.0_dp
   call mp_split_pools(qg%nk, qst, qend)
!$omp parallel do schedule(guided) default(shared) private(iq, xq, phmode, pheig, vel)
   do iq = qst, qend
      xq = num2kpt(qg%kpt(iq), qg%ndim)
      ! syp - [TODO] this line is executed twice
      !       should be optimized later
      if (ph_vel_numerical_local) then
         call solve_phonon_velocity_findiff(iq, qg%enk, qg%num_neighbors, &
             qg%neighbors(:,iq), qg%bweight, qg%bvec, vel)
         ! divided by 2*pi: solve_phonon_velocity_findiff output the phonon velocity in 2*pi/(alat*bohr)
         ! the normal unit for velocity is 1/(alat*bohr)
        !vel = vel / (2.d0 * pi)
      else
         call solve_phonon_modes(ph, xq, pheig, phmode)
         call solve_phonon_velocity(ph, xq, qg%enk(:,iq), phmode, vel)
      endif 
      qg%vnk(:,:,iq) = vel
   enddo
!$omp end parallel do
   call mp_sum(qg%vnk, inter_pool_comm)
end subroutine init_boltz_qgrid_vel


!compute enk and vel for irreducible qpoint only, 
!those of reducible k are determined from symmetry.
subroutine init_boltz_qgrid_vel_sym(qg, ph, ph_vel_numerical)
   use boltz_utils, only: num2kpt, kpt2num
   use pert_data, only: nsym, symop, at, bg, symop_inv_cart
   use qe_mpi_mod, only: mp_split_pools, mp_sum, inter_pool_comm, my_pool_id, ionode, stdout
   use force_constant, only: lattice_ifc
   use phonon_dispersion, only: solve_phonon_modes, solve_phonon_velocity, &
      solve_phonon_velocity_findiff
   implicit none
   type(lattice_ifc), intent(in) :: ph
   type(grid), intent(inout) :: qg !qgrid
   !syp - whether to adaptively smear the phonon-phonon interaction
   logical, intent(in),optional :: ph_vel_numerical
   !local variable
   integer :: iq, qst, qend, irk, i, j
   real(dp) :: xq(3), xqr(3), vel(3,ph%nm), pheig(ph%nm)
   complex(dp) :: phmode(ph%nm,ph%nm)
   logical :: ph_vel_numerical_local
   !
  !real(dp) :: symop_inv_cart(3,3,nsym)
   real(dp), allocatable :: vnk_irr(:,:,:), vabs2
   ! Identical matrix
   integer, parameter :: eye(3,3) = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3,3/))
 
   !if ph%nm /= qg%numb, invoke an error
   if( ph%nm /= qg%numb ) call errore('init_boltz_qgrid_vel_sym','ph%nm /= qg%numb',1)

   !initialize adapt_smear_phph
   ph_vel_numerical_local = .false.
   if( present(ph_vel_numerical) ) ph_vel_numerical_local = ph_vel_numerical

   allocate(qg%vnk(3, qg%numb, qg%nk))
   allocate(qg%enk_irr(qg%numb, qg%nk_irr))
   allocate( vnk_irr(3, qg%numb, qg%nk_irr) )
   vnk_irr = 0.0_dp
   !phonon velocity is not considered right now
   qg%vnk = 0.0_dp
   qg%enk_irr = 0.0_dp
   call mp_split_pools(qg%nk_irr, qst, qend)
!$omp parallel do schedule(guided) default(shared) private(iq, xq, phmode, pheig, vel)
   do iq =qst, qend
      xq = num2kpt(qg%irrk(iq), qg%ndim)
      ! syp - [TODO] this line is executed twice
      !       should be optimized later
      if (ph_vel_numerical_local) then
         call solve_phonon_velocity_findiff(qg%ir2kpt(iq), qg%enk, qg%num_neighbors, &
             qg%neighbors(:,qg%ir2kpt(iq)), qg%bweight, qg%bvec, vnk_irr(:,:,iq)) !vel)
         ! divided by 2*pi: solve_phonon_velocity_findiff output the phonon velocity in 2*pi/(alat*bohr)
         ! the normal unit for velocity is 1/(alat*bohr)
        !vnk_irr(:,:,iq) = vnk_irr(:,:,iq) / (2.d0 * pi)
      else
         call solve_phonon_modes(ph, xq, pheig, phmode)
         call solve_phonon_velocity(ph, xq, pheig, phmode, vnk_irr(:,:,iq)) !vel)
        !qg%enk_irr(:,iq) = pheig
      endif 
     !vnk_irr(:,:,iq) = vel
   enddo
!$omp end parallel do
   call mp_sum(vnk_irr, inter_pool_comm)

   !compute symop in cartesian coord.
   ! q_cart = bg * q, so q = (bg)^-1 * q_cart; and (bg)^-1 = transpose(at)
   ! if S*q = q'; then S * [(bg)^-1 * q_cart] = [(bg)^-1 * q'_cart]
   ! therefore, [bg * S * (bg)^-1] * q_cart = q'_cart
   symop_inv_cart = 0.0_dp
   do i = 1, nsym
      do j = 1, nsym
         if( all(matmul(symop(:,:,i), symop(:,:,j)) .eq. eye) ) then
            symop_inv_cart(:,:,i) = matmul( matmul(bg, symop(:,:,j)), transpose(at) )
            exit
         endif
      enddo
      !throw an error if not initialized
      if(sum( abs(symop_inv_cart(:,:,i)) ) < 1.0E-12_dp) &
         call errore('init_boltz_grid_sym','failed to compute symop_inv_cart',1)
   enddo

   !syp - map energy to irreducible k-points
   do iq = 1, qg%nk_irr
      qg%enk_irr(:, iq) = qg%enk(:, qg%ir2kpt(iq))
   enddo
   qg%vnk = 0.0_dp
   do iq = 1, qg%nk
      irk = qg%kpt2ir(iq)
      ! E(Sq) = E(q)
      !syp - symmetrization of energy
      qg%enk(:, iq) = qg%enk_irr(:, qg%kpt2ir(iq))
      !find the symmetry operation that transform irk to iq
      xqr = num2kpt( qg%irrk(irk), qg%ndim )
      do i = 1, nsym
         xq = matmul(symop(1:3, 1:3, i), xqr)
         !if find the operation, compute v(Sq)
         ! v(Sq) = v(q) * S^-1:  v(Sq)_j = \sum_{i} v(q)_i * (S^-1)_ij
         if( kpt2num(xq, qg%ndim) .eq. qg%kpt(iq) ) then
            do j = 1, qg%numb
               qg%vnk(:, j, iq) = matmul(vnk_irr(:, j, irk), symop_inv_cart(:,:,i))
            enddo
            qg%symopindex(iq) = i
            exit
         endif
      enddo
      !check
      vabs2 = dot_product(vnk_irr(:,1,irk), vnk_irr(:,1,irk))
      if(abs(dot_product(qg%vnk(:,1,iq), qg%vnk(:,1,iq)) - vabs2) > 1.E-6_dp*vabs2) &
         call errore('init_boltz_qgrid_vel_sym','failed to unfold velocity', 1)
   enddo

   !sanity check
   !the irreducible k-points representative should have the qg%symopindex = 1
   do iq = 1, qg%nk_irr
      if (qg%symopindex(qg%ir2kpt(iq)) /= 1) &
         call errore('init_boltz_qgrid_vel_sym','failed to find symopindex',1)
   enddo
   deallocate(vnk_irr)
end subroutine init_boltz_qgrid_vel_sym

!compute enk and vel for irreducible kpoint only, 
!those of reducible k are determined from symmetry.
subroutine init_boltz_grid_sym(kg, el, bmin, bmax, load)
   use pert_data, only: nsym, symop, at, bg
   use boltz_utils, only: num2kpt, kpt2num
   use qe_mpi_mod, only: mp_split_pools, mp_sum, inter_pool_comm
   use band_structure, only: electron_wann, solve_band_velocity
   implicit none
   type(electron_wann), intent(in) :: el
   !nsym: number of symmetry operation, symop: symmetry operation in crystal coord.
   integer, intent(in) :: bmin, bmax
   type(grid), intent(inout) :: kg !kgrid
   !if false, then kg is already initialized, no need to load from files.
   logical, intent(in), optional :: load
   !local variable
   logical :: load_grid
   integer :: ik, kst, kend, i, j, irk, ib
   real(dp) :: xk(3), eig(el%nb), vel(3,el%nb), symop_inv_cart(3,3,nsym), xkr(3), vabs2
   real(dp), allocatable :: vnk_irr(:,:,:)
   ! Identical matrix
   integer, parameter :: eye(3,3) = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3,3/))

   !initialize numb
   kg%numb = bmax - bmin + 1
   kg%bmin = bmin
   kg%bmax = bmax
   
   !read input files, init nk, nk_irr, kweight, irrk_k, 
   load_grid = .true.
   if( present(load) ) load_grid = load
   if( load_grid ) call boltz_grid_load(kg, 'k')
   
   !calc. eig and velocity
   allocate(kg%enk(kg%numb,kg%nk), kg%vnk(3,kg%numb, kg%nk), kg%enk_irr(kg%numb,kg%nk_irr))

   allocate( vnk_irr(3, kg%numb, kg%nk_irr) )
   kg%enk_irr = 0.0_dp;  vnk_irr = 0.0_dp
   
   call mp_split_pools(kg%nk_irr, kst, kend)
!$omp parallel do schedule(guided) default(shared) private(ik, xk, eig, vel)
   do ik = kst, kend
      xk = num2kpt( kg%irrk(ik), kg%ndim )
      call solve_band_velocity(el, xk, vel, eigvalue=eig)
      kg%enk_irr(:,ik) = eig(bmin:bmax)
      vnk_irr(:,:,ik) = vel(:, bmin:bmax)
   enddo
!$omp end parallel do
   call mp_sum(kg%enk_irr, inter_pool_comm)
   call mp_sum(vnk_irr, inter_pool_comm)

   !compute symop in cartesian coord.
   ! k_cart = bg * k, so k = (bg)^-1 * k_cart; and (bg)^-1 = transpose(at)
   ! if S*k = k'; then S * [(bg)^-1 * k_cart] = [(bg)^-1 * k'_cart]
   ! therefore, [bg * S * (bg)^-1] * k_cart = k'_cart
   symop_inv_cart = 0.0_dp
   do i = 1, nsym
      do j = 1, nsym
         if( all(matmul(symop(:,:,i), symop(:,:,j)) .eq. eye) ) then
            symop_inv_cart(:,:,i) = matmul( matmul(bg, symop(:,:,j)), transpose(at) )
            exit
         endif
      enddo
      !throw an error if not initialized
      if(sum( abs(symop_inv_cart(:,:,i)) ) < 1.0E-12_dp) &
         call errore('init_boltz_grid_sym','failed to compute symop_inv_cart',1)
   enddo

   kg%enk = 0.0_dp;  kg%vnk = 0.0_dp
   do ik = 1, kg%nk
      irk = kg%kpt2ir(ik)
      ! E(Sk) = E(k) 
      kg%enk(:, ik) = kg%enk_irr(:, irk)
      !find the symmetry operation that transform irk to ik
      xkr = num2kpt( kg%irrk(irk), kg%ndim )
      do i = 1, nsym
         xk = matmul(symop(1:3, 1:3, i), xkr)
         !if find the operation, compute v(Sk)
         ! v(Sk) = v(k) * S^-1:  v(Sk)_j = \sum_{i} v(k)_i * (S^-1)_ij
         if( kpt2num(xk, kg%ndim) .eq. kg%kpt(ik) ) then
            do ib = 1, kg%numb
               kg%vnk(:, ib, ik) = matmul(vnk_irr(:, ib, irk), symop_inv_cart(:,:,i))
            enddo
            exit
         endif
      enddo
      !check
      vabs2 = dot_product(vnk_irr(:,1,irk), vnk_irr(:,1,irk)) 
      if(abs(dot_product(kg%vnk(:,1,ik), kg%vnk(:,1,ik)) - vabs2) > 1.E-6_dp*vabs2) &
         call errore('init_boltz_grid_sym','failed to unfold velocity', 1)
   enddo
   
   deallocate(vnk_irr)
end subroutine init_boltz_grid_sym

! assume single process
subroutine symcheck_absvbs_ongrid(kg, twodarray, symmornot)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: twodarray(:,:) ! assume (iband, ikpoint)
   logical, intent(out), optional :: symmornot

   integer :: irk, ik, ib, idx, idy
   real(dp), allocatable :: onedvec(:), onedvec2(:), tmp

   if (present(symmornot)) symmornot = .true.
   !get the dimension of 2darray
   idx = size(twodarray, 1) ! number of bands
   idy = size(twodarray, 2) ! number of reducible k-points

   !allocate work array
   allocate(onedvec(idx), onedvec2(idx))

   ! for each irreducible k-point, check whether its equivalent k-points
   ! have the same absolute value for each branch
   do ik = 1, kg%nk_irr
      onedvec = twodarray(:, kg%ir2kpt(ik)) ! the baseline for ik-th irreducible k-point
      do irk = 1, kg%nk
         if (kg%kpt2ir(irk) .eq. ik) then
            onedvec2 = twodarray(:, irk)
            do ib = 1, idx
               if (abs(onedvec(ib)) < 1.0E-6_dp .and. abs(onedvec2(ib)) < 1.0E-6_dp) cycle
               tmp = max(abs(onedvec(ib)), abs(onedvec2(ib)))
               if (abs(onedvec(ib) - onedvec2(ib)/ tmp) > 1.0E-2_dp) then
                  if (present(symmornot)) symmornot = .false.
                  write(*,*) 'symcheck_absvbs_ongrid: band ', ib, ' of irreducible k-point ', ik, &
                     ' is different from its equivalent k-point ', irk
                  write(*,*) 'symcheck_absvbs_ongrid: abs(1dvec(ib) - 1dvec2(ib)) = ', abs(onedvec(ib) - onedvec2(ib))
               endif
            enddo
         endif
      enddo
   enddo
   deallocate(onedvec, onedvec2)
end subroutine symcheck_absvbs_ongrid

subroutine symmetrize_absvbs_ongrid(kg, twodarray)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(inout) :: twodarray(:,:) ! assume (iband, ikpoint)

   integer :: irk, ik, idy

   !get the dimension of 2darray
   idy = size(twodarray, 2) ! number of reducible k-points

   !sanity check
   if (idy .ne. kg%nk) call errore('symmetrize_absvbs_ongrid','inconsistent dimension of 2darray',1)

   !loop over all reducible k-points, assign its irreducible k-point's value to it.
   do ik = 1, kg%nk
      irk = kg%kpt2ir(ik)
      twodarray(:, ik) = twodarray(:, irk)
   enddo

end subroutine symmetrize_absvbs_ongrid

!compute sum_{k} v(k)*v(k) for k corresponding to the same irreducible kpoints
subroutine sum_equivalent_kpoint(kg, vvprod, ndeg)
   use boltz_utils, only: velprod
   implicit none
   type(grid), intent(in) :: kg
   !vvprod(6, numb, kg%nk_irr), ndeg(kg%nk_irr): number of k-points
   real(dp), intent(out) :: vvprod(:,:,:), ndeg(:) 
   !
   integer :: ik, irk, ib, nd(3)
   real(dp) :: vel(3)
   
   nd = shape(vvprod)
   if(nd(1).ne.6 .or. nd(2).ne.kg%numb .or. nd(3).ne.kg%nk_irr .or. size(ndeg).ne.nd(3)) &
      call errore('sum_equivalent_kpoint','inconsistent demension of array vvprod',1)

   vvprod = 0.0_dp
   ndeg   = 0.0_dp
   do ik = 1, kg%nk
      irk = kg%kpt2ir(ik)
      ndeg(irk) = ndeg(irk) + 1.0_dp
      do ib = 1, kg%numb
         vel = kg%vnk(1:3, ib, ik)
         vvprod(1:6, ib, irk) = vvprod(1:6, ib, irk) + velprod(vel, vel)
      enddo
   enddo
end subroutine sum_equivalent_kpoint

subroutine boltz_grid_generate(kg, kdim, emin, emax, bmin, bmax, nsym, symop, g, el)
   use pert_const, only: dp
   use pert_param, only: prefix
   use yaml_utils, only: output_grid_yaml
   use pert_utils, only: find_free_unit 
   use boltz_utils, only: num2kpt, kpt2num, inside_win
   use band_structure, only: electron_wann
   use qe_mpi_mod,only: ionode, ionode_id, mp_split_pools, &
      mp_bcast, mp_barrier, mp_sum, inter_pool_comm ,stdout
   implicit none
   type(grid), intent(inout) :: kg
   ! energy range
   real(dp), intent(in) :: emin, emax
   ! dimension of the k-grid; band index range; symmetry operations.
   integer, intent(in) :: kdim(3), bmin, bmax, nsym, symop(3,3,nsym)
   character(len=1), intent(in) :: g
   type(electron_wann), intent(in), optional :: el
   ! local variables
   integer :: nkr, ik, num_irk, ns, i, kdx, ist, iend, ntet, nk
   integer, allocatable:: keq(:), irrk(:), irrk_e(:), irrk_t(:)
   real(dp) :: xkg(3), xkr(3)

   call start_clock('Grid_gen')
   nkr = kdim(1)*kdim(2)*kdim(3)
   ! nkr should be smaller than 2147483647, otherwise integer overflow occurs
   if (nkr <= 0) call errore('boltz_grid_generate', &
      'illegal kdim (or integer overflow)', 1)
   kg%ndim = kdim
   kg%kweight = 1.0_dp/dble(kdim(1))/dble(kdim(2))/dble(kdim(3))
   kg%tweight = kg%kweight/6.0_dp

   !keq might be a very large array, only allocate it on root process.
if (ionode) then
   !
   allocate( keq(0:nkr-1) )
   do ik = 0, nkr-1
      keq(ik) = ik
   enddo
   ! find irreducible k-points
   i = 0 !number of irreducible k-points
   do ik = 0, nkr-1
      !skip reducible k-points
      if( keq(ik) .ne. ik ) cycle
      i = i + 1
      !get crystal coordinate of k-points
      xkg(1:3) = num2kpt(ik, kdim)
      
      ! loop over point-group opterations, time-reversal is not included.
      do ns = 1, nsym
         xkr(1:3) = matmul(symop(1:3,1:3,ns), xkg(1:3))
         !map back to index coordinate.
         kdx = kpt2num(xkr, kdim)
         ! kdx < 0 : xkr is not in the grid
         ! from all the equivalent k-poins, choose the one with the smallest kdx
         if ( kdx > ik ) then
            keq( kdx ) = ik
         !sanity check
         elseif ( kdx >= 0 .and. kdx < ik) then
            call errore('boltz_grid_generate','Error in finding irreducible k', 1)
         endif
      enddo
   enddo
   num_irk = i
   !collect irreducible k-points and check consistence
   allocate( irrk(num_irk) )
   i = 0
   do ik = 0, nkr-1
      ! this is a irreducible kpoint
      if(keq(ik) .eq. ik) then 
         i = i + 1
         !from ik, we can get the crystal coord of this kpoints
         irrk(i) = ik
         !for irreducible k, point to its position in irrk.
         keq(ik) = i
      else
         !for reducible k, point to the position of their corresponidng irr-k
         ! the implicit assumption here is that, for reducible k, 'ik' must be
         ! larger than the 'ik' of their corresponing irreducbile kpoints.
         keq(ik) = keq( keq(ik) )
      endif
      ! Up to now, irrk(i) stores the i-th irreducible points.(in index coord)
      ! and, keq(ik) point the corresponding irreducible k in irrk.
   enddo
endif
   call mp_bcast(num_irk, ionode_id, inter_pool_comm)
   if(.not. allocated(irrk) ) allocate( irrk(num_irk) )
   call mp_bcast(irrk, ionode_id, inter_pool_comm)
   
   allocate( irrk_e(num_irk) )
   if( present(el) ) then
      irrk_e(:) = 0
      !this part is time consuming, split num_irk over pools (using mpi+openmp).
      call mp_split_pools(num_irk, ist, iend)
!! openmp parallel.
!$omp parallel do schedule(guided) default(shared) private(i, xkg)
      do i = ist, iend
         xkg(1:3) = num2kpt(irrk(i), kdim)
         ! check whether this irr-k has energy levels inside [emin, emax]
         ! only bands in [bmin, bmax] are counted.
         irrk_e(i) = inside_win(el, xkg, emin, emax, bmin, bmax)
         ! >= 0: no level inside [emin, emax], all are lower than emin or larger
         !  than emax. or some lower than emin and the others larger than emax.
         !   the value is the highest band that have energy lower than emin.
         ! = -2: have levels inside [emin, emax]
         ! = -4 : error, illegal result.
         if(irrk_e(i) .eq. -4) call &
            errore('boltz_grid_generate','Error in finding bands in [emin, emax]',1)
      enddo
!$omp end parallel do
      call mp_sum(irrk_e, inter_pool_comm)
   else
      !select all the kpoints
      irrk_e(:) = -2
   endif

if (ionode) then
   call count_tetra(kdim, keq, num_irk, irrk_e, ntet)
   kg%num_tet = ntet
   allocate( kg%tetra(2, 4, ntet) )
   ! irrk_t is used to reorder irreducible-k according to the order of 
   ! their appearance in tetrahedron. 0 if not selected.
   allocate( irrk_t(num_irk) )
   call collect_tetra(kdim, keq, num_irk, irrk_e, ntet, nk, irrk_t, kg%tetra)
   kg%nk = nk
   kg%nk_irr = count(irrk_t > 0)
   allocate(kg%irrk( kg%nk_irr ))
   allocate(kg%ir2kpt( kg%nk_irr ))
   !collect selected irreducible kpts 
   kg%irrk = -1; kg%ir2kpt = -1
   do i = 1, num_irk
      if(irrk_t(i) > 0) then
         ! irrk(i) : the index coordinate of i-th irreducible-kpoint
         ! irrk_t(i) : the new location of the i-th irreducible-kpoints
         kg%irrk( irrk_t(i) ) = irrk(i)
      endif
   enddo
   if(any( kg%irrk < 0 )) call errore('boltz_grid_generate', 'illegal irrk.', 1)
   deallocate(irrk_t, keq)
   
   !recover the full reducible grid 
   allocate( kg%kpt( kg%nk ), kg%kpt2ir( kg%nk ) )
   ! collect all the reducible k-points into kg%kpt, and their corresponding irreducible k.
   ! and also update thet tetra, so that
   !  tetra(1, ic, it) -> index of irreducible k-list: kg%irrk
   !  tetra(2, ic, it) -> index of   reducible k-list: kg%kpt
   call boltz_grid_setup(kg)

   !map the first points of each irreducible representaion in the selected reducible k-points to the irreducible k-points.
   !- syp
   do i = 1, kg%nk_irr
      do ik = 1, kg%nk
         if( kg%kpt2ir(ik) .eq. i ) then
            kg%ir2kpt(i) = ik
            exit
         endif
      enddo
   enddo

   !output to hdf5 file
   if( trim(g) .eq. 'k' ) then
       call output_grid(prefix, kg, emin, emax, bmin, bmax)
   elseif ( trim(g) .eq. 'q') then
       call output_qgrid(prefix, kg, g) 
   elseif ( trim(g) .eq. 'f') then
       call output_qgrid(prefix, kg, g) 
   endif
   !some message to stdout
   write(stdout,'(/5x,a,i11 )') 'number of tetrahedra selected:', kg%num_tet
   write(stdout,'( 5x,a,i11 )') 'number of reducible ' // g // '-points:', kg%nk
   write(stdout,'( 5x,a,i11/)') 'number of irreducible ' // g // '-points:', kg%nk_irr

   call output_grid_yaml(kg%num_tet, kg%nk, kg%nk_irr, g)

endif
   !release space
   deallocate(irrk, irrk_e)

   call mp_bcast(kg%num_tet, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%nk, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%nk_irr, ionode_id, inter_pool_comm)
   !
   if( .not. associated(kg%tetra) ) allocate( kg%tetra(2, 4, kg%num_tet) )
   if( .not. associated(kg%irrk) ) allocate( kg%irrk( kg%nk_irr ) )
   call mp_bcast(kg%tetra, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%irrk, ionode_id, inter_pool_comm)
   !
   if( .not. associated(kg%kpt) ) allocate( kg%kpt( kg%nk ) )
   if( .not. associated(kg%kpt2ir) ) allocate( kg%kpt2ir( kg%nk ) )
   if( .not. associated(kg%ir2kpt) ) allocate( kg%ir2kpt( kg%nk_irr ) )
   call mp_bcast(kg%kpt, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%kpt2ir, ionode_id, inter_pool_comm)
   call mp_bcast(kg%ir2kpt, ionode_id, inter_pool_comm)
   !
   call mp_barrier(inter_pool_comm)
   !
   call stop_clock('Grid_gen')
end subroutine boltz_grid_generate


subroutine count_tetra(kdim, kidx, num_irk, irrk_e, numt)
   implicit none
   !kgrid dimension: nk1 * nk2 * nk3
   integer, intent(in) :: kdim(3)
   ! if kidx(i) = i: the i-th kpts is irreducible kpt; 
   !      otherwise kidx(i) is the irreducible kpts equivalent to the i-th kpts
   integer, intent(in) :: kidx( 0:(kdim(1)*kdim(2)*kdim(3)-1) )
   !number of irreducible kpts
   integer, intent(in) :: num_irk
   !label the status of irreducible kpts:
   ! >= 0: no level fall into [emin, emax], the value is the highest band 
   !         that have energy lower than emin.
   ! = -2: have levels fall into [emin, emax]
   integer, intent(in) :: irrk_e(num_irk)
   !number of tetrahedra selected based on the energy of vertex
   integer, intent(out) :: numt
   !
   integer :: nt, nk1, nk2, nk3, i, j, k, corn(6,4), it, ic, ckin(4)
   
   nt = 0
   nk1 = kdim(1); nk2 = kdim(2); nk3 = kdim(3)
   !careful: i,j,k convention should be in consistent with num2kpt
   do i = 0, nk1-1
   do j = 0, nk2-1
   do k = 0, nk3-1
      call get_tetra(nk1, nk2, nk3, i, j, k, corn)
      do it = 1, 6
         do ic = 1, 4
            ! whether the corner have level inside emin, emax
            ckin(ic) = irrk_e( kidx( corn(it,ic) ) )
         enddo 
      ! select tetrahedra that contributes to dos in [emin, emax]. 
      ! if: 1. any vertex of the tetrahedron has energy fall into [emin, emax]
      !     2. eig(n,vertex) have energies cross the [emin, emax] (e.g. ene in 
      !     [emin, emax] can cut this tetrahedron, important for small interval). 
      ! then the tetrahedron and all its corners are selected.
         if(any(ckin(1:4) < -1) .or. any((ckin(1:4)-ckin(1)).ne.0)) then
            nt = nt + 1
         endif
      enddo 
   enddo; enddo; enddo
   numt = nt
end subroutine


subroutine collect_tetra(kdim, kidx, num_irk, irrk_e, numt, totk, irrk_t, tetra)
   implicit none
   integer, intent(in) :: kdim(3), num_irk, irrk_e(num_irk), numt
   integer, intent(inout) :: kidx( 0:(kdim(1)*kdim(2)*kdim(3)-1) ) 
   integer, intent(out) :: totk, irrk_t(num_irk), tetra(:,:,:)
   !
   integer :: i, j, k, it, ic, map_ir(4), irr_pk, all_pk
   integer :: corn(6,4), ckin(4), nt, nk1, nk2, nk3
   !init irrk_t
   irrk_t = 0
   tetra = 0

   nt = 0;  irr_pk = 0;  all_pk = 0
   nk1 = kdim(1); nk2 = kdim(2); nk3 = kdim(3)
   !careful: i,j,k convention should be in consistent with num2kpt
   do i = 0, nk1-1
   do j = 0, nk2-1
   do k = 0, nk3-1
      call get_tetra(nk1, nk2, nk3, i, j, k, corn)
      do it = 1, 6
         do ic = 1, 4
            ! map all k-list to irr_list, 1 ~ num_irk
            map_ir(ic) =  mod( kidx( corn(it,ic) ),  2*num_irk )
            ! whether the corner have level inside emin, emax
            ckin(ic) = irrk_e( map_ir(ic) )
         enddo 
         if(any(ckin(1:4) < -1) .or. any((ckin(1:4)-ckin(1)).ne.0)) then
            nt = nt + 1
            !sanity check
            if(nt > numt) call errore('collect_tetra','nt > numt', 1)
            do ic = 1, 4
               ! re-label the irreducible-k. skip the corner that already labeled.
               if( irrk_t( map_ir(ic) ) .eq. 0 ) then
                  irr_pk = irr_pk + 1
                  irrk_t( map_ir(ic) ) = irr_pk
               endif
               ! label the selected kpoints in full kgrid. includes reducible-k
               if( kidx( corn(it, ic) ) <= num_irk ) then
                  all_pk = all_pk + 1
                  ! if kidx(ik) > num_irk, then this k-points have already been counted.
                  kidx( corn(it, ic) ) = kidx( corn(it, ic) ) + 2*num_irk
               endif
            enddo
            !collect selected tetrahedron.
            do ic = 1, 4
               tetra(1,ic,nt) = irrk_t( map_ir(ic) )
               tetra(2,ic,nt) = corn(it,ic)
            enddo
         endif
      enddo 
      enddo; enddo; enddo
      totk = all_pk
end subroutine


subroutine get_tetra(nk1, nk2, nk3, i, j, k, tetra)
   implicit none
   integer, intent(in) :: nk1, nk2, nk3, i, j, k
   integer, intent(out) :: tetra(6,4)
   !
   integer :: ip1, jp1, kp1, n1, n2, n3, n4, n5, n6, n7, n8
   
   !careful: i,j,k convention should be in consistent with num2kpt
   ! n1-n8 are the indices of k-point 1-8 forming a cube
   ip1 = mod( i+1, nk1)
   jp1 = mod( j+1, nk2)
   kp1 = mod( k+1, nk3)
   n1 =   k +   j*nk3 +   i*nk2*nk3
   n2 =   k +   j*nk3 + ip1*nk2*nk3
   n3 =   k + jp1*nk3 +   i*nk2*nk3
   n4 =   k + jp1*nk3 + ip1*nk2*nk3
   n5 = kp1 +   j*nk3 +   i*nk2*nk3
   n6 = kp1 +   j*nk3 + ip1*nk2*nk3
   n7 = kp1 + jp1*nk3 +   i*nk2*nk3
   n8 = kp1 + jp1*nk3 + ip1*nk2*nk3
   !along 1-8, this is a better choice when i, j, k are symmetric
   ! 2, 4, 1, 8;
   tetra(1,1) = n2; tetra(1,2) = n4; tetra(1,3) = n1; tetra(1,4) = n8
   ! 4, 1, 3, 8;
   tetra(2,1) = n4; tetra(2,2) = n1; tetra(2,3) = n3; tetra(2,4) = n8
   ! 2, 1, 6, 8;
   tetra(3,1) = n2; tetra(3,2) = n1; tetra(3,3) = n6; tetra(3,4) = n8
   ! 1, 3, 8, 7;
   tetra(4,1) = n1; tetra(4,2) = n3; tetra(4,3) = n8; tetra(4,4) = n7
   ! 1, 8, 5, 7;
   tetra(5,1) = n1; tetra(5,2) = n8; tetra(5,3) = n5; tetra(5,4) = n7
   ! 1, 6, 8, 5;
   tetra(6,1) = n1; tetra(6,2) = n6; tetra(6,3) = n8; tetra(6,4) = n5
end subroutine get_tetra


subroutine output_grid(prefix, kg, emin, emax, bmin, bmax)
   use pert_const, only: dp
   use boltz_utils,only: num2kpt
   use pert_utils, only: find_free_unit
   use hdf5_utils
   implicit none
   character(len=80) :: prefix
   type(grid), intent(in) :: kg
   integer, intent(in) :: bmin, bmax
   real(dp), intent(in) :: emin, emax
   !
   integer :: iunit, ir, i, brange(2)
   real(dp) :: xkg(3), erange(2)
   integer(HID_T) :: file_id
   real(dp), allocatable :: kpts(:,:)
   
   !output irreducible k-list
   iunit = find_free_unit()
   open(iunit, file=trim(prefix)//'_tet.kpt',status='replace',form='formatted')
   write(iunit,'(1x,i8,a, 3(2x,i4), 2x,a,i8, 2x,a,i8)') kg%nk_irr, " crystal", &
      (kg%ndim(i), i=1,3), ' #.tetra ', kg%num_tet, ' #.tot.k ', kg%nk
   do ir = 1, kg%nk_irr
      xkg(1:3) = num2kpt( kg%irrk(ir), kg%ndim)
      write(iunit,'(3(f12.8,2x), f16.12)') xkg(1:3), 1.d0/real(kg%nk_irr, dp)
   enddo
   close(iunit)

   allocate( kpts(3, kg%nk) )
   do ir = 1, kg%nk
      kpts(:,ir) = num2kpt(kg%kpt(ir), kg%ndim)
   enddo
   
   erange = (/emin, emax/)
   brange = (/bmin, bmax/)
   !ouput tetra 
   call hdf_open_file(file_id, trim(prefix)//'_tet.h5', status='NEW')
   !
   call hdf_write_dataset(file_id, 'kgrid_dim', kg%ndim(1:3))
   call hdf_write_dataset(file_id, 'num_kpts',  kg%nk)
   call hdf_write_dataset(file_id, 'num_irr_kpts', kg%nk_irr)
   call hdf_write_dataset(file_id, 'num_tetra', kg%num_tet)
   call hdf_write_dataset(file_id, 'band_window', brange)
   call hdf_write_dataset(file_id, 'energy_window', erange)
   !
   call hdf_write_dataset(file_id, 'kpts_irr',  kg%irrk( 1:kg%nk_irr ) )
   call hdf_write_dataset(file_id, 'tetra', kg%tetra(1:2, 1:4, 1:kg%num_tet) )
   call hdf_write_dataset(file_id, 'kpts_all',  kg%kpt( 1:kg%nk ) )
   call hdf_write_dataset(file_id, 'kpt2ir', kg%kpt2ir(1:kg%nk ) )
   call hdf_write_dataset(file_id, 'ir2kpt', kg%kpt2ir(1:kg%nk_irr ) )
   
   call hdf_write_dataset(file_id, 'kpts_all_crys_coord', kpts)
   call hdf_close_file(file_id)

   deallocate( kpts )
end subroutine output_grid

! output q-grid without energy windows
subroutine output_qgrid(prefix, kg, g)
   use pert_const, only: dp
   use boltz_utils,only: num2kpt
   use pert_utils, only: find_free_unit
   use hdf5_utils
   implicit none
   character(len=80) :: prefix         !< prefix of file name
   type(grid), intent(in) :: kg        !< q-grid
   character(len=1), intent(in) :: g
   integer :: iunit                    !< free integer for output
   integer :: ir                       !< irreducible index
   integer :: i                        !< grid spatial dimension
   integer :: brange(2)
   real(dp) :: xkg(3)                  !< one irreducible q vector in lattice coordinates
   real(dp) :: erange(2)
   integer(HID_T) :: file_id           !< file ID
   real(dp), allocatable :: kpts(:,:)  !< all irreducible q vectors in lattice coordinates

   !output irreducible k-list
   iunit = find_free_unit()
   if( trim(g) .eq. 'q' ) then
      open(iunit, file=trim(prefix)//'_qg_tet.kpt',status='replace',form='formatted')
   elseif ( trim(g) .eq. 'f') then
      open(iunit, file=trim(prefix)//'_qgf_tet.kpt',status='replace',form='formatted')
   endif
   write(iunit,'(1x,i8,a, 3(2x,i4), 2x,a,i8, 2x,a,i8)') kg%nk_irr, " crystal", &
      (kg%ndim(i), i=1,3), ' #.tetra ', kg%num_tet, ' #.tot.k ', kg%nk
   do ir = 1, kg%nk_irr
      xkg(1:3) = num2kpt( kg%irrk(ir), kg%ndim)
      write(iunit,'(3(f12.8,2x), f16.12)') xkg(1:3), 1.d0/real(kg%nk_irr, dp)
   enddo
   close(iunit)

   !syp - for phonon velocity calculation input file
   iunit = find_free_unit()
   open(iunit, file=trim(prefix)//'_fullgrid.qpt',status='replace',form='formatted')
   write(iunit,'(1x,i8,a,i8,a, 3(2x,i4),a, 2x, i8,a,2x,i8,a)')kg%nk,"  #.FBZ.k", kg%nk_irr, " #.irr.k",(kg%ndim(i), i=1,3),&
      " crystal", kg%num_tet, ' #.tetra ', kg%nk, ' #.tot.k '
   do ir = 1, kg%nk
      xkg(1:3) = num2kpt( kg%kpt(ir), kg%ndim)
      write(iunit,'(3(f12.8,2x), 2x, f16.12)') xkg(1:3), 1.d0/real(kg%nk, dp)
   enddo
   close(iunit)
   !- syp
   allocate( kpts(3, kg%nk) )
   do ir = 1, kg%nk
      kpts(:,ir) = num2kpt(kg%kpt(ir), kg%ndim)
   enddo

   if( trim(g) .eq. 'q' ) then
      call hdf_open_file(file_id, trim(prefix)//'_qg_tet.h5', status='NEW')
   elseif ( trim(g) .eq. 'f') then
      call hdf_open_file(file_id, trim(prefix)//'_qgf_tet.h5', status='NEW')
   endif
   !
   call hdf_write_dataset(file_id, 'kgrid_dim', kg%ndim(1:3))
   call hdf_write_dataset(file_id, 'num_kpts',  kg%nk)
   call hdf_write_dataset(file_id, 'num_irr_kpts', kg%nk_irr)
   call hdf_write_dataset(file_id, 'num_tetra', kg%num_tet)
   call hdf_write_dataset(file_id, 'kpts_irr',  kg%irrk( 1:kg%nk_irr ) )
   call hdf_write_dataset(file_id, 'tetra', kg%tetra(1:2, 1:4, 1:kg%num_tet) )
   call hdf_write_dataset(file_id, 'kpts_all',  kg%kpt( 1:kg%nk ) )
   call hdf_write_dataset(file_id, 'kpt2ir', kg%kpt2ir(1:kg%nk ) )
   call hdf_write_dataset(file_id, 'ir2kpt', kg%ir2kpt(1:kg%nk_irr ) )
   call hdf_write_dataset(file_id, 'kpts_all_crys_coord', kpts)
   call hdf_close_file(file_id)

   deallocate( kpts )
end subroutine output_qgrid


! Read input files: tet.kpt, tet.tetra
subroutine boltz_grid_load(kg, g)
   use pert_const, only: dp
   use yaml_utils, only: output_grid_yaml
   use pert_utils, only: find_free_unit
   use boltz_utils, only: kpt2num
   use pert_param, only: prefix, boltz_kdim, boltz_qdim
   use qe_mpi_mod, only: stdout, ionode, ionode_id, mp_bcast, inter_pool_comm
   use hdf5_utils
   implicit none
   type(grid), intent(inout) :: kg
   character(len=1), intent(in) :: g ! grid type
   ! local variables
   integer :: kunit, tunit, ios, ik, i, irec, ikpt !, i, ib, ierr, sunit
   integer :: ir_nk, itmp(3), num_t, tot_k, ttmp(8), nkr, tmp !, nept, 
   character(len=120) :: fname, fname1, ctmp1, ctmp2, ctmp3
   real(dp) :: ktmp(3) !, wtmp, rtmp
   logical :: has_file
   integer(HID_T) :: file_id
   integer :: boltz_grid_dim(3) ! general dimension
   
   !sanity check
   if( associated(kg%tetra) .and. associated(kg%irrk) ) then
      write(stdout, '(5x,a)') "Warn: kgrid is already loaded, skip the duplicate loading."
      return
   endif

   if (ionode) then
      write(stdout,'(5x,a)') ">loading kpt, tetra from " // trim(prefix) // '_tet.h5'
      !
      select case (g)
         case ('k')
            fname = trim(prefix)//"_tet.h5"
            fname1 = trim(prefix)//"_tet.kpt"
            boltz_grid_dim = boltz_kdim
         case ('q')
            fname = trim(prefix)//"_qg_tet.h5"
            fname1 = trim(prefix)//"_qg_tet.kpt"
            boltz_grid_dim = boltz_qdim
         case ('f')
                 fname = trim(prefix)//"_qgf_tet.h5"
                 fname1 = trim(prefix)//"_qgf_tet.kpt"
                 boltz_grid_dim = boltz_kdim
         case default
            call errore('boltz_grid_load','illegal grid name', 1)
      end select

      inquire(file=trim(fname), exist=has_file)
      if(.not. has_file) call errore('boltz_grid_load','missing '// trim(fname), 1)
      !
      call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
      call hdf_read_dataset(file_id, 'kgrid_dim', itmp(1:3))
      call hdf_read_dataset(file_id, 'num_kpts',  tot_k)
      call hdf_read_dataset(file_id, 'num_irr_kpts', ir_nk)
      call hdf_read_dataset(file_id, 'num_tetra', num_t)
   
      if(any( (boltz_grid_dim(1:3)-itmp(1:3)) .ne. 0) ) &
         call errore('boltz_grid_load','boltz_kdim does not match '// trim(fname),1)
      
      kg%nk  = tot_k
      kg%nk_irr  = ir_nk
      kg%num_tet = num_t
      kg%ndim(1:3) = boltz_grid_dim(1:3)
      allocate( kg%irrk( kg%nk_irr ) )
      allocate( kg%tetra(2, 4, kg%num_tet) )
      allocate( kg%kpt( kg%nk ) )
      allocate( kg%kpt2ir( kg%nk ) )
      allocate( kg%ir2kpt( kg%nk_irr ) )
      !
      call hdf_read_dataset(file_id, 'kpts_irr', kg%irrk)
      call hdf_read_dataset(file_id, 'tetra', kg%tetra)
      call hdf_read_dataset(file_id, 'kpts_all', kg%kpt)
      call hdf_read_dataset(file_id, 'kpt2ir', kg%kpt2ir)
      call hdf_read_dataset(file_id, 'ir2kpt', kg%ir2kpt)
      call hdf_close_file(file_id)
      
      !if prefix_tet.kpt exist, check consistency
      inquire(file=trim(fname1), exist=has_file)
      if(has_file) then
         !read irreducible k-points  
         kunit = find_free_unit()
         open(kunit,file=trim(fname1),status='old',form='formatted',err=100,iostat=ios)
         read(kunit, *) ir_nk, ctmp1, itmp(1:3), ctmp2, num_t, ctmp3, tot_k
      
         if((kg%nk_irr .ne. ir_nk).or.(kg%num_tet .ne. num_t).or.(kg%nk .ne. tot_k) &
            .or. any(boltz_grid_dim-itmp .ne. 0)) &
         call errore('boltz_grid_load', 'inconsistent ' // trim(fname) // ', ' // trim(fname1), 1)
      
         do ik = 1, kg%nk_irr
            read(kunit,*) ktmp(1:3)
            ikpt = kpt2num( ktmp(1:3), boltz_grid_dim )
            if( kg%irrk(ik) .ne. ikpt ) call errore('boltz_grid_load', &
               'inconsistent kpoints in ' // trim(fname) // ' and ' // trim(fname1), 1)
         enddo
         
         close(kunit)
      endif
      
      !some message to stdout
      write(stdout,'(/5x,a,i11 )') 'number of tetrahedra selected:', kg%num_tet
      write(stdout,'( 5x,a,i11 )') 'number of reducible ' // g // '-points:', kg%nk
      write(stdout,'( 5x,a,i11/)') 'number of irreducible ' // g // '-points:', kg%nk_irr

      call output_grid_yaml(kg%num_tet, kg%nk, kg%nk_irr, g)

   endif


   call mp_bcast( kg%nk,  ionode_id, inter_pool_comm)
   call mp_bcast( kg%nk_irr,  ionode_id, inter_pool_comm)
   call mp_bcast( kg%num_tet, ionode_id, inter_pool_comm)
   call mp_bcast( kg%ndim, ionode_id, inter_pool_comm)
   !
   if(.not. associated(kg%irrk) )  allocate( kg%irrk( kg%nk_irr ) )
   if(.not. associated(kg%tetra))  allocate( kg%tetra(2, 4, kg%num_tet) )
   call mp_bcast(kg%irrk, ionode_id, inter_pool_comm)
   call mp_bcast(kg%tetra, ionode_id, inter_pool_comm)
   !
   if( .not. associated(kg%kpt) ) allocate( kg%kpt( kg%nk ) )
   if( .not. associated(kg%kpt2ir) ) allocate( kg%kpt2ir( kg%nk ) )
   if( .not. associated(kg%ir2kpt) ) allocate( kg%ir2kpt( kg%nk_irr ) )
   call mp_bcast(kg%kpt, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%kpt2ir, ionode_id, inter_pool_comm)
   call mp_bcast(kg%ir2kpt, ionode_id, inter_pool_comm)

   !weight of the k-points
   kg%kweight = 1.0_dp / dble(kg%ndim(1)) / dble(kg%ndim(2)) / dble(kg%ndim(3))
   kg%tweight = kg%kweight / 6.0_dp
   
   return
100 call errore('boltz_grid_load','opening file '//trim(fname1),abs(ios))
end subroutine boltz_grid_load


!! read irreducible k-list and tetra, and recover the full grid
subroutine boltz_grid_setup(kg)
   implicit none
   type(grid), intent(inout) :: kg
   ! local variable
   integer, parameter :: nsym = 48 ! max number of symmetry operations
   integer :: it, ic, ir, n, i, ik, ierr
   integer, allocatable :: fktmp(:,:), kcol(:)
   
   !initial work array
   allocate( fktmp(nsym, kg%nk_irr), kcol(kg%nk),  stat=ierr)
   if(ierr /= 0) call errore('boltz_grid_setup','allocating space failed',1)

   fktmp(:,:) = -1
   ! collect kpoint corresponding to the same irreducible k into fktmp(:,ir)
   do it = 1, kg%num_tet
   do ic = 1, 4
      ! ir: index of irreducible-k list; n: kpt coordinate, num2kpt(n, kdim)
      ir = kg%tetra(1, ic, it)
      n  = kg%tetra(2, ic ,it)
      ! scan all the element of fktmp(:,ir)
      do i = 1, nsym
         ! fktmp(:,ir) store all the possible kpoints (n) appears in tetra
         ! that equivalent to irreducible k-list ir.
         ! fktmp(i,ir) == 0: means this position is empty.
         ! scan all the non-empty elements in fktmp(:,ir), if found the same n, 
         ! which means the k-points already appears before, then move to next. 
         if( n .eq. fktmp(i, ir) ) then
            kg%tetra(2, ic, it) = i
            exit
         ! if kpoint n never appears before, add to empty position of fktmp(:,ir)
         elseif(fktmp(i,ir) < 0) then
            fktmp(i,ir) = n
            kg%tetra(2,ic,it) = i
            exit
         ! if fktmp(i,ir) .ne. -1 and .ne. n, then continue.
         ! but after i reach nsym, fktmp still .ne. -1 and .ne. n, raise error.
         elseif(i .eq. nsym) then
            call errore('boltz_grid_setup','error in setup full k-grid', 1)
         endif
      enddo
   enddo; enddo
   
   !check and collect all the kpt in kcol
   kcol(:) = -1;  ik = 0
   do ir = 1, kg%nk_irr
   do i = 1, nsym
      if(fktmp(i,ir) < 0) cycle
      ik = ik + 1
      ! sanity check
      if(ik > kg%nk) call errore('boltz_grid_setup',"ik > kg%nk",1)
      kcol(ik) = fktmp(i, ir)
   enddo; enddo
   ! sanity check
   if(ik .ne. kg%nk) call errore('boltz_grid_setup', "ik .ne. kg%nk.",1)
   
   !sort and store in kg%kpt
   !initialize kg%kpt with kpoints stored in fktmp in ascending order
   kg%kpt2ir(:) = 0;  kg%kpt(:) = -1
   do ir = 1, kg%nk_irr
   do i = 1, nsym
      if(fktmp(i,ir) < 0) cycle
      ! find position of this kpt in the ordered kgrid.
      ic = 1
      do ik = 1, kg%nk
         if( fktmp(i,ir) > kcol(ik) ) ic = ic + 1
      enddo
      ! if this kpt is larger then n kpts, its position is n+1.
      kg%kpt( ic ) = fktmp(i, ir)
      ! this k-points equivalent to irreducible kpoints-ir
      kg%kpt2ir( ic ) = ir
      ! reassign fktmp(i,ir) to its position in kgrid
      fktmp(i,ir) = ic
   enddo; enddo
   !debug & test
   if( any(kg%kpt(:) < 0) .or. any(kg%kpt2ir(:) .eq. 0) ) &
      call errore('boltz_grid_setup','kgrid sorting failed.',1)

   ! update tetra, afterward:
   !  tetra(1, ic, it) -> index of irreducible k-list: kg%irrk
   !  tetra(2, ic, it) -> index of   reducible k-list: kg%kpt
   do it = 1, kg%num_tet
   do ic = 1, 4
      ir = kg%tetra(1, ic, it)
      i  = kg%tetra(2, ic ,it)
      kg%tetra(2, ic, it) = fktmp(i, ir)
   enddo; enddo
   !release memory
   deallocate(fktmp, kcol)
end subroutine boltz_grid_setup

end module boltz_grid
