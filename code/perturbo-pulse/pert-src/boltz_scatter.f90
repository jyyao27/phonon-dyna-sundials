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

#include "cdyn_config.h"

module boltz_scatter
   use kinds, only: i8b
   use boltz_grid, only: grid
   use boltz_utils, only: idx2band
   use yaml_utils, only: ymlout, output_scatter_info_yaml
   use pert_const, only: dp, czero
   use pert_data,  only: epr_fid, qc_dim, kc_dim
   use pert_param, only: phfreq_cutoff, phfreq_cutoff_ph, delta_smear, delta_smear_ph, use_mem, &
      load_scatter_eph, load_scatter_ph3, dyna_phph, scat_impl
   use band_structure, only: electron_wann, solve_eigenvalue_vector
   use phonon_dispersion, only: lattice_ifc, lattice_ifct, solve_phonon_modes
   use parallel_include
   use qe_mpi_mod, only: stdout, ionode, mp_barrier, inter_pool_comm, npool, my_pool_id, mp_sum, &
      mp_split_pools
   use elphon_coupling_matrix, only: elph_mat_wann, release_elph_mat_wann, &
      init_elph_mat_wann, eph_fourier_el_para, eph_fourier_elph, eph_transform_fast
   use boltz_scatter_sizest
#if defined(STORE_G2DT)
   use pert_const, only: pi
   use pert_utils, only: gauss
#endif
   use force_constant, only: lattice_ifc

   implicit none
   public

   type, private :: channel_ph
      integer :: ik
      integer :: nkq 
      integer, pointer :: iq(:) => null()
      integer, pointer :: nchl(:) => null()
   end type
   !=========================================================================
   ! Types for kq_grid data structure

   type, private :: kq_channel
      integer :: ik
      integer :: nkq
      !
      integer, pointer :: iq(:) => null()
      integer, pointer :: ikq(:) => null()
      integer, pointer :: nchl(:) => null()
   end type

   type :: phonon_grid
      integer :: ndim(3)  !size of the grid
      !
      integer :: npts  !number of points selected
      integer, pointer :: list(:) => null()
      real(dp), pointer :: freq(:,:) => null()
      !
      real(dp):: weight  ! 1/(ndim(1) * ndim(2) * ndim(3))
   end type

   type :: t_scatt_compact
      integer :: scat_index
      integer :: bands_index
      real(dp) :: g2
   end type
   !=========================================================================
   ! Types for hybrid scattering data structures

   !user-defined data structure for scattering channel. based on kgrid
   type :: t_scatter
      !location of k and kq in kgrid%kpt: xk(1:3) = num2kpt(kgrid%kpt(ik), ndim)
      !NOTE: (k, kq) and (kq, k) are equivalent, only one of them are stored.
      integer :: ik
      integer :: ikq
      !location of q = (k+q) - k in qgrid%list: xq(1:3) = num2kpt(qgrid%list(iq), ndim)
      integer :: iq

#if defined(SCAT_FWD)
      ! It's possible to eliminate one of these two integer values, and just
      ! store the end-index of each kq-pair's sequence of channels.  But,
      ! this approach is just simpler to parallelize with OpenMP/OpenACC.

      integer(kind=i8b) :: chl_start_index

      !number of (mkq, nk, mu) terms, mu is phonon modes index
      !only g_m,n,mu of (k, kq) satisfy the following condition are stored.
      !abs( abs(e_nk - e_mkq) - w_qmu ) < delta_E. (delta_E = e_thr*broadening)
      integer :: nchl
#endif
   end type t_scatter

   ! A "flattened" and denormalized scattering-channel array for rapid,
   ! (hopefully) OpenMP/OpenACC-friendly traversal.  The main difference
   ! between t_scatt and t_flatscat is that t_flatscat doesn't have any
   ! sub-arrays to access.  The downside is that the (ik, ikq, iq) values
   ! are repeated across multiple rows.
   type :: t_scatter_channel
#if defined(SCAT_REV)
      !location of k and kq in kgrid%kpt: xk(1:3) = num2kpt(kgrid%kpt(ik), ndim)
      !NOTE: (k, kq) and (kq, k) are equivalent, only one of them are stored.
      integer :: scat_index
#endif

      integer :: bands_index   ! (m, n, im) packed into a single integer

      !sypeng: to be compatible with Kelly's additional treatment of dt for ph-e coupling
      !sypeng: TODO
!#if !defined(STORE_G2DT)
      real(dp) :: eph_g2       ! g2
!#else
      real(dp) :: g2_dt1       ! g2*dt1
      real(dp) :: g2_dt2       ! g2*dt2
!#endif
   end type t_scatter_channel

   !=========================================================================
   ! Data structures for dynamics simulation

   integer, private, save :: nk_loc, nq_loc, nq_loc_fine, iq_st1, iq_end1
   type(kq_channel), allocatable, private, save :: kq_pair(:)
   type(channel_ph), allocatable, private, save :: qq_pair(:)
   integer, allocatable, save :: qq_pair_arr(:,:)
   integer, allocatable, save :: q_start_list(:)
   !
   type(phonon_grid), save :: qgrid
   !
   integer, save :: num_kq_pairs         ! replaces num_scat
   integer, save :: num_qq_pair
   
   integer(kind=i8b), save :: num_scatter_channels ! replaces num_flatscat, the same as num_scat in ph-ph
   integer(kind=i8b), save :: ph3_num_scat ! replaces num_flatscat
   !integer(kind=8), save :: ph3_num_scat ! replaces num_flatscat
   integer, save :: num_qq_pair_total ! replaces num_flatscat
   !
   integer, allocatable, save :: scatter_ph(:,:)
   !
   type(t_scatter), allocatable, save :: scatter(:)
   type(t_scatter_channel), allocatable, save :: scatter_channels(:)
   type(t_scatt_compact), allocatable, save :: tscat_ph(:)
   real(dp), allocatable, save :: tsmear_ph_list(:)
   real(dp), allocatable, save :: tsmear_list(:)

contains


subroutine boltz_scatter_setup(kg, el, ph, qdim)
   implicit none
   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   !dimension for the q-grid. q-should be commensurate with the k-grid
   integer, intent(in) :: qdim(3)
   !local
   integer :: nkpt, i, ik, collect_pools(npool), nq_col(npool)
   integer :: prev_kqp
   integer(kind=i8b) :: prev_sch
   logical :: overflow
   real(dp) :: wcut, e_thr
   type(elph_mat_wann) :: elph

   !init
   nkpt = kg%nk
   wcut  = phfreq_cutoff
   e_thr = delta_smear*3.0_dp  !exp(-3*3) ~ 10^-4

   !interleave distribution of kpoints over pools, for better load balance
   nk_loc = nkpt/npool + merge(1, 0, my_pool_id < mod(nkpt, npool))
   allocate( kq_pair(nk_loc) )
   !store the allocated kpts (location in kg%kpt) of the current pool
   i = 0
   do ik = my_pool_id+1, nkpt, npool
      i = i + 1
      kq_pair(i)%ik = ik
   enddo
   !debug
   if(i .ne. nk_loc) call errore('boltz_scatter_setup', 'i .ne. nk_loc', 1)
   ! only compute qgrid if not dyna_phph
   qgrid%ndim(1:3) = qdim(1:3)
   qgrid%weight = 1.0_dp / dble(qdim(1)) / dble(qdim(2)) / dble(qdim(3))
   !generate (k,q) pair (initialized kq_pair, qgrid)
   call setup_kq_pair(kg, el, ph, wcut, e_thr, nk_loc, kq_pair, qgrid)

   !
   ! compute and save g2 to disk
   if(.not. load_scatter_eph) then
      call init_elph_mat_wann(epr_fid, kc_dim, qc_dim, elph)
      !compute g2 and sove them to temperary files
      call compute_scatter_eph_g2 &
         (kg, el, ph, elph, qgrid, wcut, e_thr, use_mem)
      !release space
      call release_elph_mat_wann(elph)
   endif

   !!! END KQ-GRID SETUP !!!

   !
   ! Compute the total number of (k,q) pairs and the total number of scatter
   ! channels on this MPI node.
   num_kq_pairs = 0
   num_scatter_channels = 0
   overflow = .false.
   do i = 1, nk_loc
      if (kq_pair(i)%nkq .eq. 0) cycle

      ! For detecting overflows
      prev_kqp = num_kq_pairs
      prev_sch = num_scatter_channels

      num_kq_pairs = num_kq_pairs + kq_pair(i)%nkq
      num_scatter_channels = num_scatter_channels + sum( kq_pair(i)%nchl(:) )

      if (prev_kqp >= 0 .and. num_kq_pairs < 0) then
         write(stdout, *) 'ERROR:  Detected overflow on num_kq_pairs!'
         flush(stdout)
         overflow = .true.
      endif

      if (prev_sch >= 0 .and. num_scatter_channels < 0) then
         write(stdout, *) 'ERROR:  Detected overflow on num_scatter_channels!'
         flush(stdout)
         overflow = .true.
      endif
   enddo

   if (overflow) call errore('setup_scatter', 'Overflows detected!', 1)

   !
   ! Collect the number of q-points and (k,q)-pairs in each pool so we can
   ! print it and save it to the YAML file.
   nq_col = 0
   if (dyna_phph) then
      nq_col( my_pool_id+1 ) = qdim(1)*qdim(2)*qdim(3)
   else
      nq_col( my_pool_id+1 ) = qgrid%npts
   endif
   call mp_sum( nq_col, inter_pool_comm )
   !
   collect_pools = 0
   collect_pools( my_pool_id+1 ) = num_kq_pairs
   call mp_sum( collect_pools, inter_pool_comm )
   !output info
   write(stdout, '(5x, a)') "No. of q-points and (k,q) pairs on each processor:"
   do i = 1, npool
      write(stdout, '(5x, a, i4, a, i9, a, i12)')  &
         "Proc.", i, ':  #.q ', nq_col(i), ' |  #.(k,q)', collect_pools(i)
   enddo

   write(stdout, '(/, 5x, a, i23)') "Total No. of q-points:", sum(nq_col)
   write(stdout, '(5x, a, i20, /)') "Total No. of (k,q) pairs:", sum(INT8(collect_pools))
   flush(stdout)

   ! output to the YAML file
   call output_scatter_info_yaml(npool, nq_col, collect_pools)

   call setup_scatter(kg, el, ph%nm)

   if (ionode) then
      ! Log data-structure allocation details into the YAML output file.
      write(ymlout, '(/,a)') 'boltz_scatter memory:'
      call output_kq_pair_alloc_stats_to_yaml()
      call output_scatter_alloc_stats_to_yaml()
      call output_scatter_channels_alloc_stats_to_yaml()
   endif

   do i = 1, nk_loc
      if( kq_pair(i)%nkq == 0 ) cycle
      deallocate( kq_pair(i)%iq, kq_pair(i)%ikq, kq_pair(i)%nchl )
   enddo
   deallocate( kq_pair )
   !syp - 20241007: commendted by sypeng, this is needed by rates_scat_int
   !      which will be called in coupled method when dyna_phph is true.
!  if ( dyna_phph ) then
!     ! qgrid bascially not needed, deallocate to save space
!     deallocate(qgrid%list)
!     deallocate(qgrid%freq)
!  endif

end subroutine boltz_scatter_setup


! Generate the initial kq-pair data structure that drives the rest of the
! dynamics simulation.
subroutine setup_kq_pair(kg, el, ph, wcut, e_thr, nk_loc, kq_pair, qg)
   use pert_utils, only: match_table
   use pert_param, only: dyna_phph
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus
   implicit none
   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   real(dp), intent(in) :: wcut, e_thr
   !
   integer, intent(in) :: nk_loc
   type(kq_channel), intent(inout) :: kq_pair(:)
   !
   type(phonon_grid), intent(inout) :: qg
   !
   !local variables
   integer :: nkpt, nmod, maxq, nkq, nq_tmp, i, ik, ikq, j, qidx, iq, k
   logical :: ltable(kg%numb, kg%numb, ph%nm)
   integer, allocatable :: qcollect(:)
   real(dp) :: xq(3)
   !
   integer,  pointer :: ql_tmp(:), q_check(:)
   real(dp), pointer :: qw_tmp(:,:)
   type(kq_channel), pointer :: kq_tmp(:)

   call start_clock('set_kq_pair')
   !init
   nkpt = kg%nk
   nmod = ph%nm

   allocate( kq_tmp(nk_loc) )
   !find iq
   !maxq = qg%ndim(1)*qg%ndim(2)*qg%ndim(3) - 1
   if ( dyna_phph) then
      ! consider all possible q's 
      maxq = qg%ndim(1)*qg%ndim(2)*qg%ndim(3) - 1
   else
      maxq = qg%ndim(3) * (1 +  qg%ndim(2) * (1 + qg%ndim(1)) ) / 2
   endif
   allocate(qcollect(0:maxq))
   qcollect = 0
   !(k, k') and (k', k) are connected, only one of them is stored.
   !to balance load, each kpoint is associated with equal number of k+q.
   !so for each ik in grid%kpt, the (k+q) are choose as ik+1, ik+2, ...., ik+(nkpt-1)/2
   !  if (ik + i) is larger than nkpt, restart from 1. for example, nkpt + 5 -> 5
   !
   !initialize and pre-allocate space
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      kq_tmp(i)%ik = ik
      !assign (k+q) evenly for each k, for load balance
      nkq = (nkpt-1)/2 + merge(1, 0, ((mod(nkpt-1,2)>0) .and. ik<=nkpt/2))
      kq_tmp(i)%nkq = nkq

      allocate( kq_tmp(i)%iq(nkq) )
      kq_tmp(i)%iq(:) = -1
   enddo
   !
   !collect q-points and store them in qg%list.
!$omp parallel do schedule(guided) default(shared) private(i, ik, nkq, j, ikq, qidx)
   do i = 1, nk_loc
      ik = kq_tmp(i)%ik
      nkq = kq_tmp(i)%nkq

      !loop over (k+q) points assigned to this current k-points
      do j = 1, nkq
         ikq = mod(ik+j-1, nkpt) + 1
         !store (k+q) - k in kq_pair
         !(since q and -q has the same freqencies, we only store one of them)
         ! if xq = (ikq - iq) is not on the q-grid (qg), then qidx = -1.
         !qidx = min( kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim), &
         if ( dyna_phph ) then
            qidx = kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim)
         else 
            !syp - for e-ph, they address both absorption and emission simultaneously
            !    - the first kpts_minus is for absorption, the second is for emission
            !syp - 20241011 [TOCHECK] I suggest the group to revisit here, together with the 
            !      comment in boltz_utiles.f90/num2kpt
            qidx = min( kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim), &
               kpts_minus(kg%kpt(ik), kg%kpt(ikq), kg%ndim, qg%ndim) )
         endif
         !debug
         !if(qidx > maxq) call errore('setup_kq_pair', 'qidx > maxq', 1)
         !
         !if this xq is on the q-grid qg.
         if(qidx > -1) then
            kq_tmp(i)%iq(j)  = qidx
!$omp atomic
            qcollect(qidx) = qcollect(qidx) + 1
         endif
      enddo
   enddo
!$omp end parallel do

   !collect q-points
   iq = 0
   do qidx = 0, maxq
      if(qcollect(qidx) > 0) then
         iq = iq + 1
         !syp - qcollect(qidx) store the serial number of each serial q-point (kind of full list)
         !    - in the selected q-point list.
         !    - for exam. iq means qidx-th q-point is No. iq in the selected q-point list.
         qcollect(qidx) = iq
      endif
   enddo
   nq_tmp = iq

   allocate( ql_tmp(nq_tmp) )
   do qidx = 0, maxq
      iq = qcollect(qidx)
      !syp - inverse of qcollect, but similar with kq_tmp(i)%iq(j) above.
      !    - store the iq-th selected q-point's representation in the full q-point list. 
      if( iq > 0) ql_tmp( iq ) = qidx
   enddo

   !update kq_pair%iq
   do i = 1, nk_loc
      do j = 1, kq_tmp(i)%nkq
         !syp - qidx is j-th iq's representation in full q-point list.
         qidx = kq_tmp(i)%iq(j)
         !qidx = -1 if xq is not on the q-grid.
         if(qidx > -1) then
            iq = qcollect(qidx)
            !syp - Now change kq_tmp(i)%iq(j)'s value from the representation in full q-point list
            !    - to the representation in the selected q-point list.
            kq_tmp(i)%iq(j) = iq
         endif
      enddo
   enddo
   !release space
   deallocate(qcollect)

   !allocate space
   allocate(qw_tmp(nmod, nq_tmp))
   do i = 1, nk_loc
      nkq = kq_tmp(i)%nkq
      allocate( kq_tmp(i)%ikq(nkq)  )
      allocate( kq_tmp(i)%nchl(nkq) )
      !initialize
      kq_tmp(i)%nchl(:) = 0
   enddo

   !compute phonon freq for qpts
!$omp parallel default(shared) private(iq, xq, i, ik, j, ikq, ltable)
!$omp do schedule(guided)
   do iq = 1, nq_tmp
      xq  = num2kpt( ql_tmp(iq),  qg%ndim )
      call solve_phonon_modes(ph, xq, qw_tmp(:,iq))
   enddo
!$omp end do

!$omp do schedule(guided)
   do i = 1, nk_loc
      ik = kq_tmp(i)%ik
      !nkq = (nkpt-1)/2 + merge(1, 0, ((mod(nkpt-1,2)>0) .and. ik<=nkpt/2))
      do j = 1, kq_tmp(i)%nkq
         ikq = mod(ik+j-1, nkpt) + 1
         kq_tmp(i)%ikq(j) = ikq
         !
         iq = kq_tmp(i)%iq(j)
         if(iq > 0) then
            call match_table(kg%enk(:,ikq), kg%enk(:,ik), qw_tmp(:,iq), wcut, e_thr, ltable)
            kq_tmp(i)%nchl(j) = count(ltable)
         endif
      enddo
   enddo
!$omp end do
!$omp end parallel

   !update qg%list
   allocate( q_check(nq_tmp) )
   q_check = 0
   !collect (k,q) pair that has at least one valid scattering channel (nchl > 0).
   do i= 1, nk_loc
      nkq = count( kq_tmp(i)%nchl(:) > 0 )
      kq_pair(i)%nkq = nkq
      ! only needed when nkq > 0
      if( nkq > 0) then !!TODO kellyyao - added this back
         allocate( kq_pair(i)%iq(nkq), kq_pair(i)%ikq(nkq), kq_pair(i)%nchl(nkq) )
         !
         k = 0
         do j = 1, kq_tmp(i)%nkq
            if( kq_tmp(i)%nchl(j) > 0 ) then
               k = k + 1
               !
               kq_pair(i)%ikq(k) = kq_tmp(i)%ikq(j)
               kq_pair(i)%nchl(k) = kq_tmp(i)%nchl(j)
               iq = kq_tmp(i)%iq(j)
               kq_pair(i)%iq(k) = iq
               !
               q_check(iq) = q_check(iq) + 1
            endif
         enddo
         !debug
         if(k .ne. nkq) call errore('setup_kq_pair', 'k .ne. nkq', 1)
      endif
      !
      deallocate(kq_tmp(i)%iq, kq_tmp(i)%ikq, kq_tmp(i)%nchl)
   enddo
   deallocate( kq_tmp )

   !if number of q-points is reduced.
   if( count(q_check > 0) .ne. nq_tmp ) then
      k = 0
      do i = 1, nq_tmp
         if(q_check(i) > 0) then
            k = k + 1
            q_check(i) = k
         endif
      enddo
      qg%npts = k

      !update qg%list
      allocate( qg%list( qg%npts ), qg%freq(nmod, qg%npts) )
      do i = 1, nq_tmp
         if(q_check(i) > 0) then
            !syp - qg%list(i) means the representation of i-th selected q-point in the 
            !    - whole q-point list.
            qg%list( q_check(i) ) = ql_tmp(i)
            qg%freq(:, q_check(i)) = qw_tmp(:, i)
         endif
      enddo
      deallocate( ql_tmp, qw_tmp )

      !update kq_pair
      do i = 1, nk_loc
         do j = 1, kq_pair(i)%nkq
            iq = kq_pair(i)%iq(j)
            !syp - kq_pair(i)%iq(j) store the place of j-th q-point in the
            !    - whole selected q-point list, not the whole list.
            kq_pair(i)%iq(j) = q_check(iq)
         enddo
      enddo
   else
      qg%npts = nq_tmp
      qg%list => ql_tmp
      qg%freq => qw_tmp
   endif
   deallocate( q_check )

   call stop_clock('set_kq_pair')
end subroutine setup_kq_pair


! Generate an HDF5 file for the current MPI pool (i.e. MPI node), and generate
! the eph and g2 values into the HDF5 file.
!
! Each MPI node will generate its own HDF5 file, so the working directory will
! end up with N files for N nodes.  These files can then be reloaded on
! subsequent perturbo runs, as long as the number of MPI nodes doesn't change.
subroutine compute_scatter_eph_g2(kg, el, ph, elph, qg, wcut, e_thr, use_mem)
   use pert_param, only: prefix, tmp_dir, dyna_phph, adapt_smear_eph
   use pert_utils, only: abs2, match_table
   use polar_correction,  only: eph_wan_longrange
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus, &
      calc_adaptive_smear
   use pert_output, only: progressbar_init, progressbar
   use hdf5_utils
   implicit none
   type(grid), intent(in) :: kg
   type(phonon_grid), intent(in) :: qg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   type(elph_mat_wann), intent(in) :: elph
   !
   logical, intent(in) :: use_mem
   real(dp), intent(in) :: wcut, e_thr
   !
   integer :: nmod, qidx, numq
   integer :: i, j, ik, iq, ikq, it, im, n, m, nkq, tot_chnl, bands(3)
   real(dp) :: xq(3), xk(3), xkq(3), wq(ph%nm), ek(el%nb), ekq(el%nb)
   complex(dp) :: uqt(ph%nm, ph%nm), uk(el%nb, el%nb), ukq(el%nb, el%nb), gpol(ph%nm)
   logical :: ltable(kg%numb, kg%numb, ph%nm)
   !
   integer, allocatable :: ist(:), bnd_idx(:)
   real(dp), allocatable :: eph_g2(:), g2(:,:,:)
   complex(dp), allocatable :: uq(:,:,:), gp(:,:), g_kerp(:,:,:,:,:), gkq(:,:,:)
   !
   integer(HID_T) :: file_id
   character(len=120) :: fname, dset_name, msg
   character(len=6), external :: int_to_char
   ! adaptive smear   
   real(dp), allocatable :: adsmear(:)
   real(dp) :: adsmear_tmp
   real(dp) :: vk(3), vkq(3)

   call start_clock('compute_g2')
   !sanity check
   if(elph%nb .ne. el%nb .or. elph%na .ne. ph%na) &
      call errore('compute_scatter_eph_g2','array dimension mismatch',1)
   nmod = ph%nm
   numq = qg%npts

   if( use_mem ) allocate( uq(nmod, nmod, numq) )
   if(elph%lpol) allocate( gp(nmod, numq) )
   if(elph%lpol .or. use_mem) then
!$omp parallel do schedule(guided) default(shared) private(iq, xq, wq, uqt)
      do iq = 1, numq
         xq = num2kpt(qg%list(iq), qg%ndim)
         call solve_phonon_modes(ph, xq, wq, uqt)
         if(elph%lpol) call eph_wan_longrange(elph%pol, xq, gp(:, iq), uqt)
         if(use_mem)  uq(:,:, iq) = uqt(:,:)
      enddo
!$omp end parallel do
   endif

   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_eph_g2_p"
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   call hdf_open_file(file_id, trim(fname), status='NEW')

   ! allocate work space
   allocate( g_kerp(3, elph%max_nrp, elph%nb, elph%nb, elph%na) )
   ! print out progress info
   call progressbar_init('Computing EPhMatrix:')
   ! main loop
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      nkq = kq_pair(i)%nkq
      ! skip if nkq = 0
      if( nkq < 1 ) cycle

      !kpts in crystal coordinate
      xk = num2kpt( kg%kpt(ik), kg%ndim )
      call solve_eigenvalue_vector(el, xk, ek, uk)
      call eph_fourier_el_para(elph, xk, g_kerp)
      !
      allocate( ist(nkq) )
      ist(1) = 0
      do j = 1, nkq-1
         ist(j+1) = ist(j) + kq_pair(i)%nchl(j)
      enddo
      !allocate space
      tot_chnl = sum(kq_pair(i)%nchl(:))
      allocate( bnd_idx(tot_chnl), eph_g2(tot_chnl) )
      bnd_idx = 0;   eph_g2 = 0.0_dp

      ! implement adaptive smearing 
      if ( adapt_smear_eph ) then
         allocate( adsmear(tot_chnl) )
         adsmear = 0.0_dp
      endif

!$omp parallel default(shared) private(j, iq, ikq, xkq, ekq, ukq, xq, wq, uqt, &
!$omp& qidx, gpol, gkq, g2, ltable, it, im, n, m, bands, adsmear_tmp, vk, vkq)
      allocate( gkq(el%nb, el%nb, nmod), g2(el%nb, el%nb, nmod) )
!$omp do schedule(guided)
      do j = 1, nkq
         iq = kq_pair(i)%iq(j)
         ikq = kq_pair(i)%ikq(j)
         !
         xkq = num2kpt(kg%kpt(ikq), kg%ndim)
         call solve_eigenvalue_vector(el, xkq, ekq, ukq)
         !
         xq = num2kpt( qg%list(iq), qg%ndim )
         if(use_mem) then
            uqt(:,:) = uq(:,:, iq)
         else
            call solve_phonon_modes(ph, xq, wq, uqt)
         endif
         !
         gpol = (0.0_dp,0.0_dp)
         !NOTE: qg%list(iq) could be (k+q)-k or k-(k+q)
         !  Here we need to check which case it is
         qidx= kpts_minus( kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim)
         ! iq -> (k+q) - k
         ! only invert if it is carriers only
         if ( .not. dyna_phph) then
            if( qidx .eq. qg%list(iq) ) then
               !uqt(:,:) = uq(:,:, iq)
               if(elph%lpol) gpol(:) = gp(:,iq)
            ! iq -> k - (k+q) .e.g -q
            else if( kpts_plus(qidx, qg%list(iq), qg%ndim, qg%ndim) .eq. 0 ) then
               !NOTE: it's important to use xq = -xq; alternatively we can
               ! use xq'=num2kpt(qidx), but xq' might be equal to xq + G
               ! (reciprocal lattice vectors). a small numerical noise due
               ! to G might break the relation u(q)* = u(-q);
               ! syp - explained in boltz_utils.f90/num2kpt
               xq  = -xq
               !u(q)* = u(-q);  g^L(q)* = g^L(-q)
               uqt(:,:) = dconjg( uqt(:,:) )
               if(elph%lpol) gpol(:) = dconjg( gp(:,iq) )
            else
               !write(stdout,*) kpts_plus(qidx, qg%list(iq), qg%ndim, qg%ndim), qg%list(iq), qidx
               write(msg,'(2(3f12.6,2x))') num2kpt(qidx, qg%ndim), xq
               call errore('compute_scatter_eph_g2','failed: ' // trim(msg), 1)
            endif
         endif

         !get e-ph matrix elements in wannier gauge and cartesian coord.
         call eph_fourier_elph(elph, xq, g_kerp, gkq)
         call eph_transform_fast(elph, uqt, uk, ukq, gkq, gpol)
         !compute |g|^2
         g2 = abs2( gkq )
         !
         call match_table(kg%enk(:,ikq), kg%enk(:,ik), qg%freq(:,iq), wcut, e_thr, ltable)
         !sanity check
         if( kq_pair(i)%nchl(j) .ne. count(ltable) ) &
            call errore('compute_scatter_eph_g2','kq_pair(i)%nchl(j) .ne. count(ltable)',1)

         it = ist(j)
         do im = 1, nmod
         do n = 1, kg%numb
         do m = 1, kg%numb
            if( ltable(m,n,im) ) then
               it = it + 1
               bands = (/m, n, im/)
               bnd_idx(it) = band2idx(bands, kg%numb)
               !syp - the denominator is 2*omega_q following Eq. (24) of Perturbo CPC paper
               eph_g2(it) = g2(m+kg%bmin-1, n+kg%bmin-1, im) / (2.0_dp*qg%freq(im,iq))
               if ( adapt_smear_eph ) then
                  vk = kg%vnk(:,n,ik) ! try vq first
                  vkq = kg%vnk(:,m,ikq)
                  adsmear_tmp = sqrt(2.0_dp) * calc_adaptive_smear(vk,vkq,maxval(kg%ndim))
                  adsmear_tmp = min(adsmear_tmp, delta_smear)
                  ! min smear = max smear * 0.05
                  adsmear_tmp = max(adsmear_tmp, delta_smear*0.05_dp)
                  adsmear(it) = adsmear_tmp
               endif
            endif
         enddo; enddo; enddo
      enddo
!$omp end do
      deallocate(gkq, g2)
!$omp end parallel
      !
      ! output to hdf5 files
      dset_name = "bands_index_" // trim( int_to_char(i) )
      call hdf_write_dataset(file_id, trim(dset_name), bnd_idx)
      dset_name = "eph_g2_" // trim( int_to_char(i) )
      call hdf_write_dataset(file_id, trim(dset_name), eph_g2)
      if ( adapt_smear_eph ) then
         dset_name = "adapt_smear_" // trim( int_to_char(i) )
         call hdf_write_dataset(file_id, trim(dset_name), adsmear)
      endif
      ! show progress
      call progressbar(i, nk_loc)
      !
      deallocate(bnd_idx, eph_g2, ist)
      if ( adapt_smear_eph ) deallocate(adsmear)
   enddo
   call hdf_close_file(file_id)

   deallocate(g_kerp)
   if( use_mem ) deallocate( uq )
   if(elph%lpol) deallocate( gp )

   call mp_barrier(inter_pool_comm)
   call stop_clock('compute_g2')
end subroutine compute_scatter_eph_g2

! kyao not finished
subroutine boltz_ph3scat_setup_fast(qg, ph, pht, isfine_in)
   use boltz_utils, only : num2kpt, kpt2num
   use pert_param, only: symm 
   implicit none
   type(grid), intent(in) :: qg
   type(lattice_ifc), intent(in) :: ph
   type(lattice_ifct), intent(in) :: pht
   logical, intent(in), optional :: isfine_in
   real(dp) :: wcut, e_thr
   integer :: nkpt, i, j, ik, collect_pools(npool), nq_col(npool), k, collect_pools_chl(npool)
   !dimension for the q-grid. q-should be commensurate with the k-grid
   integer :: iq_st, iq_end
   integer :: qndim(3), nkpt_c
   integer :: ierr
   logical :: isfine

   call start_clock('phscat_setup')
   write(stdout,'(5x,a20)') '>start ph3scat_setup'
   isfine = .false.
   if (present(isfine_in)) then
      isfine = isfine_in
      write(stdout,'(5x,a)') 'interpolate from a coarser grid:'
   endif

   !!!! TEMP CHANGE by kyao
   do i = 1, 3
      qndim(i) = qg%ndim(i) / 3
   enddo
   !init
   nkpt = qg%nk
   if (isfine) then
      !!!! TEMP CHANGE by kyao
      !nkpt_c = nkpt / 8
      nkpt_c = qg%nk_irr
   elseif(symm)then
      nkpt_c = qg%nk_irr
   else
      nkpt_c = nkpt
   endif
   wcut  = phfreq_cutoff_ph
   e_thr = delta_smear_ph*3.0_dp  !exp(-3*3) ~ 10^-4

   !syp@20241005 [TOCHECK] 
   call mp_split_pools(nkpt_c, iq_st, iq_end, nq_loc)
   !interleave distribution of qpoints over pools, for better load balance
   !nq_loc = nkpt/npool + merge(1, 0, my_pool_id < mod(nkpt, npool))
   !write iq_st, iq_end, nq_loc
   allocate( qq_pair(nq_loc) )
   i = 0
   do ik = iq_st, iq_end
      i = i + 1
      if (isfine) then
         !!!! TEMP CHANGE by kyao
         qq_pair(i)%ik = qg%ir2kpt(ik)
         !qq_pair(i)%ik = kpt2num(num2kpt(ik-1, qndim),qg%ndim) + 1
      elseif(symm)then
         ! here we need to assign the place of each irreducible k-point 
         ! in the qg%nk kpoints list (the selected reducible kpoints)
         ! to the qq_pair(i)%ik
         qq_pair(i)%ik = qg%ir2kpt(ik)
      else
         qq_pair(i)%ik = ik
      endif
   enddo
   !debug
   if(i .ne. nq_loc) call errore('boltz_ph3scat_setup', 'i .ne. nq_loc', 1)

   !generate (q,q) pair (initialized qq_pair)
   !else
   if (isfine) then
      call setup_qq_pair(qg, ph, wcut, e_thr, nq_loc, qq_pair, .true.)
   else
      call setup_qq_pair(qg, ph, wcut, e_thr, nq_loc, qq_pair)
   endif
   !
   num_qq_pair = 0
   ph3_num_scat = 0
   do i = 1, nq_loc
      if (qq_pair(i)%nkq < 1) cycle
      num_qq_pair = num_qq_pair + qq_pair(i)%nkq
      ph3_num_scat = ph3_num_scat + sum( qq_pair(i)%nchl(:) )
   enddo

   nq_col = 0                                                                                                                 
   nq_col( my_pool_id+1 ) = nq_loc                                                                                            
   call mp_sum( nq_col, inter_pool_comm ) 

   collect_pools = 0
   collect_pools( my_pool_id+1 ) = num_qq_pair
   call mp_sum( collect_pools, inter_pool_comm )

   !sypeng added
   collect_pools_chl = 0
   collect_pools_chl( my_pool_id+1 ) = ph3_num_scat
   call mp_sum( collect_pools_chl, inter_pool_comm )

   !output info
   write(stdout, '(5x, a)') "No. of q-points and (q,q) pairs on each processor:"
   num_qq_pair_total = 0
   do i = 1, npool
     !num_qq_pair_total = num_qq_pair_total + collect_pools(i) 
     !write(stdout, '(5x, a, i4, a, i12)')  &
     !   "Proc.", i, ' |  #.(q,q)', collect_pools(i)
      write(stdout, '(5x, a, i4, a, i9, a, i12, a, i12)')  &
         "Proc.", i, ':  #.irr-q ', nq_col(i), " |  #.(q,q')", collect_pools(i), " |  #.channels-q", collect_pools_chl(i)
   enddo
   ! Log data-structure allocation details into the YAML output file.
   if (ionode) then
      write(ymlout, '(/,a)') 'boltz_ph3_scatter memory:'
      call output_qq_pair_alloc_stats_to_yaml()
   endif

   if(.not. load_scatter_ph3) then
      !compute psi2 and sove them to temperary files
      call compute_scatter_ph3_psi2(qg, ph, pht, wcut, e_thr, use_mem)
   endif

   ! setup 2d array for scatter_ph
   call load_scatter_channels_ph()
   if (ionode) then
      write(ymlout, '(/,a)') 'boltz_ph3_scatter memory:'
      call output_scatter_ph_alloc_stats_to_yaml()
      call output_tscat_ph_alloc_stats_to_yaml()
   endif

   write(stdout,'(5x,a20)') '>ph3scat_setup done.'
   call stop_clock('phscat_setup')
end subroutine boltz_ph3scat_setup_fast


!syp - [finished]
subroutine boltz_ph3scat_setup(qg, ph, pht, isfine_in)
   use boltz_utils, only : num2kpt, kpt2num
   use pert_param, only: symm 
   implicit none
   type(grid), intent(in) :: qg
   type(lattice_ifc), intent(in) :: ph
   type(lattice_ifct), intent(in) :: pht
   logical, intent(in), optional :: isfine_in
   real(dp) :: wcut, e_thr
   integer :: nkpt, i, j, ik, collect_pools(npool), nq_col(npool), k, collect_pools_chl(npool)
   !dimension for the q-grid. q-should be commensurate with the k-grid
   integer :: iq_st, iq_end
   integer :: qndim(3), nkpt_c
   integer :: ierr
   logical :: isfine

   call start_clock('phscat_setup')
   write(stdout,'(5x,a20)') '>start ph3scat_setup'
   isfine = .false.
   if (present(isfine_in)) then
      isfine = isfine_in
      write(stdout,'(5x,a)') 'interpolate from a coarser grid:'
   endif

   do i = 1, 3
      qndim(i) = qg%ndim(i) / 2
   enddo
   !init
   nkpt = qg%nk
   if (isfine) then
      nkpt_c = nkpt / 8
   elseif(symm)then
      nkpt_c = qg%nk_irr
   else
      nkpt_c = nkpt
   endif
   wcut  = phfreq_cutoff_ph
   e_thr = delta_smear_ph*3.0_dp  !exp(-3*3) ~ 10^-4

   !syp@20241005 [TOCHECK] 
   call mp_split_pools(nkpt_c, iq_st, iq_end, nq_loc)
   !interleave distribution of qpoints over pools, for better load balance
   !nq_loc = nkpt/npool + merge(1, 0, my_pool_id < mod(nkpt, npool))
   !write iq_st, iq_end, nq_loc
   allocate( qq_pair(nq_loc) )
   i = 0
   do ik = iq_st, iq_end
      i = i + 1
      if (isfine) then
         qq_pair(i)%ik = kpt2num(num2kpt(ik-1, qndim),qg%ndim) + 1
      elseif(symm)then
         ! here we need to assign the place of each irreducible k-point 
         ! in the qg%nk kpoints list (the selected reducible kpoints)
         ! to the qq_pair(i)%ik
         qq_pair(i)%ik = qg%ir2kpt(ik)
      else
         qq_pair(i)%ik = ik
      endif
   enddo
   !debug
   if(i .ne. nq_loc) call errore('boltz_ph3scat_setup', 'i .ne. nq_loc', 1)

   !generate (q,q) pair (initialized qq_pair)
   !else
   if (isfine) then
      call setup_qq_pair(qg, ph, wcut, e_thr, nq_loc, qq_pair, .true.)
   else
      call setup_qq_pair(qg, ph, wcut, e_thr, nq_loc, qq_pair)
   endif
   !
   num_qq_pair = 0
   ph3_num_scat = 0
   do i = 1, nq_loc
      if (qq_pair(i)%nkq < 1) cycle
      num_qq_pair = num_qq_pair + qq_pair(i)%nkq
      ph3_num_scat = ph3_num_scat + sum( qq_pair(i)%nchl(:) )
   enddo

   nq_col = 0                                                                                                                 
   nq_col( my_pool_id+1 ) = nq_loc                                                                                            
   call mp_sum( nq_col, inter_pool_comm ) 

   collect_pools = 0
   collect_pools( my_pool_id+1 ) = num_qq_pair
   call mp_sum( collect_pools, inter_pool_comm )

   !sypeng added
   collect_pools_chl = 0
   collect_pools_chl( my_pool_id+1 ) = ph3_num_scat
   call mp_sum( collect_pools_chl, inter_pool_comm )

   !output info
   write(stdout, '(5x, a)') "No. of q-points and (q,q) pairs on each processor:"
   num_qq_pair_total = 0
   do i = 1, npool
     !num_qq_pair_total = num_qq_pair_total + collect_pools(i) 
     !write(stdout, '(5x, a, i4, a, i12)')  &
     !   "Proc.", i, ' |  #.(q,q)', collect_pools(i)
      write(stdout, '(5x, a, i4, a, i9, a, i12, a, i12)')  &
         "Proc.", i, ':  #.irr-q ', nq_col(i), " |  #.(q,q')", collect_pools(i), " |  #.channels-q", collect_pools_chl(i)
   enddo
   ! Log data-structure allocation details into the YAML output file.
   if (ionode) then
      write(ymlout, '(/,a)') 'boltz_ph3_scatter memory:'
      call output_qq_pair_alloc_stats_to_yaml()
   endif

   if(.not. load_scatter_ph3) then
      !compute psi2 and sove them to temperary files
      call compute_scatter_ph3_psi2(qg, ph, pht, wcut, e_thr, use_mem)
   endif

   ! setup 2d array for scatter_ph
   call load_scatter_channels_ph()
   if (ionode) then
      write(ymlout, '(/,a)') 'boltz_ph3_scatter memory:'
      call output_scatter_ph_alloc_stats_to_yaml()
      call output_tscat_ph_alloc_stats_to_yaml()
   endif

   write(stdout,'(5x,a20)') '>ph3scat_setup done.'
   call stop_clock('phscat_setup')
end subroutine boltz_ph3scat_setup


! Open the HDF5 file generated for this MPI pool (i.e. MPI node) that
! contains the eph_g2 and bands_index values.  The file-handle is
! returned back in the file_id out-parameter.
!
! The subroutine doesn't report any errors because if the file can't
! be opened then the subroutine halts the program with an error.
subroutine open_pool_eph_g2_file(file_id)
   use hdf5_utils
   use pert_param, only: prefix, tmp_dir
   implicit none
   !
   ! For reading the HDF5 file data
   integer(HID_T), intent(out) :: file_id
   logical :: has_file
   character(len=120) :: fname
   character(len=6), external :: int_to_char

   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_eph_g2_p"
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   ! check if file exist
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('load_scatter_eph_g2', "Missing file: "//trim(fname), 1)
   !
   ! open file
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
end subroutine open_pool_eph_g2_file

! Read the eph_g2 array and bands_index array for a given k-grid location
! from the HDF5 file generated for this MPI pool (i.e. MPI node).  The
! file-handle and k-grid location are passed in as input-parameters.  The
! bnd_idx and eph_g2 pointers are expected to point to arrays that are
! large enough to receive the values stored in the file.
subroutine read_pool_eph_g2_data(file_id, i_kloc, bnd_idx, eph_g2, adsmear)
   use hdf5_utils
   use pert_param, only: adapt_smear_eph
   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: i_kloc
   integer, allocatable, intent(inout) :: bnd_idx(:)
   real(dp), allocatable, intent(inout) :: eph_g2(:)
   real(dp), allocatable, intent(inout), optional :: adsmear(:)
   !
   character(len=120) :: dset_name
   character(len=6), external :: int_to_char

   dset_name = "bands_index_" // trim( int_to_char(i_kloc) )
   call hdf_read_dataset(file_id, trim(dset_name), bnd_idx)
   dset_name = "eph_g2_" // trim( int_to_char(i_kloc) )
   call hdf_read_dataset(file_id, trim(dset_name), eph_g2)
   if ( adapt_smear_eph ) then
      dset_name = "adapt_smear_" // trim( int_to_char(i_kloc) )
      call hdf_read_dataset(file_id, trim(dset_name), adsmear)
   endif
end subroutine


subroutine setup_scatter(kg, el, nmodes)
   use hdf5_utils
   use boltz_utils, only: idx2band
   use qe_mpi_mod, only: ionode, stdout
   use pert_param, only: prefix, tmp_dir, ph_mode_exclude_ranges, &
      adapt_smear_eph
   implicit none

   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   integer, intent(in) :: nmodes

   !local
   integer :: i_scat, tot_chl
   integer(kind=i8b) :: i_scat_chl
   integer :: n, i, ik, nkq, ikq, iq, j, nchl, it, ib
#if defined(STORE_G2DT)
   ! Variables needed to compute -g2*dt1 and -g2*dt2
   integer :: m
   real(dp) :: wq, enk, emkq, g2, dt1, dt2
#endif
   !
   ! For reading the HDF5 file data
   integer(HID_T) :: file_id
   integer, allocatable :: bnd_idx(:) ! For reading band-indexes from HDF5 file
   real(dp), allocatable :: eph_g2(:) ! For reading eph-g2 values from HDF5 file
   real(dp), allocatable :: adsmear(:)
   !
   ! Exclude ph mode variables
   integer :: im, jm, im1, im2
   integer :: bands(3)
   logical :: include_phonon_modes(nmodes)

   call start_clock('scat_setup')
   write(stdout,'(5x,a20)') '>start scatter_setup'

   allocate( scatter(num_kq_pairs) )
   allocate( scatter_channels(num_scatter_channels) )
   if ( adapt_smear_eph) then
      allocate(tsmear_list(num_scatter_channels))
   endif

   ! TODO(donnie):  Remove this output?
#if defined(SCAT_FWD)
   write (*, '(5x, a)') ' * scatter will support forward-traversal'
#endif

#if defined(SCAT_REV)
   write (*, '(5x, a)') ' * scatter will support reverse-traversal'
#endif

   write (*, '(5x, a, i12)') ' * Number of kq-pairs:         ', num_kq_pairs
   write (*, '(5x, a, i12)') ' * Number of scatter-channels: ', num_scatter_channels

   !read g2 from files

   call start_clock("load_eph_g2")

   call open_pool_eph_g2_file(file_id)

   i_scat_chl = 0
   i_scat = 0
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      nkq = kq_pair(i)%nkq
      !skip if nkq = 0
      if(nkq < 1) cycle

      ! Arrays for reading values from HDF5 file
      tot_chl = sum( kq_pair(i)%nchl(:) )
      allocate( bnd_idx(tot_chl), eph_g2(tot_chl) )
      bnd_idx = 0
      eph_g2 = 0.0_dp

      ! Read bnd_idx and eph_g2 from HDF5 file
      if ( adapt_smear_eph ) then
         allocate( adsmear(tot_chl) )
         adsmear = 0.0_dp
         call read_pool_eph_g2_data(file_id, i, bnd_idx, eph_g2, adsmear)
      else
         call read_pool_eph_g2_data(file_id, i, bnd_idx, eph_g2)
      endif

      it = 0
      do j = 1, nkq
         i_scat = i_scat + 1
         iq = kq_pair(i)%iq(j)
         ikq = kq_pair(i)%ikq(j)
         !
         scatter(i_scat)%ik = ik
         scatter(i_scat)%iq = iq
         scatter(i_scat)%ikq = ikq

         nchl = kq_pair(i)%nchl(j)

#if defined(SCAT_FWD)
         ! must add 1 since i_scat_chl is incremented in the innermost loop
         scatter(i_scat)%chl_start_index = i_scat_chl + 1
         scatter(i_scat)%nchl = nchl
#endif

         do ib = 1, nchl
            i_scat_chl = i_scat_chl + 1
            it = it + 1
#if defined(SCAT_REV)
            scatter_channels(i_scat_chl)%scat_index = i_scat
#endif
            scatter_channels(i_scat_chl)%bands_index = bnd_idx(it)
            if (adapt_smear_eph ) then
               tsmear_list(i_scat_chl) = adsmear(it)
            endif
#if defined(STORE_G2DT)
            ! Precompute and store the -g2*dt1 and -g2*dt2 terms, since they
            ! don't change from iteration to iteration.

            bands = idx2band( bnd_idx(it), kg%numb )
            m = bands(1);  n = bands(2);  im = bands(3)

            ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
            g2 = eph_g2(it)

            ! eigenvalue
            wq   = qgrid%freq(im, iq)
            enk  = kg%enk(n, ik)
            emkq = kg%enk(m, ikq)

            ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
            ! for dyna_phph, need qg%kweight here
            ! we compute dt1, dt2, e-ph
            ! for ph-e *kg%kweight/qg%kweight (in boltz_scatter_integral kelly unsure) 
            ! here is fine because type(phonon_grid) qgrid
            ! is generated even for dyna_phph = true   
            ! and qgrid%weight = qg%kweight 
            if ( adapt_smear_eph ) then                
               dt1 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq + wq)
               dt2 = qgrid%weight * pi * 2.0_dp * gauss(tsmear_list(i), enk - emkq - wq)
            else
               ! delta(Enk - Emkq + wq)
               dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
               ! delta(Enk - Emkq -wq)
               dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
            endif

            scatter_channels(i_scat_chl)%g2_dt1 = g2 * dt1
            scatter_channels(i_scat_chl)%g2_dt2 = g2 * dt2
            !sypeng: to be compatible with Kelly's additional treatment of dt for ph-e coupling
            !sypeng: TODO
            scatter_channels(i_scat_chl)%eph_g2 = eph_g2(it)
#else
            scatter_channels(i_scat_chl)%eph_g2 = eph_g2(it)
#endif
         enddo
      enddo
      deallocate(bnd_idx, eph_g2)
      if ( adapt_smear_eph ) deallocate( adsmear )
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

      ! Iterate over all scatter-channels to zero out eph_g2 values for
      ! phonon modes to be excluded.
!$omp parallel do schedule(guided) default(shared) private(i_scat_chl, bands, im)
      do i_scat_chl = 1, num_scatter_channels
         bands = idx2band( scatter_channels(i_scat_chl)%bands_index, kg%numb )
         im = bands(3)

         ! Set the g2 elements to zero for the phonon modes to exclude
         if (.not. include_phonon_modes(im)) then
#if !defined(STORE_G2DT)
            scatter_channels(i_scat_chl)%eph_g2 = 0.0_dp
#else
            scatter_channels(i_scat_chl)%g2_dt1 = 0.0_dp
            scatter_channels(i_scat_chl)%g2_dt2 = 0.0_dp
            !sypeng: to be compatible with Kelly's additional treatment of dt for ph-e coupling
            !sypeng: TODO
            scatter_channels(i_scat_chl)%eph_g2 = 0.0_dp
#endif
         endif
      enddo
!$omp end parallel do

      call stop_clock("exclude_ph_modes")
   end if

   !
   write(stdout,'(5x,a20)') '>scatter_setup done.'
   call stop_clock('scat_setup')
   flush(stdout)

end subroutine setup_scatter


!> Based on the array from the input file ph_mode_exclude_ranges
!! build the include_phonon_modes bollean
!! array. If an element of this array is false, this phonon mode
!! will be omitted in the scat: g2 will be set to zero.
subroutine interpret_ph_exclude_ranges(nmodes, include_phonon_modes)
   use yaml_utils, only: ymlout
   use qe_mpi_mod, only: ionode, stdout
   use pert_param, only: ph_mode_exclude_ranges
   implicit none
   integer, intent(in) :: nmodes
   logical, intent(inout) :: include_phonon_modes(nmodes)
   !
   integer :: im, jm, im1, im2
   character(len=6), external :: int_to_char
   character(len=120) :: msg

   include_phonon_modes(:) = .true.

   do im = 1, size(ph_mode_exclude_ranges, 1)

      im1 = ph_mode_exclude_ranges(im, 1)
      if( im1 > 0 ) then

         im2 = ph_mode_exclude_ranges(im, 2)
         if ( im1 > nmodes .or. im2 > nmodes ) then
            write(msg, '(a,I4)') &
            'exclude mode indices must be smaller than the number of phonon modes:', nmodes
            call errore('load_scatter_eph_g2', trim(msg), 1)

         else if (im2 < im1) then
            write(msg, '(a,I2,a,I2,a)') &
            'ph_mode_exclude_ranges(',im, &
            ', 2) is smaller than ph_mode_exclude_ranges(',im, ', 1)'
            call errore('load_scatter_eph_g2', trim(msg), 1)
         end if

         do jm = im1, im2
            include_phonon_modes(jm) = .false.
         end do

      end if
   end do

   ! output which phonon modes will be included or excluded
   if (ionode) then
      write(stdout, '(/,8x,a)') 'mode number | included/excluded (+/-)'
      write(stdout, '(8x,a)')   '-------------------------------------'
      do im = 1, nmodes
         if(include_phonon_modes(im)) then
            write(stdout, '(10x,I8,2x,a,12x,a)') im, '|', '+'
         else
            write(stdout, '(10x,I8,2x,a,12x,a)') im, '|', '-'
         end if
      end do
      write(stdout, '(8x,a,/)')   '-------------------------------------'

      ! output phonon mode include into YAML
      write(ymlout, '(/,a)') 'include phonon modes:'
      do im = 1, nmodes
         if(include_phonon_modes(im)) then
            write(ymlout, '(6x,a,a,2x,a)') trim(int_to_char(im)), ':', 'True' 
         else
            write(ymlout, '(6x,a,a,2x,a)') trim(int_to_char(im)), ':', 'False' 
         end if
      end do
   end if
end subroutine

subroutine load_scatter_channels_ph()
   use hdf5_utils
   use pert_param, only: prefix, tmp_dir, adapt_smear_phph
   implicit none

   integer(kind=i8b) :: i_scat_chl
   integer :: i_scat, i, ik, nkq, j, tot_chnl, nchl, it, ib
   integer, allocatable :: bnd_idx(:)
   real(dp), allocatable :: ph3_psi2(:)
   real(dp), allocatable :: adsmear(:)
   !
   logical :: has_file
   integer(HID_T) :: file_id
   character(len=120) :: fname, dset_name
   character(len=6), external :: int_to_char
   
   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_ph3_psi2_p" 
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   ! check if file exist
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('load_scatter_channels_ph', "Missing file: "//trim(fname), 1)
   ! open file
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')

   ! initialize an array direct array 
   allocate(scatter_ph(2, num_qq_pair))
   allocate(tscat_ph(ph3_num_scat))
   if ( adapt_smear_phph) then
      allocate(tsmear_ph_list(ph3_num_scat))
   endif
   i_scat_chl = 0
   i_scat = 0
   ! for each ik
   do i = 1, nq_loc
      nkq = qq_pair(i)%nkq
      if (nkq < 1) cycle
      ik = qq_pair(i)%ik

      tot_chnl = sum( qq_pair(i)%nchl(:) )
      if (tot_chnl < 1) then
         write(*,*) 'no channels'
      endif
      allocate( bnd_idx(tot_chnl), ph3_psi2(tot_chnl) )
      bnd_idx = 0;   ph3_psi2 = 0.0_dp
      if ( adapt_smear_phph ) then
         allocate( adsmear(tot_chnl) )
         adsmear = 0.0_dp
      endif

      !read bnd_idx and ph3_g2 from hdf5 file
      dset_name = "bands_index_" // trim( int_to_char(i) )
      call hdf_read_dataset(file_id, trim(dset_name), bnd_idx)
      dset_name = "ph3_psi2_" // trim( int_to_char(i) )
      call hdf_read_dataset(file_id, trim(dset_name), ph3_psi2)
      if ( adapt_smear_phph ) then
         dset_name = "adapt_smear_" // trim( int_to_char(i) )
         call hdf_read_dataset(file_id, trim(dset_name), adsmear)
      endif

      it = 0
      ! for each iq
      do j = 1, nkq
         i_scat = i_scat + 1
         nchl = qq_pair(i)%nchl(j)

         scatter_ph(1,i_scat) = ik
         scatter_ph(2,i_scat) = qq_pair(i)%iq(j)
         ! for each channel
         do ib = 1, nchl
            i_scat_chl = i_scat_chl + 1
            it = it + 1
            !write(*,*) ik, it, i_scat, i_scat_chl 
            ! fill scat_index, band_index, and matrix element g
            tscat_ph(i_scat_chl)%scat_index = i_scat
            tscat_ph(i_scat_chl)%bands_index = bnd_idx(it) 
            tscat_ph(i_scat_chl)%g2 = ph3_psi2(it) 
      
            if (adapt_smear_phph ) then
               tsmear_ph_list(i_scat_chl) = adsmear(it)
            endif
         enddo
      enddo
      ! deallocate right after 
      deallocate( qq_pair(i)%iq, qq_pair(i)%nchl )
      deallocate(bnd_idx, ph3_psi2)
      if ( adapt_smear_phph ) deallocate( adsmear )
      !
   enddo
   deallocate( qq_pair )
   !
   call hdf_close_file(file_id)
   call mp_barrier( inter_pool_comm )
end subroutine load_scatter_channels_ph


!------
! Functions for reporting memory usage of data structures

! Output details regarding the memory allocation of the kq-pair
! data structure.
subroutine output_kq_pair_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! kq_pair

   bytes_per_entry = STORAGE_SIZE(kq_pair(1)) / 8
   total_bytes = bytes_per_entry * SIZE(kq_pair, kind=8)

   ! subarrays of kq_pair

   total_allocs = 1
   do i = 1, nk_loc
      if (kq_pair(i)%nkq == 0) cycle

      ! Each element in kq_pair has three subarrays
      total_allocs = total_allocs + 3
      total_bytes = total_bytes + &
         SIZE(kq_pair(i)%iq  , kind=8) * STORAGE_SIZE(kq_pair(i)%iq(1)  ) / 8 + &
         SIZE(kq_pair(i)%ikq , kind=8) * STORAGE_SIZE(kq_pair(i)%ikq(1) ) / 8 + &
         SIZE(kq_pair(i)%nchl, kind=8) * STORAGE_SIZE(kq_pair(i)%nchl(1)) / 8
   enddo

   write(ymlout,'(3x, a)') 'kq_pair [intermediate]:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', nk_loc
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i12/)') 'allocations: ', total_allocs

   call boltzscat_accum_cpu_bytes(total_bytes)
end subroutine output_kq_pair_alloc_stats_to_yaml

! Output details regarding the memory allocation of the scatter
!------
! Functions for reporting memory usage of data structures

! Output details regarding the memory allocation of the qq-pair
! data structure.
subroutine output_qq_pair_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! qq_pair

   bytes_per_entry = STORAGE_SIZE(qq_pair(1)) / 8
   total_bytes = bytes_per_entry * SIZE(qq_pair, kind=8)

   ! subarrays of qq_pair

   total_allocs = 1
   do i = 1, nq_loc
      if (qq_pair(i)%nkq == 0) cycle

      ! Each element in qq_pair has three subarrays
      total_allocs = total_allocs + 2
      total_bytes = total_bytes + &
         SIZE(qq_pair(i)%iq  , kind=8) * STORAGE_SIZE(qq_pair(i)%iq(1)  ) / 8 + &
         SIZE(qq_pair(i)%nchl, kind=8) * STORAGE_SIZE(qq_pair(i)%nchl(1)) / 8
   enddo

   write(ymlout,'(3x, a)') 'qq_pair [intermediate]:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', nq_loc
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i8/)') 'allocations: ', total_allocs
end subroutine output_qq_pair_alloc_stats_to_yaml

! Output details regarding the memory allocation of the scat
! data structure.
!subroutine output_ph3_scat_alloc_stats_to_yaml()
!   implicit none
!   integer :: i, total_allocs
!   integer(8) :: bytes_per_entry, total_bytes
!
!   ! scat
!
!   bytes_per_entry = STORAGE_SIZE(ph3_scat(1)) / 8
!   total_bytes = bytes_per_entry * SIZE(ph3_scat, kind=8)
!
!   ! subarrays of scat
!   total_allocs = 1
!   do i = 1, num_qq_pair
!      if (ph3_scat(i)%nchl == 0) cycle
!
!      ! Each element in scat has two subarrays
!      total_allocs = total_allocs + 2
!      total_bytes = total_bytes + &
!         SIZE(ph3_scat(i)%bands_index, kind=8) * STORAGE_SIZE(ph3_scat(i)%bands_index(1)) / 8 + &
!         SIZE(ph3_scat(i)%eph_g2     , kind=8) * STORAGE_SIZE(ph3_scat(i)%eph_g2(1)     ) / 8
!   enddo
!
!   write(ymlout,'(3x, a)') 'ph3_scat:'
!   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
!   write(ymlout,'(6x, a, i12)') 'entries: ', num_qq_pair
!   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
!   write(ymlout,'(6x, a, i8/)') 'allocations: ', total_allocs
!end subroutine output_ph3_scat_alloc_stats_to_yaml

! Output details regarding the memory allocation of the scatter_ph
! data structure.
subroutine output_scatter_ph_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! scat

   bytes_per_entry = (STORAGE_SIZE(scatter_ph(1,1)) + STORAGE_SIZE(scatter_ph(2,1))) / 8
   total_bytes = bytes_per_entry * SIZE(scatter_ph, kind=8)

   write(ymlout,'(3x, a)') 'scatter_ph:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_qq_pair
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
end subroutine output_scatter_ph_alloc_stats_to_yaml

! Output details regarding the memory allocation of the t_scat_ph
! data structure.
subroutine output_scatter_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes, subelems

   ! scatter

   bytes_per_entry = STORAGE_SIZE(scatter(1)) / 8
   total_bytes = bytes_per_entry * SIZE(scatter, kind=8)

   write(ymlout,'(3x, a)') 'scatter:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_kq_pairs
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i5/)') 'allocations: ', 1

   call boltzscat_accum_cpugpu_bytes(total_bytes)
end subroutine output_scatter_alloc_stats_to_yaml

! Output details regarding the memory allocation of the scatter_channels
! data structure.
subroutine output_scatter_channels_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes, subelems

   ! scatter_channels

   bytes_per_entry = STORAGE_SIZE(scatter_channels(1)) / 8
   total_bytes = bytes_per_entry * SIZE(scatter_channels, kind=8)

   write(ymlout,'(3x, a)') 'scatter_channels:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_scatter_channels
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i5/)') 'allocations: ', 1

   call boltzscat_accum_cpugpu_bytes(total_bytes)
end subroutine output_scatter_channels_alloc_stats_to_yaml

! Output details regarding the memory allocation of the t_scat_ph
! data structure.
subroutine output_tscat_ph_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! scat

   bytes_per_entry = STORAGE_SIZE(tscat_ph(1)) / 8
   total_bytes = bytes_per_entry * SIZE(tscat_ph, kind=8)

   write(ymlout,'(3x, a)') 'tscat_ph:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', ph3_num_scat
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
end subroutine output_tscat_ph_alloc_stats_to_yaml

subroutine compute_scatter_ph3_psi2(qg, ph, pht, wcut, e_thr, use_mem)
   use pert_param, only: prefix, tmp_dir, delta_smear_ph, adapt_smear_phph, symm
   use pert_utils, only: abs2, match_table, match_table_all, get_exp_ikr, pure_exp_ikr, &
      match_table_absorp_emit_ph, match_table_absorp_emit_ph_dyna
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus, &
      calc_adaptive_smear, kpts_plus_invert
   use pert_output, only: progressbar_init, progressbar
   use phonon_dispersion, only: solve_phi3mat, solve_phonon_velocity_findiff
   use qe_mpi_mod, only: stdout, ionode, world_comm, mp_barrier
   use force_constant_thirdorder, only: expikr_set, expikr, lattice_ifct
   use triplet_cell, only : trip_cell

   !
   !use pert constant
   use pert_const, only: ryd2mev, czero
   use hdf5_utils
   use pert_data, only: nat, tau

   use third_order_fc 

   implicit none
   type(grid), intent(in) :: qg
   type(lattice_ifc), intent(in) :: ph
   type(lattice_ifct), intent(in) :: pht
   !
   logical, intent(in) :: use_mem
   real(dp), intent(in) :: wcut, e_thr
   !
   type(trip_cell), pointer :: tc
   integer :: ierr, myrank, nprocs
   type(expikr_set) :: expikr_set1, expikr_set2, expikr_set_tot

   integer :: nmod, qidx, numq1
   integer :: i, j, ik, iq, ikq, it, im, n, m, nkq, tot_chnl, bands(3)
   integer :: ia, ja, ka, iexp, idx
   real(dp) :: xq(3), xk(3), xkq(3), wq(ph%nm), wk(ph%nm), wkq(ph%nm)
   complex(dp) :: uk(ph%nm, ph%nm), uq(ph%nm, ph%nm), ukq(ph%nm, ph%nm),ph3
   logical :: ltable(ph%nm, ph%nm, ph%nm)

   ! ph-ph-svd
   type(fct_matrix) :: fct
   complex(dp), allocatable :: phiqr(:,:,:,:)
   complex(dp), allocatable :: phi3_qq(:,:,:)
   !
   integer, allocatable :: ist(:), bnd_idx(:)
   integer, allocatable :: qq_idx(:,:)
   real(dp), allocatable :: ph3_psi2(:)
   real(dp), allocatable :: adsmear(:)
   complex(dp), allocatable :: uqgrid(:,:,:)
   complex(dp), allocatable :: exp_ikr1(:,:), exp_ikr2(:,:)
   complex(dp) :: expikr_tmp
   !
   integer(HID_T) :: file_id
   character(len=120) :: fname, dset_name, msg
   character(len=6), external :: int_to_char

   real(dp), allocatable :: cryst_tau(:,:)
   real(dp) :: vk(3), vkq(3)
   real(dp), allocatable :: velocity(:,:,:)
   real(dp) :: adsmear_tmp
   integer :: ongrid, ix
   real(dp) :: fweights(3)

   call start_clock('comp_ph3psi2')

   call load_fct(prefix, fct)
   allocate(phi3_qq(ph%nm,ph%nm,ph%nm), phiqr(fct%nr,ph%nm,ph%nm,ph%nm))
   phi3_qq = czero
   phiqr = czero


   fweights(1) = 1.0_dp/qg%ndim(1)
   fweights(2) = 1.0_dp/qg%ndim(2)
   fweights(3) = 1.0_dp/qg%ndim(3)
   allocate( cryst_tau(3, nat) )
   cryst_tau(:,:) = tau(:,:)
   !sanity check
   ! disable sanity check for now
   if(ph%nm .ne. pht%nm .or. ph%na .ne. pht%na .or. ph%nm .ne.qg%numb) &
      call errore('compute_scatter_ph3_psi2','array dimension mismatch',1)
   nmod = ph%nm

   numq1 = qg%nk;
   if( use_mem ) then
      allocate( uqgrid(nmod, nmod, qg%nk))
   endif
  
!$omp parallel do schedule(guided) default(shared) private(iq, xq, wq, uq, ia, ja, ka, idx, tc)
      do iq = 1, numq1
         xq = num2kpt(qg%kpt(iq), qg%ndim)
         call solve_phonon_modes(ph, xq, wq, uq)
         if(use_mem)  uqgrid(:,:, iq) = uq(:,:)
      enddo
!$omp end parallel do


   
   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_ph3_psi2_p" 
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   call hdf_open_file(file_id, trim(fname), status='NEW')

   ! print out progress info
   call progressbar_init('Computing Ph3Matrix:')
   ! main loop
   do i = 1, nq_loc
      ik = qq_pair(i)%ik
      nkq = qq_pair(i)%nkq
      if (nkq < 1) cycle

      !kpts in crystal coordinate
      xk = num2kpt( qg%kpt(ik), qg%ndim )
      !
      allocate( ist(nkq) )
      ist(1) = 0
      do j = 1, nkq-1
         ist(j+1) = ist(j) + qq_pair(i)%nchl(j)
      enddo
      !allocate space
      tot_chnl = sum(qq_pair(i)%nchl(:))
      allocate( bnd_idx(tot_chnl), ph3_psi2(tot_chnl))
      !allocate( qq_idx(3, tot_chnl))
      bnd_idx = 0;   ph3_psi2 = 0.0_dp;
      !qq_idx(:,:) = 0

      ! syp - phiqr(ir1,im2,im1,im): im2 for R=0, im1 for ir1, im for ik
      call phirr_to_phiqr(fct, xk, phiqr)
      
      if ( adapt_smear_phph ) then
         allocate( adsmear(tot_chnl) )
         adsmear = 0.0_dp
      endif

!$omp parallel default(shared) private(j, iq, ikq, xq, xkq, uk, uq, ukq, wq, &
!$omp& qidx, ltable, it, im, n, m, bands, vk, vkq, adsmear_tmp, ph3, expikr_set_tot, ia, ja, ka, idx, ix, ongrid, phi3_qq)
!$omp do schedule(guided)
      do j = 1, nkq
         !allocate ( expikr_set_tot%exps_storage(1, pht%na*pht%na*pht%na) )
         iq = qq_pair(i)%iq(j)
         !ikq = qq_pair(i)%ikq(j)
         ikq= kpts_plus_invert( qg%kpt(ik), qg%kpt(iq), qg%ndim) + 1
         !
         xq = num2kpt(qg%kpt(iq), qg%ndim)
         xkq = num2kpt(qg%kpt(ikq), qg%ndim)

         !
         if(use_mem) then
            uk(:,:) = uqgrid(:,:,ik)
            uq(:,:) = uqgrid(:,:,iq)
            ukq(:,:) = uqgrid(:,:,ikq)
         else
            call solve_phonon_modes(ph, xk, wq, uk)
            call solve_phonon_modes(ph, xq, wq, uq)
            call solve_phonon_modes(ph, xkq, wq, ukq)
         endif

         if (symm) then
            call match_table_absorp_emit_ph_dyna(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
            !call match_table_absorp_emit_ph(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
         else
            call match_table_absorp_emit_ph(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
            !call match_table_all(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
         endif
         !sanity check
         if( qq_pair(i)%nchl(j) .ne. count(ltable) ) &
            call errore('compute_scatter_ph3_psi2','qq_pair(i)%nchl(j) .ne. count(ltable)',1)

         ! syp - phiqr(im2,im1,im): im2 for R=0, im1 for iq, im for ik
         !            equivlently,  im2 for ikq, im1 for iq, im for ik
         call phiqr_to_phiqq(fct, xq, phiqr, phi3_qq)
         call phi3_transform(nmod, ukq, uq, uk, phi3_qq)


         it = ist(j)
         do im = 1, nmod
            if (qg%enk(im,iq) .le. wcut) cycle
         do n = 1, nmod
            if (qg%enk(n,ik) .le. wcut) cycle
         do m = 1, nmod
            if (qg%enk(m,ikq) .le. wcut) cycle

            if( ltable(m,n,im) ) then
               it = it + 1
               bands = (/m, n, im/)
               bnd_idx(it) = band2idx(bands, nmod)
               !qq_idx(:,it) = (/ik, iq, ikq/)

               if ( adapt_smear_phph ) then
                  vk = qg%vnk(:,im,iq) ! try vq first
                  vkq = qg%vnk(:,m,ikq)
                  adsmear_tmp = sqrt(2.0_dp) * calc_adaptive_smear(vk,vkq,maxval(qg%ndim))
                  ! test 
                  adsmear_tmp = min(adsmear_tmp, delta_smear_ph)
                  ! min smear = max smear * 0.05
                  adsmear_tmp = max(adsmear_tmp, delta_smear_ph*0.05_dp)
                  adsmear(it) = adsmear_tmp
               endif

               ph3 = czero;
               ! only solve for phi3mat if adpative smear is not used or 
               ! if adaptive is used and energy is less than that
               if ( adapt_smear_phph) then
                  if (abs(qg%enk(m, ikq) - qg%enk(n,ik) - qg%enk(im, iq)) > 3.0*adsmear_tmp) then
                     ph3_psi2(it) = abs2(ph3)
                     cycle
                  endif
               endif
               
               ! for # of iq and ikq
               ! disable calculating g2 in reality
           !!!!call solve_phi3mat(pht, xk, xq, xkq, qg%enk(:,ik), qg%enk(:,iq), qg%enk(:,ikq), &
           !!!!   uk, uq, ukq, n, im, m, ph3)
               
               !call solve_phi3mat(pht, expikr_set_tot, qg%enk(n,ik), qg%enk(im,iq), qg%enk(m,ikq), &
               !   uk(:,n), uq(:,im), ukq(:,m), ph3)
           !!!!ph3_psi2(it) = abs2(ph3)
               
               ph3_psi2(it) = abs2(phi3_qq(m,im,n)/sqrt(qg%enk(n,ik)*qg%enk(im,iq)*qg%enk(m,ikq))/sqrt(2.0_dp)**3/6.0_dp)


            endif
         enddo; enddo; enddo
      
      enddo
!$omp end do
!$omp end parallel
      !
      ! output to hdf5 files
      dset_name = "bands_index_" // trim( int_to_char(i) )
      call hdf_write_dataset(file_id, trim(dset_name), bnd_idx)
      dset_name = "ph3_psi2_" // trim( int_to_char(i) )
      call hdf_write_dataset(file_id, trim(dset_name), ph3_psi2)
      if ( adapt_smear_phph ) then
         dset_name = "adapt_smear_" // trim( int_to_char(i) )
         call hdf_write_dataset(file_id, trim(dset_name), adsmear)
      endif
      !dset_name = "qq_index_" // trim( int_to_char(i) )
      !call hdf_write_dataset(file_id, trim(dset_name), qq_idx)
      ! show progress
      call progressbar(i, nq_loc)
      !
      deallocate(bnd_idx, ph3_psi2, ist)
      !deallocate(qq_idx)
      if ( adapt_smear_phph ) deallocate(adsmear)
   enddo
   call hdf_close_file(file_id)
   deallocate(phi3_qq,phiqr)
   
   if( use_mem ) deallocate( uqgrid )
   
   deallocate( cryst_tau)

   call mp_barrier(inter_pool_comm)
   call stop_clock('comp_ph3psi2')
end subroutine compute_scatter_ph3_psi2

subroutine setup_qq_pair(qg, ph, wcut, e_thr, nq_loc, qq_pair, isfine_input)
   use pert_utils, only: match_table, match_table_all, match_table_absorp_emit_ph, &
      match_table_absorp_emit_ph_dyna
   use pert_param, only: symm
   use boltz_utils, only: kpts_minus, num2kpt, kpt2num, band2idx, kpts_plus_invert
   use qe_mpi_mod, only: world_comm, my_pool_id
   implicit none
   type(grid), intent(in) :: qg
   type(lattice_ifc), intent(in) :: ph
   real(dp), intent(in) :: wcut, e_thr
   !
   integer, intent(in) :: nq_loc
   type(channel_ph), intent(inout) :: qq_pair(:)
   logical, intent(in), optional :: isfine_input
   logical :: isfine
   !integer, intent(in), optional :: reduce_dims(3)
   !
   !
   !local variables
   integer :: nkpt, nmod, maxq, nkq, nq_tmp, i, ik, ikq, j, qidx, iq, k
   integer :: qidx_tmp
   logical :: ltable(qg%numb, qg%numb, ph%nm)
   real(dp) :: xq(3)
   !
   type(channel_ph), pointer :: kq_tmp(:)
   integer :: qndim(3), ikc, iqc, nkpt_c
   integer :: ierr, myrank, nprocs

   call start_clock('set_qq_pair')
   isfine = .false.
   if ( present(isfine_input)) isfine = isfine_input 
   nkpt = qg%nk
   nmod = ph%nm
   !!! TEMP CHANGE to factor of 3
   qndim = qg%ndim / 3 
   nkpt_c = qndim(1) * qndim(2) * qndim(3)

   allocate( kq_tmp(nq_loc) )
   !find iq
   maxq = qg%ndim(1)*qg%ndim(2)*qg%ndim(3)
   !(k, k') and (k', k) are connected, only one of them is stored.
   !to balance load, each kpoint is associated with equal number of k+q. 
   !so for each ik in grid%kpt, the (k+q) are choose as ik+1, ik+2, ...., ik+(nkpt-1)/2
   !  if (ik + i) is larger than nkpt, restart from 1. for example, nkpt + 5 -> 5
   !
   !initialize and pre-allocate space
   do i = 1, nq_loc
      !write(*,*) 'allocating iq for i=', i
      ik = qq_pair(i)%ik
      kq_tmp(i)%ik = ik
      !assign (k+q) evenly for each k, for load balance
      ! find coarse ikc
      if (isfine) then
         !!!! TEMP CHANGE by kyao
         !ikc = kpt2num(num2kpt(ik-1, qg%ndim), qndim) + 1
         !nkq = int((nkpt_c+1)/2) + merge(1, 0, ((mod(nkpt_c+1,2)>0) .and. ikc<=nkpt_c/2))
         nkq = nkpt_c
      elseif(symm)then
          !syp - syp's proposal: for q, we use irreducible q-points
          !                      for q', we use all q-points
          !                      for q'', obtained by momentum conservation
         nkq = nkpt
      else
         nkq = int((nkpt+1)/2) + merge(1, 0, ((mod(nkpt+1,2)>0) .and. ik<=nkpt/2))
      endif
      kq_tmp(i)%nkq = nkq

      allocate( kq_tmp(i)%iq(nkq) )
      kq_tmp(i)%iq(:) = -1
   enddo
   !
!$omp parallel do schedule(guided) default(shared) private(i, ik, nkq, j, qidx, &
!$omp& qidx_tmp, iqc, ikc)
   do i = 1, nq_loc
      ik = kq_tmp(i)%ik
      nkq = kq_tmp(i)%nkq
      !write(*,*) 'computing iq for i= ', i, ' nkq=', nkq
      !loop over (k+q) points assigned to this current k-points
      do j = 1, nkq
         if (isfine) then
            ! find coarse grid index for ik
            ikc = kpt2num(num2kpt(ik-1, qg%ndim), qndim) + 1
            iqc = mod(ikc+j-2, nkpt_c) + 1
            qidx = kpt2num(num2kpt(iqc-1, qndim), qg%ndim) + 1
         else
            qidx = mod(ik+j-2,nkpt) + 1
         endif
         !debug
         if(qidx > maxq) call errore('setup_qq_pair', 'qidx > maxq', 1)
         
         !
         if(qidx > 0) then
            kq_tmp(i)%iq(j)  = qidx
         endif
      enddo 
   enddo
!$omp end parallel do
   
   !allocate space
   do i = 1, nq_loc
      nkq = kq_tmp(i)%nkq
      !allocate( kq_tmp(i)%ikq(nkq)  )
      allocate( kq_tmp(i)%nchl(nkq) )
      !initialize
      kq_tmp(i)%nchl(:) = 0
   enddo


!$omp parallel default(shared) private(i, j, ik, iq, ikq, ltable) 
!$omp do schedule(guided)
   do i = 1, nq_loc
      !write(*,*) 'matching tables for', i
      ik = kq_tmp(i)%ik
      do j = 1, kq_tmp(i)%nkq
         iq = kq_tmp(i)%iq(j)
         !
         if(iq > 0) then
            ! generate ikq (starting from 1) from iq = -1 or starting from 1
            ! ik needs to be from 0 to 1 because this is kpts_plus, 
            ! iq also needs to be from 0 to 1
            ! invert makes ikq = -ik - iq
            ikq = kpts_plus_invert(qg%kpt(ik), qg%kpt(iq), qg%ndim) + 1
             
            !kq_tmp(i)%ikq(j) = ikq
            if (symm) then
               call match_table_absorp_emit_ph_dyna(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
               !call match_table_absorp_emit_ph(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
            else
               !call match_table_all(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
               call match_table_absorp_emit_ph(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
            endif
            kq_tmp(i)%nchl(j) = count(ltable)

         endif
      enddo
   enddo
!$omp end do
!$omp end parallel
   !write(*,*) 'collecting qq pairs'
   !collect (q,q) pair that has at least one valid scattering channel (nchl > 0).
   do i= 1, nq_loc
      nkq = count( kq_tmp(i)%nchl(:) > 0 )
      qq_pair(i)%nkq = nkq
      !write(*,*) 'new number of nkqs =', nkq
      !allocate( qq_pair(i)%iq(nkq), qq_pair(i)%ikq(nkq), qq_pair(i)%nchl(nkq) )
      if( nkq > 0) then !!TODO kellyyao - added this back
         allocate( qq_pair(i)%iq(nkq), qq_pair(i)%nchl(nkq) )
         !
         k = 0
         do j = 1, kq_tmp(i)%nkq
            if( kq_tmp(i)%nchl(j) > 0 ) then
               k = k + 1
               qq_pair(i)%iq(k) = kq_tmp(i)%iq(j)
               !qq_pair(i)%ikq(k) = kq_tmp(i)%ikq(j)
               qq_pair(i)%nchl(k) = kq_tmp(i)%nchl(j)
               !write(*,*) 'nchl for ', j, ' is ', qq_pair(i)%nchl(k)
            endif
         enddo
         !debug
         if(k .ne. nkq) call errore('setup_qq_pair', 'k .ne. nkq', 1)
      endif
      deallocate(kq_tmp(i)%iq, kq_tmp(i)%nchl)
   enddo


   deallocate( kq_tmp )
   
   call stop_clock('set_qq_pair')
end subroutine setup_qq_pair



!#include "setup_qq_pair.f90"

!#include "compute_scatter_ph3_psi2.f90"

end module boltz_scatter
