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

module boltz_scatter
   use boltz_grid, only: grid
   use pert_const, only: dp, czero
   use pert_data,  only: epr_fid, qc_dim, kc_dim
   use pert_param, only: phfreq_cutoff, phfreq_cutoff_ph, delta_smear, delta_smear_ph, use_mem, &
      load_scatter_eph, load_scatter_ph3, dyna_phph
   use band_structure, only: electron_wann, solve_eigenvalue_vector
   use phonon_dispersion, only: lattice_ifc, lattice_ifct, solve_phonon_modes
   use parallel_include
   use qe_mpi_mod, only: stdout, mp_barrier, inter_pool_comm, npool, my_pool_id, &
      mp_sum, ionode_id, mp_gather, mp_split_pools, mp_bcast
   use elphon_coupling_matrix, only: elph_mat_wann, release_elph_mat_wann, &
      init_elph_mat_wann, eph_fourier_el_para, eph_fourier_elph, eph_transform_fast
   use pert_utils, only: find_free_unit
   use yaml_utils, only: ymlout, output_scatter_info_yaml
   implicit none
   public
   !user-defined data strucutre for scatering channel. based on kgrid
   type :: t_scatt
      !location of k and kq in kgrid%kpt: xk(1:3) = num2kpt(kgrid%kpt(ik), ndim)
      !NOTE: (k, kq) and (kq, k) are equivalent, only one of them are stored.
      integer :: ik
      integer :: iq
      integer :: ikq
      !number of (mkq, nk, mu) terms, mu is phonon modes index
      !only g_m,n,mu of (k, kq) satisfy the following condition are stored. 
      !abs( abs(e_nk - e_mkq) - w_qmu ) < delta_E. (delta_E = e_thr*broadening)
      integer :: nchl
      
      !use logical table to determine nt and select (mkq, nk, mu) terms:
      ! if table_bands(mkq, nk, mu) is true, then
      !   matrix elements:  |g_{ mu, iq}(nk,ik; mkq,ikq)|^2 is stored
      !use any() or count() to check if there is .true. in table_bands and nt
      !
      !index and g2 are stored for selected (mkq, nk, mu). 
      ! map (mkq, nk, mu) to a integer 'index' and store it in bands_index
      ! (mu-1)*numb*numb + (nk-1)*numb + (mkq-1)  + 1 = index
      integer, pointer :: bands_index(:)=>null()
      !integer, pointer :: iqf(:,:)=>null()
      !eph_g2(nchl), nchl = size( bands_index ).
      real(dp), pointer :: eph_g2(:)=>null()
      !eph_g2(it) : |g_{ mu, iq}(nk,ik; mkq,ikq)|^2
   end type t_scatt

   !user-defined data strucutre for scatering channel. based on kgrid
   type :: t_scatt_ph
      !location of k and kq in kgrid%kpt: xk(1:3) = num2kpt(kgrid%kpt(ik), ndim)
      !NOTE: (k, kq) and (kq, k) are equivalent, only one of them are stored.
      integer :: ik
      integer :: iq
      !number of (mkq, nk, mu) terms, mu is phonon modes index
      !only g_m,n,mu of (k, kq) satisfy the following condition are stored. 
      !abs( abs(e_nk - e_mkq) - w_qmu ) < delta_E. (delta_E = e_thr*broadening)
      integer :: nchl
      
      !use logical table to determine nt and select (mkq, nk, mu) terms:
      ! if table_bands(mkq, nk, mu) is true, then
      !   matrix elements:  |g_{ mu, iq}(nk,ik; mkq,ikq)|^2 is stored
      !use any() or count() to check if there is .true. in table_bands and nt
      !
      !index and g2 are stored for selected (mkq, nk, mu). 
      ! map (mkq, nk, mu) to a integer 'index' and store it in bands_index
      ! (mu-1)*numb*numb + (nk-1)*numb + (mkq-1)  + 1 = index
      integer, pointer :: bands_index(:)=>null()
      !integer, pointer :: iqf(:,:)=>null()
      !eph_g2(nchl), nchl = size( bands_index ).
      real(dp), pointer :: eph_g2(:)=>null()
      !eph_g2(it) : |g_{ mu, iq}(nk,ik; mkq,ikq)|^2
   end type t_scatt_ph

   type :: t_scatt_smear
      integer :: nchl
      
      real(dp), pointer :: adsmear(:)=>null()
      ! adsmear(it) : adaptive gaussian smearing
   end type t_scatt_smear

   type, private :: channel_ph
      integer :: ik
      integer :: nkq 
      integer, pointer :: iq(:) => null()
      integer, pointer :: nchl(:) => null()
   end type

   type, private :: channel
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
      real(dp) :: dscat
   end type

   integer, private, save :: nk_loc, nq_loc, nq_loc_fine, iq_st1, iq_end1
   type(channel), allocatable, private, save :: kq_pair(:)

   type(channel_ph), allocatable, private, save :: qq_pair(:), qq_pair_fine(:)
   integer, allocatable, save :: qq_pair_arr(:,:)
   integer, allocatable, save :: q_start_list(:)
   !type(qq_interp_list), allocatable, private, save :: qq_interp(:)
   !
   type(phonon_grid), save :: qgrid
   !
   integer(kind=8), save :: num_kq_pair, num_qq_pair, num_scat, ph3_num_scat
   !integer, save :: num_kq_pair, num_qq_pair, num_scat, ph3_num_scat
   integer(kind=8), save :: num_qq_pair_total
   type(t_scatt), allocatable, save :: scat(:)
   type(t_scatt_ph), allocatable, save :: ph3_scat(:)
   type(t_scatt_smear), allocatable, save :: tsmear_ph(:)

   integer, allocatable, save :: scatter_ph(:,:)
   integer, allocatable, save :: scatter_e(:,:)

   type(t_scatt_compact), allocatable, save :: tscat_ph(:)
   type(t_scatt_compact), allocatable, save :: tscat(:)
   real(dp), allocatable, save :: tsmear_ph_list(:)
   real(dp), allocatable, save :: tsmear_list(:)
   
contains

subroutine boltz_scatter_setup(kg, el, ph, qdim)
   use yaml_utils, only: ymlout, output_scatter_info_yaml
   implicit none
   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   !dimension for the q-grid. q-should be commensurate with the k-grid
   integer, intent(in) :: qdim(3)
   !local 
   integer :: nkpt, i, ik, collect_pools(npool), nq_col(npool)
   real(dp) :: wcut, e_thr
   type(elph_mat_wann) :: elph
   
   call start_clock('scat_setup')
   write(stdout,'(5x,a20)') '>start scatter_setup'
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
   num_kq_pair = 0
   num_scat = 0
   do i = 1, nk_loc
      if (kq_pair(i)%nkq < 1) cycle
      num_kq_pair = num_kq_pair + kq_pair(i)%nkq
      num_scat = num_scat + sum( kq_pair(i)%nchl(:) )
   enddo
   !
   nq_col = 0
   nq_col( my_pool_id+1) = qdim(1)*qdim(2)*qdim(3) 

   call mp_sum( nq_col, inter_pool_comm )
   !
   collect_pools = 0
   collect_pools( my_pool_id+1 ) = num_kq_pair
   call mp_sum( collect_pools, inter_pool_comm )
   !output info
   write(stdout, '(5x, a)') "No. of q-points and (k,q) pairs on each processor:"
   do i = 1, npool
      write(stdout, '(5x, a, i4, a, i9, a, i12)')  &
         "Proc.", i, ':  #.q ', nq_col(i), ' |  #.(k,q)', collect_pools(i)
   enddo

   ! output to the YAML file
   call output_scatter_info_yaml(npool, nq_col, collect_pools)
   !
   ! Log data-structure allocation details into the YAML output file.
   write(ymlout, '(/,a)') 'boltz_scatter memory:'
   call output_kq_pair_alloc_stats_to_yaml()
   ! compute and save g2 to disk
   if(.not. load_scatter_eph) then
      call init_elph_mat_wann(epr_fid, kc_dim, qc_dim, elph)
      !compute g2 and sove them to temperary files
      call compute_scatter_eph_g2 &
         (kg, el, ph, elph, qgrid, wcut, e_thr, use_mem)
      !release space
      call release_elph_mat_wann(elph)
   endif
   !
   !read g2 from files
   call load_scatter_channels(kg, ph%nm)
   call output_scatter_e_alloc_stats_to_yaml()
   call output_tscat_alloc_stats_to_yaml()
   !call load_scatter_eph_g2(nk_loc, kq_pair, num_kq_pair, scat, tsmear_e, kg, ph%nm)
   !call output_scat_alloc_stats_to_yaml()
   !do i = 1, nk_loc
   !   if( kq_pair(i)%nkq == 0 ) cycle
   !   deallocate( kq_pair(i)%iq, kq_pair(i)%ikq, kq_pair(i)%nchl )
   !enddo
   if ( dyna_phph ) then
      ! qgrid bascially not needed, deallocate to save space
      deallocate(qgrid%list)
      deallocate(qgrid%freq)
   endif


   !
   write(stdout,'(5x,a20)') '>scatter_setup done.'
   call stop_clock('scat_setup')
end subroutine boltz_scatter_setup

subroutine boltz_ph3scat_setup(qg, ph, pht, isfine_in)
   use boltz_utils, only : num2kpt, kpt2num
   implicit none
   type(grid), intent(in) :: qg
   type(lattice_ifc), intent(in) :: ph
   type(lattice_ifct), intent(in) :: pht
   logical, intent(in), optional :: isfine_in
   real(dp) :: wcut, e_thr
   integer :: nkpt, i, j, ik, collect_pools(npool), nq_col(npool)
   integer(kind=8) :: itotal
   !dimension for the q-grid. q-should be commensurate with the k-grid
   integer, allocatable :: recvcount(:), displs(:)
   integer :: nword, iq_st, iq_end
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
   else
      nkpt_c = nkpt
   endif
   wcut  = phfreq_cutoff_ph
   e_thr = delta_smear_ph*3.0_dp  !exp(-3*3) ~ 10^-4

   
   call mp_split_pools(nkpt_c, iq_st, iq_end, nq_loc)
   !interleave distribution of qpoints over pools, for better load balance
   !nq_loc = nkpt/npool + merge(1, 0, my_pool_id < mod(nkpt, npool))
   allocate( qq_pair(nq_loc) )
   i = 0
   do ik = iq_st, iq_end
      i = i + 1
      if (isfine) then
         qq_pair(i)%ik = kpt2num(num2kpt(ik-1, qndim),qg%ndim) + 1
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
   collect_pools = 0
   collect_pools( my_pool_id+1 ) = num_qq_pair
   call mp_sum( collect_pools, inter_pool_comm )
   !output info
   write(stdout, '(5x, a)') "No. of q-points and (q,q) pairs on each processor:"
   num_qq_pair_total = 0
   do i = 1, npool
      num_qq_pair_total = num_qq_pair_total + collect_pools(i) 
      write(stdout, '(5x, a, i4, a, i12)')  &
         "Proc.", i, ' |  #.(q,q)', collect_pools(i)
   enddo
   ! Log data-structure allocation details into the YAML output file.
   write(ymlout, '(/,a)') 'boltz_scatter memory:'
   call output_qq_pair_alloc_stats_to_yaml()

   if(.not. load_scatter_ph3) then
      !compute psi2 and sove them to temperary files
      call compute_scatter_ph3_psi2(qg, ph, pht, wcut, e_thr, use_mem)
      !release space
      !call release_elph_mat_wann(elph)
   endif

   ! setup 2d array for scatter_ph
   call load_scatter_channels_ph()
   call output_scatter_ph_alloc_stats_to_yaml()
   call output_tscat_ph_alloc_stats_to_yaml()

   !else ! no saving space with array
      !allocate( ph3_scat(num_qq_pair) )
      !call load_scatter_ph3_psi2(nq_loc, qq_pair, num_qq_pair, ph3_scat, tsmear_ph)
      !call output_ph3_scat_alloc_stats_to_yaml()
      !do i = 1, nq_loc
      !   if( qq_pair(i)%nkq == 0 ) cycle
      !   deallocate( qq_pair(i)%iq, qq_pair(i)%nchl )
      !enddo
      !deallocate( qq_pair)
   !endif


   ! broadcast 2d array for scatter_ph
   !if (isfine) then
   !   allocate( scatter_ph(2, num_qq_pair_total) )
   !   ! gather q_start_list
   !   allocate( recvcount(npool), displs(npool) )
   !   nword = 2
   !   displs = 0
   !   recvcount = 0
   !   recvcount( my_pool_id+1 ) = num_qq_pair * nword
   !   do i = 2, npool
   !      displs(i) = displs(i - 1) + recvcount(i - 1) 
   !   enddo 
   !   call MPI_AllGATHERV(scatter_ph_loc, recvcount(my_pool_id+1), MPI_INTEGER, &
   !      scatter_ph, recvcount, displs, MPI_INTEGER, inter_pool_comm, ierr)
   !   deallocate( recvcount, displs)
   !   call output_scatter_ph_alloc_stats_to_yaml()
   !   !deallocate( recvcount, displs, scatter_ph_loc)
   !endif

   write(stdout,'(5x,a20)') '>ph3scat_setup done.'
   call stop_clock('phscat_setup')
end subroutine boltz_ph3scat_setup

subroutine boltz_match_grid(qg, ph)
   implicit none
   type(grid), intent(in) :: qg
   type(lattice_ifc), intent(in) :: ph
   real(dp) :: wcut, e_thr
   integer :: nkpt, i, ik, iq, collect_pools(npool)
end subroutine boltz_match_grid

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

subroutine load_scatter_channels(kg, nmodes)
   use hdf5_utils
   use boltz_utils, only: idx2band
   use pert_param, only: prefix, tmp_dir, adapt_smear_eph, ph_mode_exclude_ranges
   implicit none

   type(grid), intent(in) :: kg
   integer, intent(in) :: nmodes
   integer :: i_scat, i_scat_chl, i, ik, nkq, j, tot_chnl, nchl, it, ib
   integer, allocatable :: bnd_idx(:)
   real(dp), allocatable :: eph_g2(:)
   real(dp), allocatable :: adsmear(:)
   !
   integer(HID_T) :: file_id
   ! Exclude ph mode variables
   integer :: im, jm, im1, im2, iscat, ns
   integer :: bands(3)
   logical :: include_phonon_modes(nmodes)
   
   call start_clock("load_eph_g2")
   call open_pool_eph_g2_file(file_id)

   if( any(ph_mode_exclude_ranges(:, :) > 0) .or. (.not. dyna_phph)) then
      allocate(scatter_e(3, num_kq_pair))
   else
      allocate(scatter_e(2, num_kq_pair))
   endif

   allocate(tscat(num_scat))

   if ( adapt_smear_eph) then
      allocate(tsmear_list(num_scat))
   endif
   i_scat_chl = 0
   i_scat = 0
   ! for each ik
   do i = 1, nk_loc
      nkq = kq_pair(i)%nkq
      if (nkq < 1) cycle
      ik = kq_pair(i)%ik

      tot_chnl = sum( kq_pair(i)%nchl(:) )
      allocate( bnd_idx(tot_chnl), eph_g2(tot_chnl) )
      bnd_idx = 0;  eph_g2 = 0.0_dp
      if ( adapt_smear_eph ) then
         allocate( adsmear(tot_chnl) )
         adsmear = 0.0_dp
         call read_pool_eph_g2_data(file_id, i, bnd_idx, eph_g2, adsmear)
      else
         call read_pool_eph_g2_data(file_id, i, bnd_idx, eph_g2)
      endif

      it = 0
      ! for each iq
      do j = 1, nkq
         i_scat = i_scat + 1
         nchl = kq_pair(i)%nchl(j)
         scatter_e(1,i_scat) = ik
         scatter_e(2,i_scat) = kq_pair(i)%ikq(j)
         if( any(ph_mode_exclude_ranges(:, :) > 0) .or. (.not. dyna_phph)) then
            scatter_e(3,i_scat) = kq_pair(i)%iq(j)
         endif
         ! for each channel
         do ib = 1, nchl
            i_scat_chl = i_scat_chl + 1
            it = it + 1
            !write(*,*) ik, it, i_scat, i_scat_chl 
            ! fill scat_index, band_index, and matrix element g
            tscat(i_scat_chl)%scat_index = i_scat
            tscat(i_scat_chl)%bands_index = bnd_idx(it) 
            tscat(i_scat_chl)%g2 = eph_g2(it) 
      
            if (adapt_smear_eph ) then
               tsmear_list(i_scat_chl) = adsmear(it)
            endif
         enddo
      enddo
      ! deallocate right after 
      deallocate( kq_pair(i)%ikq, kq_pair(i)%nchl )
      deallocate( kq_pair(i)%iq)
      deallocate(bnd_idx, eph_g2)
      if ( adapt_smear_eph ) deallocate( adsmear )
      !
   enddo
   deallocate( kq_pair )
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
!$omp parallel do schedule(guided) default(shared) private(i_scat_chl, bands, im)
      do i_scat_chl = 1, num_scat
         bands = idx2band( tscat(i_scat_chl)%bands_index, kg%numb )
         im = bands(3)
         ! Set the g2 elements to zero for the phonon modes to exclude
         if (.not. include_phonon_modes(im)) then
            tscat(i_scat_chl)%g2 = 0.0_dp
         endif
      enddo
!$omp end parallel do

      call stop_clock("exclude_ph_modes")
   end if

end subroutine load_scatter_channels

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

   integer :: i_scat, i_scat_chl, i, ik, nkq, j, tot_chnl, nchl, it, ib
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
   if(.not. has_file) call errore('load_scatter_ph3_psi2', "Missing file: "//trim(fname), 1)
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
   write(ymlout,'(6x, a, i8/)') 'allocations: ', total_allocs
end subroutine output_kq_pair_alloc_stats_to_yaml

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
subroutine output_ph3_scat_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! scat

   bytes_per_entry = STORAGE_SIZE(ph3_scat(1)) / 8
   total_bytes = bytes_per_entry * SIZE(ph3_scat, kind=8)

   ! subarrays of scat
   total_allocs = 1
   do i = 1, num_qq_pair
      if (ph3_scat(i)%nchl == 0) cycle

      ! Each element in scat has two subarrays
      total_allocs = total_allocs + 2
      total_bytes = total_bytes + &
         SIZE(ph3_scat(i)%bands_index, kind=8) * STORAGE_SIZE(ph3_scat(i)%bands_index(1)) / 8 + &
         SIZE(ph3_scat(i)%eph_g2     , kind=8) * STORAGE_SIZE(ph3_scat(i)%eph_g2(1)     ) / 8
   enddo

   write(ymlout,'(3x, a)') 'ph3_scat:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_qq_pair
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i8/)') 'allocations: ', total_allocs
end subroutine output_ph3_scat_alloc_stats_to_yaml

! Output details regarding the memory allocation of the scatter_ph
! data structure.
subroutine output_scatter_e_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! scat

   bytes_per_entry = (STORAGE_SIZE(scatter_e(1,1)) + STORAGE_SIZE(scatter_e(2,1))) / 8
   total_bytes = bytes_per_entry * SIZE(scatter_e, kind=8)

   write(ymlout,'(3x, a)') 'scatter_e:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_kq_pair
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
end subroutine output_scatter_e_alloc_stats_to_yaml

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
subroutine output_tscat_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! scat

   bytes_per_entry = STORAGE_SIZE(tscat(1)) / 8
   total_bytes = bytes_per_entry * SIZE(tscat, kind=8)

   write(ymlout,'(3x, a)') 'tscat:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_scat
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
end subroutine output_tscat_alloc_stats_to_yaml

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


! Output details regarding the memory allocation of the scat
! data structure.
subroutine output_scat_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! scat

   bytes_per_entry = STORAGE_SIZE(scat(1)) / 8
   total_bytes = bytes_per_entry * SIZE(scat, kind=8)

   ! subarrays of scat
   total_allocs = 1
   do i = 1, num_kq_pair
      if (scat(i)%nchl == 0) cycle

      ! Each element in scat has two subarrays
      total_allocs = total_allocs + 2
      total_bytes = total_bytes + &
         SIZE(scat(i)%bands_index, kind=8) * STORAGE_SIZE(scat(i)%bands_index(1)) / 8 + &
         SIZE(scat(i)%eph_g2     , kind=8) * STORAGE_SIZE(scat(i)%eph_g2(1)     ) / 8
   enddo

   write(ymlout,'(3x, a)') 'scat:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_kq_pair
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i8/)') 'allocations: ', total_allocs
end subroutine output_scat_alloc_stats_to_yaml

#include "setup_kq_pair.f90"

#include "compute_scatter_eph_g2.f90"

#include "load_scatter_eph_g2.f90"

#include "setup_qq_pair.f90"
#include "setup_qq_pair_fine.f90"

#include "compute_scatter_ph3_psi2.f90"

end module boltz_scatter
