!==============================================================================
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
!  solve iterative linearized Boltzmann equation, and real-time dynamics
!
! Maintenance:
!===============================================================================

#include "cdyn_config.h"


module boltz_scatter_integral_tgt
   use kinds, only: i8b
   use pert_const, only: dp, pi
   use boltz_grid, only: grid
   use boltz_utils, only: idx2band
   use pert_utils, only: bose, fermi, gauss
   use pert_param, only: delta_smear
   use pert_output, only: stopwatch
   use yaml_utils, only: ymlout
   use boltz_scatter, only: num_kq_pairs, scatter, num_scatter_channels, scatter_channels, qgrid
   use qe_mpi_mod, only: stdout, ionode, mp_barrier, mp_sum, inter_pool_comm
   use iso_c_binding, only: c_funloc, c_funptr, c_int, c_loc, c_size_t
   use sorting, only: c_qsort
   use boltz_scatter_sizest
   implicit none
   private


   ! This is the element-type for an intermediate data structure for accumulating
   ! the indexes of sources that affect specific elements in the target array.
   ! It is a simple growable vector implementation.
   type :: t_scattgt_accum
      integer :: capacity = 0 ! Elements that we have room to store
      integer :: size     = 0 ! Elements that we are actually storing

      ! Growable array of source-array indexes that affect this target element.
      ! The sign of the index indicates add or subtract.
      integer(kind=i8b), allocatable :: sign_sources(:)
   end type t_scattgt_accum


   ! This structure represents a single target (i1, i2) element in the target
   ! array that is updated by a collection of source entries in the flat
   ! scattering array.  Source entries are represented by integer indexes.
   ! Since some updates are subtracted, the sign of the source-entry index is
   ! used to indicate addition (nonnegative) or subtraction (negative).
   type :: t_scatter_target
      ! The (n, ik) or (m, ikq) indexes packed into a single integer.
      integer :: tgt_index

      ! Number of sources in the scattgts_sources array that correspond to
      ! this scatter-target.  They start at the src_start_index and are
      ! contiguous from that point in the array.
      integer :: nsrc

      ! Starting index into the scattgts_sources table for this target index.
      ! The table can be very large, so we need 8-byte integers.  This value
      ! is after nsrc to optimize packing of the data structure.
      integer(kind=i8b) :: src_start_index

   end type t_scatter_target


   ! This is the result computed for each time-step.  It is the same number of
   ! elements as scatter_channels.
   real(dp), allocatable, save :: scatchan_mkq2nk(:)

   ! Mapping from scatter_channels elements (or scatchan_mkq2nk elements)
   ! to elements in the epcol array to be updated.
   integer, save :: num_scatter_targets
   type(t_scatter_target), allocatable, target, save :: scatter_targets(:)
   !
   integer, allocatable, save :: scattgts_sources(:)

   public :: cdyna_scatter_target_setup, cdyna_scat_int_tgt
contains

! map index to (n, ik) / (m, ikq) values
pure function idx2tgt(idx) result(tgt)
#if defined(_OPENACC)
!$acc routine seq
#endif
   implicit none
   integer, intent(in) :: idx
   integer :: tgt(2)
   tgt(1) = ibits(idx, 0,  8) + 1
   tgt(2) = ibits(idx, 8, 23) + 1
end function idx2tgt

! map (n, ik) / (m, ikq) values to index
pure function tgt2idx(n_or_m, ik_or_ikq) result(idx)
   implicit none
   integer, intent(in) :: n_or_m, ik_or_ikq
   integer :: idx
   idx = 0
   call mvbits(n_or_m    - 1, 0,  8, idx, 0)
   call mvbits(ik_or_ikq - 1, 0, 23, idx, 8)
end function tgt2idx


! Accumulate one more scatter-target and source-index into the intermediate
! scattgts_accum data structure.  This function essentially just implements
! the growable-array logic.
subroutine accum_scattgt(scattgts_accum, n, ik, src_i)
   implicit none
   type(t_scattgt_accum), intent(inout) :: scattgts_accum(:,:)
   integer, intent(in) :: n, ik
   integer(kind=i8b), intent(in) :: src_i

   integer :: sz, cap
   integer(kind=i8b), allocatable :: tmp_array(:)

   sz = scattgts_accum(n, ik)%size
   cap = scattgts_accum(n, ik)%capacity
   if (sz == cap) then
      ! We don't have capacity to store the new value.
      ! Need to allocate or reallocate
      if (cap == 0) then
         ! Need to do an initial allocation.  Start with a capacity of 16 elements.
         allocate( scattgts_accum(n, ik)%sign_sources(16) )
         scattgts_accum(n, ik)%capacity = 16
      else
         ! Need to grow the allocation.  Double it each time.
         ! In Fortran we must do a new allocation and copy the existing data over.
         ! The MOVE_ALLOC() routine is supposed to deallocate the target allocation
         ! before moving the source allocation into the target variable, but this
         ! doesn't seem to happen with nvfortran.  So, deallocate before moving.
         allocate( tmp_array(cap * 2) )
         tmp_array(1:cap) = scattgts_accum(n, ik)%sign_sources(1:cap)
         tmp_array(cap+1:) = 0
         deallocate( scattgts_accum(n, ik)%sign_sources )
         call MOVE_ALLOC(tmp_array, scattgts_accum(n, ik)%sign_sources)
         scattgts_accum(n, ik)%capacity = cap * 2
      endif
   endif

   ! Now that we know capacity > size, we store the new value
   scattgts_accum(n, ik)%size = sz + 1
   scattgts_accum(n, ik)%sign_sources(sz + 1) = src_i
end subroutine accum_scattgt


! Order the source-indexes used for a particular target-index in increasing
! absolute order, to ensure the source-array is scanned in sequence.  (We
! assume this will improve cache locality.)
integer(c_int) function compare_src_indexes(p1, p2) bind(C)
   use iso_c_binding
   implicit none
   type(c_ptr), intent(in), value :: p1, p2
   integer(kind=i8b), pointer :: src1, src2
   integer(c_int) :: v

   ! C qsort() passes the items to compare as void* pointers.
   call c_f_pointer(p1, src1)
   call c_f_pointer(p2, src2)

   ! Intermediate value is 64-bits, and we need to convert it to 32-bits.
   ! This will preserve the sign, which is all we really care about for
   ! the comparison function.
   v = rshift(abs(src1) - abs(src2), 32)

   compare_src_indexes = v
end function compare_src_indexes


integer(c_int) function compare_scatter_target_lengths_desc(p1, p2) bind(C)
   use iso_c_binding
   implicit none
   type(c_ptr), intent(in), value :: p1, p2
   type(t_scatter_target), pointer :: tgt1, tgt2
   integer(c_int) :: v

   ! C qsort() passes the items to compare as void* pointers.
   call c_f_pointer(p1, tgt1)
   call c_f_pointer(p2, tgt2)

   v = tgt2%nsrc - tgt1%nsrc
   compare_scatter_target_lengths_desc = v
end function compare_scatter_target_lengths_desc


! Generate an "inverted" scattering-channel data structure that is oriented
! around the target locations in the probability array, rather than the
! source locations.
subroutine cdyna_scatter_target_setup(kg)
   implicit none
   type(grid), intent(in) :: kg
   !local
   integer :: j, ik, iq, ikq, bands(3), m, n, im
   integer(kind=i8b) :: i, i_start, ns

   type(t_scattgt_accum), allocatable, target :: scattgts_accum(:,:)

   integer :: tgt_index, i_scat_target
   integer, parameter :: BITS_IN_BYTE = 8

!   real(dp), pointer :: kg_enk(:,:)

   call start_clock('target_setup')

   write(stdout,'(5x,a)') '>start target_setup'

   allocate( scatchan_mkq2nk(num_scatter_channels) )

   ! Populate the "scatter-targets accumulator" so we can figure out what
   ! source elements affect each target element.

   ! scattgts_accum has the same dimensions as dist0 etc.  Each element in this
   ! structure records the source elements that affect it.
   allocate( scattgts_accum(kg%numb, kg%nk) )

! Prefer SCAT_REV for this code.
#if defined(SCAT_REV)
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
#elif defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! Details for update to epcol(n, ik)
         call accum_scattgt(scattgts_accum, n, ik, ns) ! addition

         ! Details for update to epcol(m, ikq)
         call accum_scattgt(scattgts_accum, m, ikq, -ns) ! subtraction
      enddo
#if !defined(SCAT_REV)
   enddo
#endif

   ! For each scatter-target element, sort the array of source-indexes so we
   ! hopefully have good data-locality in the incremental calculation code.
   ! (We could also count up how many scatter-targets we have in this loop,
   ! but that would make the code harder for the compiler to parallelize.)
!$omp parallel do schedule(guided) collapse(2) default(private)
   do n = 1, kg%numb
      do ik = 1, kg%nk
         if (scattgts_accum(n, ik)%size == 0) cycle

         call c_qsort(c_loc(scattgts_accum(n, ik)%sign_sources(1)), &
            elem_count = int(scattgts_accum(n, ik)%size, kind=c_size_t), &
            elem_size = int(storage_size(scattgts_accum(n, ik)%sign_sources(1))/BITS_IN_BYTE, kind=c_size_t), &
            f_cmp = c_funloc(compare_src_indexes))
      enddo
   enddo

   ! Count up how many scattering targets we have.
   ! This will be at most the number of elements in epcol array.
   num_scatter_targets = 0
   do n = 1, kg%numb
      do ik = 1, kg%nk
         if (scattgts_accum(n, ik)%size > 0) then
            num_scatter_targets = num_scatter_targets + 1
         endif
      enddo
   enddo

   write (stdout, '(5X,A,I0)') ' * Number of scattering target elements:  ', num_scatter_targets

   ! Allocate the scatter_targets array, then populate its contents.
   allocate( scatter_targets(num_scatter_targets) )
   allocate( scattgts_sources(2 * num_scatter_channels) )

   i_scat_target = 0
   i_start = 1  ! Where the next run of source-indexes starts in scattgts_sources
   do n = 1, kg%numb
      do ik = 1, kg%nk
         if (scattgts_accum(n, ik)%size == 0) cycle

         ! Initialize the next scatter_targets record
         i_scat_target = i_scat_target + 1
         scatter_targets(i_scat_target)%tgt_index = tgt2idx(n, ik)
         scatter_targets(i_scat_target)%nsrc = scattgts_accum(n, ik)%size

         ! Set the source-start-index
         scatter_targets(i_scat_target)%src_start_index = i_start

         ! Copy across the source-indexes
         i = i_start + scattgts_accum(n, ik)%size
         scattgts_sources(i_start:i) = scattgts_accum(n, ik)%sign_sources
         i_start = i
      enddo
   enddo

   if (ionode) then
      ! Log data-structure allocation details into the YAML output file.
      call output_scattgts_accum_alloc_stats_to_yaml(scattgts_accum)
      call output_scatter_targets_alloc_stats_to_yaml()
   endif

   ! Deallocate this after outputting the memory usage stats,
   ! even though we could have freed the internal arrays earlier.
   do n = 1, kg%numb
      do ik = 1, kg%nk
         if (scattgts_accum(n, ik)%size > 0) then
            ! Deallocate the array of source-indexes in this scattgts_accum element
            deallocate( scattgts_accum(n, ik)%sign_sources )
         endif
      enddo
   enddo
   deallocate( scattgts_accum )

   ! EXPERIMENTAL:  Sort scatter-targets based on array lengths
#if defined(CDYN_SORT_SCAT_TGTS)
   write (stdout, '(5X,A)') ' * Sorting scatter_targets on descending subarray lengths'
   call c_qsort(c_loc(scatter_targets(1)), &
                elem_count = int(num_scatter_targets, kind=c_size_t), &
                elem_size = int(storage_size(scatter_targets(1))/BITS_IN_BYTE, kind=c_size_t), &
                f_cmp = c_funloc(compare_scatter_target_lengths_desc))
#endif

#if defined(_OPENACC)
! Tell OpenACC what data structures we want moved to the GPU,
! so it only happens once.
!$acc enter data copyin(scatter(:), scatter_channels(:)), &
!$acc&           create(scatchan_mkq2nk(:)), &
!$acc&           copyin(scatter_targets(:), scattgts_sources(:)), &
!$acc&           copyin(kg%enk(:,:), qgrid%freq(:,:))
#endif

   call stop_clock('target_setup')

   write(stdout,'(5x,a)') '>target_setup done.'
   flush(stdout)

end subroutine cdyna_scatter_target_setup


! Real time dynamics, Refer to eq.17. 18 in Eur. Phys. J. B (2016) 89: 239
! Bernardi, First-principles dynamics of electrons and phonons
subroutine cdyna_scat_int_tgt(kg, eptemp, dist, epcol)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, dist(kg%numb, kg%nk)
   ! e-ph collision term computed from given distribution function dist(:,:)
   real(dp), intent(out) :: epcol(kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, im, n, m, bands(3)
   integer(kind=i8b) :: ns
   real(dp) :: wq, enk, emkq, ph_n, fnk, fmkq, fabs, fem
#if !defined(STORE_G2DT)
   real(dp) :: qgrid_weight_2pi, g2, dt1, dt2
#endif
   integer :: tgt(2), sign_src_i
   integer(kind=i8b) :: j, j_start, j_end
   real(dp) :: total, delta

   integer :: kg_numb
   real(dp), pointer :: kg_enk(:,:), qgrid_freq(:,:)

   kg_numb = kg%numb
   kg_enk => kg%enk
   qgrid_freq => qgrid%freq

#if !defined(STORE_G2DT)
   ! Copy this value into a local scalar variable, or else OpenACC thinks the
   ! entire struct must be copied to the GPU.  It also gives us a chance to
   ! extract some math from the main loop.
   qgrid_weight_2pi = qgrid%weight * pi * 2.0_dp
#endif

#if defined(_OPENACC)
!$acc data present(scatter(:), scatter_channels(:), &
!$acc&             scatter_targets(:), scattgts_sources(:), scatchan_mkq2nk(:)) &
!$acc&     present(kg_enk(:,:), qgrid_freq(:,:)) &
!$acc&     copyout(epcol(:,:))
#endif

   ! Clear the entire target arrays.
   ! scatchan_mkq2nk(:) = 0.0_dp ! TODO(donnie):  Unnecessary?
   epcol(:,:) = 0.0_dp

#if !defined(_OPENACC)
#if !defined(STORE_G2DT)
!$omp parallel do schedule(guided) default(shared) private(i, ns, iq, ik, ikq, &
!$omp&  bands, im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fnk, fmkq, fabs, fem)
#else
!$omp parallel do schedule(guided) default(shared) private(i, ns, iq, ik, ikq, &
!$omp&  bands, im, n, m, wq, enk, emkq, ph_n, fnk, fmkq, fabs, fem)
#endif
#endif
   ! Compute the values used to update the target array
! Prefer SCAT_REV for this code.
#if defined(SCAT_REV)
#if defined(_OPENACC)
!$acc parallel loop private(bands) async
#endif
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
#elif defined(SCAT_FWD)
#if defined(_OPENACC)
!$acc parallel loop private(iq, ik, ikq, bands) async
#endif
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      ! (note:  can't seem to put an "acc loop" directive here because of the
      ! preceding operations, but I don't want to do those every iteration)
#if defined(_OPENACC)
!$acc loop seq
#endif
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg_numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! eigenvalue
         wq   = qgrid_freq(im, iq)
         enk  = kg_enk(n, ik)
         emkq = kg_enk(m, ikq)

         ! phonon occupation
         ph_n = bose(eptemp, wq) ! N_q\mu

#if !defined(STORE_G2DT)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2

         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq)
         dt1 = qgrid_weight_2pi * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid_weight_2pi * gauss(delta_smear, enk - emkq - wq)
#endif

         ! electron distributaion
         fnk  = dist(n, ik)
         fmkq = dist(m, ikq)

         ! scattering process:  formula 1
         !fabs = fnk*(one - fmkq)*wgq - fmkq*(one - fnk)*(one + wgq)
         !fem  = fnk*(one - fmkq)*(one + wgq) - fmkq*(one - fnk)*wgq
         ! formula 2
         fabs = fnk*(ph_n + fmkq) - fmkq*(ph_n + 1.0_dp)
         fem  = fnk*(ph_n + 1.0_dp) - fmkq*(ph_n + fnk)

         ! Record contribution from k+q to update array in second pass
#if !defined(STORE_G2DT)
         scatchan_mkq2nk(ns) = - g2 * (dt1*fabs + dt2*fem)
#else
         scatchan_mkq2nk(ns) = -(scatter_channels(ns)%g2_dt1 * fabs + &
                                 scatter_channels(ns)%g2_dt2 * fem)
#endif
      enddo
#if !defined(SCAT_REV)
   enddo
#endif
#if !defined(_OPENACC)
!$omp end parallel do
#endif

   ! Compute the totals for each target element in epcol array
#if defined(_OPENACC)
!$acc parallel loop private(tgt) async
#else
!$omp parallel do schedule(guided) default(shared) private(i, tgt, &
!$omp&  j, j_start, j_end, sign_src_i, total, delta)
#endif
   do i = 1, num_scatter_targets
      ! tgt_index has the (n, ik) or (m, ikq) indexes packed into
      ! a single integer.  tgt(1) = n or m, tgt(2) = ik or ikq.
      tgt = idx2tgt(scatter_targets(i)%tgt_index)
      total = 0.0_dp

      j_start = scatter_targets(i)%src_start_index
      j_end = j_start + scatter_targets(i)%nsrc - 1
#if defined(_OPENACC)
!$acc loop reduction(+:total)
#endif
      do j = j_start, j_end
         sign_src_i = scattgts_sources(j)
         if (sign_src_i < 0) then
            delta = -scatchan_mkq2nk(-sign_src_i)
         else
            delta = scatchan_mkq2nk(sign_src_i)
         endif
         total = total + delta
      enddo

      epcol(tgt(1), tgt(2)) = total
   enddo
#if !defined(_OPENACC)
!$omp end parallel do
#endif

#if defined(_OPENACC)
!$acc wait
#endif

   ! combine contributions from all the pools/processes.

   ! If OpenACC is being used, and also if CUDA-aware MPI is available, we
   ! can use the device-pointer with the MPI call.

#if defined(_OPENACC) && defined(USE_CUDA_AWARE_MPI)
   ! CUDA-aware MPI:  Need to do this before the $acc end data directive.

!$acc host_data use_device(epcol)
   call mp_sum(epcol, inter_pool_comm)
!   call mp_barrier(inter_pool_comm)
!$acc end host_data
!$acc end data
#else
   ! For regular MPI we need the $acc end data before the MPI call, to ensure
   ! the host-data has been updated properly.
!$acc end data

   call mp_sum(epcol, inter_pool_comm)
!   call mp_barrier(inter_pool_comm)
#endif

end subroutine cdyna_scat_int_tgt

!------
! Functions for reporting memory usage of data structures

! Output details regarding the memory allocation of the scatter_targets
! data structure.
subroutine output_scatter_targets_alloc_stats_to_yaml()
   implicit none
   integer(8) :: bytes_per_entry, total_bytes, subelems, total_subelems

   ! scatter_targets

   bytes_per_entry = STORAGE_SIZE(scatter_targets(1)) / 8
   total_bytes = bytes_per_entry * SIZE(scatter_targets, kind=8)

   write(ymlout,'(3x, a)') 'scatter_targets:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_scatter_targets
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a/)') 'allocations: 1'

   call boltzscat_accum_cpugpu_bytes(total_bytes)

   ! scattgts_sources

   bytes_per_entry = STORAGE_SIZE(scattgts_sources(1)) / 8
   total_bytes = bytes_per_entry * SIZE(scattgts_sources, kind=8)

   write(ymlout,'(3x, a)') 'scattgts_sources:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', 2 * num_scatter_channels
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a/)') 'allocations: 1'

   call boltzscat_accum_cpugpu_bytes(total_bytes)
end subroutine output_scatter_targets_alloc_stats_to_yaml

! Output details regarding the memory allocation of the scattgts_accum
! intermediate data structure.
subroutine output_scattgts_accum_alloc_stats_to_yaml(scattgts_accum)
   implicit none
   type(t_scattgt_accum), intent(in) :: scattgts_accum(:,:)
   integer(8) :: entries, bytes_per_entry, total_bytes
   integer :: allocations, n, ik

   ! scattgts_accum

   entries = SIZE(scattgts_accum, kind=8)
   bytes_per_entry = STORAGE_SIZE(scattgts_accum(1, 1)) / 8
   total_bytes = bytes_per_entry * entries

   ! Sub-arrays
   allocations = 1  ! Include top-level allocation
   do n = 1, SIZE(scattgts_accum, 1)
      do ik = 1, SIZE(scattgts_accum, 2)
         if (scattgts_accum(n, ik)%size == 0) cycle

         allocations = allocations + 1
         total_bytes = total_bytes + &
            SIZE(scattgts_accum(n, ik)%sign_sources, kind=8) * &
            STORAGE_SIZE(scattgts_accum(n, ik)%sign_sources(1)) / 8
      enddo
   enddo

   write(ymlout,'(3x, a)') 'scattgts_accum [intermediate]:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', entries
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i8/)') 'allocations: ', allocations

   call boltzscat_accum_cpu_bytes(total_bytes)
end subroutine output_scattgts_accum_alloc_stats_to_yaml

end module boltz_scatter_integral_tgt
