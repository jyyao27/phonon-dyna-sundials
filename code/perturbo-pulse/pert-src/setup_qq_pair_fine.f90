!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
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
subroutine setup_qq_pair_fine(qg, ph, num_qq_pair, iq_range, ph3_scat, total_scat, qq_pair_arr, q_start_list)
   use pert_utils, only: match_table, match_table_all
   use qe_mpi_mod, only: stdout, mp_barrier, inter_pool_comm, npool, my_pool_id, &
      mp_sum, ionode_id, mp_gather, mp_split_pools, ionode, mp_bcast
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus_invert, &
      find_q_neighbor, kpt2num
   implicit none
   type(grid), intent(in) :: qg
   type(lattice_ifc), intent(in) :: ph
   integer, intent(in) :: num_qq_pair, total_scat
   integer, intent(in) :: iq_range(2)
   type(t_scatt_ph), intent(inout) :: ph3_scat(num_qq_pair)
   integer, intent(in) :: qq_pair_arr(2,total_scat)
   integer, intent(in) :: q_start_list(:)
   !
   !
   !local variables
   integer :: nkpt, nmod, maxq, nkq, nq_tmp, i, ik, ikq, j, qidx, iq, k
   integer :: is, qs, qe, ongrid, ix, ongrid2 
   integer :: iq1, iq2, iq_nb1, iq_nb2, inb1, inb2
   integer, allocatable :: iq_nb1a(:), iq_nb2a(:)
   integer, allocatable :: iqf_tmp(:,:)
   integer :: iqf_tmp_size
   logical :: ltable(qg%numb, qg%numb, ph%nm)
   real(dp) :: xq(3), qpt1(3), qpt2(3), fweights(3)
   integer :: qndim(3)
   integer :: qc_eqv1, qc_eqv2, nkpt_c, ic, icr, ic0
   !
   integer, allocatable :: iqlist(:,:)
   logical :: found, found_reverse
   integer :: icp(4)

   call start_clock('set_qq_pair_fine')
   
   do i = 1, 3
      fweights(i) = 1.0_dp / qg%ndim(i)
      qndim(i) = qg%ndim(i)/2
   enddo
   nkpt = qg%nk
   nkpt_c = qndim(1) * qndim(2) * qndim(3)
   nmod = ph%nm
   !find iq
   !
   ic = 0
   icr = 0
   ic0 = 0
   write(*,*) my_pool_id, num_qq_pair, total_scat
   write(*,*) my_pool_id, ph3_scat(1)%ik, iq_range(1)
   write(*,*) size(q_start_list)
   do i = 1, size(q_start_list)
      write(*,*) i, q_start_list(i)
   enddo
   do i = 1, total_scat
      if (mod(i, 1000) .eq. 1) write(*,*) i
      if (allocated(iq_nb1a)) deallocate(iq_nb1a)
      if (allocated(iq_nb2a)) deallocate(iq_nb2a)
      iq1 = qq_pair_arr(1,i)
      iq2 = qq_pair_arr(2,i)
      ! if on grid, no use for neighbors
      qpt1 = num2kpt(iq1-1, qg%ndim)
      qpt2 = num2kpt(iq2-1, qg%ndim)
      ongrid = 0
      do ix = 1, 3
         if (mod(abs(nint(qpt1(ix)/fweights(ix))),2) == 1) then 
            ongrid = ongrid + 1
         endif
      enddo
      do ix = 1, 3
         if (mod(abs(nint(qpt2(ix)/fweights(ix))),2) == 1) then 
            ongrid = ongrid + 1
         endif
      enddo
      if ( ongrid == 0) then
         found = .true.
         ic0 = ic0 + 1
         cycle
      endif
      call find_q_neighbor(qpt1, qg%ndim, fweights, iq_nb1a)
      call find_q_neighbor(qpt2, qg%ndim, fweights, iq_nb2a)
      found = .false.
      found_reverse = .false.
      ! now this is a loop, need to change to random
      do inb1 = 1, size(iq_nb1a)
         iq_nb1 = iq_nb1a(inb1)
         ! convert to coarse grid idx
         qc_eqv1 = kpt2num(num2kpt(iq_nb1-1, qg%ndim), qndim) + 1
         ! if our pool's num_qq_pair doesnt even have iq_nb1, we just skip
         !do is = 1, num_qq_pair
         do inb2 = 1, size(iq_nb2a)
            iq_nb2 = iq_nb2a(inb2)
            qc_eqv2 = kpt2num(num2kpt(iq_nb2-1, qg%ndim), qndim) + 1
            ! search for iq_nb1, iq_nb2 correspoding is from ph3_scat
            !if ((iq_range(1) > qc_eqv1) .or. (iq_range(2) < qc_eqv1)) then
               !if ((iq_range(1) > qc_eqv2) .or. (iq_range(2) < qc_eqv2)) then
                  !write(*,*) qc_eqv1, qc_eqv2, 'cycled'
                  !cycle
               !endif
            !endif
            !if ((iq_range(1) <= qc_eqv1) .and. (iq_range(2) >= qc_eqv1)) then
               !qs = q_start_list(qc_eqv1)
               !if (qs <= num_qq_pair) then
                  !if (ph3_scat(qs)%ik .eq. iq_nb1) then
                     !qe = merge(q_start_list(qc_eqv1+1), num_qq_pair, qc_eqv1 < nkpt_c)
                     do is = 1, num_qq_pair
                     !do is = qs, qe 
                        if (ph3_scat(is)%iq == iq_nb2) then
                        if (ph3_scat(is)%ik == iq_nb1) then
                           found = .true.
                           ic = ic + 1
                           exit
                        endif; endif
                     enddo
                  !endif
               !else
               !   write(*,*) qc_eqv1, qs, qe, num_qq_pair
               !endif
            !endif
            if ((found .eqv. .false.) .and. (iq_nb1 .ne. iq_nb2)) then
               if ((iq_range(1) <= qc_eqv2) .and. (iq_range(2) >= qc_eqv2)) then
                  qs = q_start_list(qc_eqv2)
                  !if (qs <= num_qq_pair) then
                  !write(*,*) num_qq_pair, size(ph3_scat), iq_nb2
                     if (ph3_scat(qs)%ik .eq. iq_nb2) then
                        qe = merge(q_start_list(qc_eqv2+1), num_qq_pair, qc_eqv2 < nkpt_c)
                        !do is = 1, num_qq_pair
                        do is = qs, qe
                           if (ph3_scat(is)%iq == iq_nb1) then
                           !if (ph3_scat(is)%ik == iq_nb2) then
                              found_reverse = .true.
                              icr = icr + 1
                              exit
                           endif!; endif
                        enddo
                     endif
                  !else
                  !   write(*,*) qc_eqv1, qc_eqv2, qs, qe, num_qq_pair
                  !endif
               endif
            endif
            if (( found ) .or. (found_reverse )) then
               exit
            endif
         enddo;
         !if ( found .eqv. .true.) then
         if (( found ) .or. (found_reverse )) then
            exit
         endif
      enddo;
      if (( found ) .or. (found_reverse )) then
      !if ( found .eqv. .true.) then
         ! check if ph3_scat(is)%iqf is allocated
         !write(*,'(3(f7.3,",",2x))') num2kpt(iq_nb1+1, qg%ndim)
         !write(*,'(3(f7.3,",",2x))') num2kpt(iq_nb2+1, qg%ndim)
         !do ik = 1, 64
         !   if (ph3_scat(is)%iqf(1,ik) .eq. 0) then
         !      ph3_scat(is)%iqf(1,ik) = iq1
         !      ph3_scat(is)%iqf(2,ik) = iq2
         !   endif
         !enddo
         !if ( associated(ph3_scat(is)%iqf) ) then
         !   ! read (2,n)
         !   iqf_tmp_size = size(ph3_scat(is)%iqf, 2)
         !   allocate(iqf_tmp(2,iqf_tmp_size))
         !   iqf_tmp = ph3_scat(is)%iqf(:,:)
         !   deallocate(ph3_scat(is)%iqf)
         !   allocate(ph3_scat(is)%iqf(2,iqf_tmp_size+1))
         !   ph3_scat(is)%iqf(:,1:iqf_tmp_size) = iqf_tmp(:,:) 
         !   ph3_scat(is)%iqf(:,iqf_tmp_size+1) = (/iq1, iq2/)
         !   if ( found_reverse ) then
         !      ph3_scat(is)%iqf(:,iqf_tmp_size+1) = (/iq2, iq1/)
         !   endif
         !   deallocate(iqf_tmp)
         !else
         !   allocate(ph3_scat(is)%iqf(2,1))
         !   ph3_scat(is)%iqf(:,1) = (/iq1, iq2/)
         !   if ( found_reverse ) then
         !      ph3_scat(is)%iqf(:,1) = (/iq2, iq1/)
         !   endif
         !endif
      endif
      ! what to do if nothing is found - no worries
      ! record is, we also need a 
   enddo
   icp(my_pool_id+1) = ic
   call mp_sum(ic, inter_pool_comm)
   call mp_sum(icr, inter_pool_comm)
   call mp_sum(icp, inter_pool_comm)
   write(*,*) 'icp=', icp
   call mp_barrier(inter_pool_comm)
   write(*,*) 'ic = ', ic0, ic, icr, icr+ic+ic0
   call stop_clock('set_qq_pair_fine')
end subroutine setup_qq_pair_fine
