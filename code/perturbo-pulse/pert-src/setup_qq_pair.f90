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

subroutine setup_qq_pair(qg, ph, wcut, e_thr, nq_loc, qq_pair, isfine_input)
   use pert_utils, only: match_table, match_table_all, match_table_absorp_emit_ph
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
   qndim = qg%ndim / 2 
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
         ikc = kpt2num(num2kpt(ik-1, qg%ndim), qndim) + 1
         nkq = int((nkpt_c+1)/2) + merge(1, 0, ((mod(nkpt_c+1,2)>0) .and. ikc<=nkpt_c/2))
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
               call match_table_absorp_emit_ph(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
            else
               call match_table_all(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
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
      !deallocate(kq_tmp(i)%iq, kq_tmp(i)%ikq, kq_tmp(i)%nchl)
   enddo


   deallocate( kq_tmp )
   

   call stop_clock('set_qq_pair')
end subroutine setup_qq_pair
