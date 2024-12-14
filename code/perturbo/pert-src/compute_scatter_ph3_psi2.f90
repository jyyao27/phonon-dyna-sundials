subroutine compute_scatter_ph3_psi2(qg, ph, pht, wcut, e_thr, use_mem)
   use pert_param, only: prefix, tmp_dir, delta_smear_ph, adapt_smear_phph
   use pert_utils, only: abs2, match_table, match_table_all, get_exp_ikr, pure_exp_ikr
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus, &
      calc_adaptive_smear, kpts_plus_invert
   use pert_output, only: progressbar_init, progressbar
   use phonon_dispersion, only: solve_phi3mat, solve_phonon_velocity_findiff
   use qe_mpi_mod, only: stdout, ionode, world_comm, mp_barrier
   use force_constant_thirdorder, only: expikr_set, expikr
   use triplet_cell, only : trip_cell
   !
   !use pert constant
   use pert_const, only: ryd2mev
   use hdf5_utils
   use pert_data, only: nat, tau
   implicit none
   integer :: ierr, myrank, nprocs
   type(grid), intent(in) :: qg
   type(lattice_ifc), intent(in) :: ph
   type(lattice_ifct), intent(in) :: pht
   type(trip_cell), pointer :: tc
   !
   logical, intent(in) :: use_mem
   real(dp), intent(in) :: wcut, e_thr
   !
   type(expikr_set) :: expikr_set1, expikr_set2, expikr_set_tot

   integer :: nmod, qidx, numq1
   integer :: i, j, ik, iq, ikq, it, im, n, m, nkq, tot_chnl, bands(3)
   integer :: ia, ja, ka, iexp, idx
   real(dp) :: xq(3), xk(3), xkq(3), wq(ph%nm), wk(ph%nm), wkq(ph%nm)
   complex(dp) :: uk(ph%nm, ph%nm), uq(ph%nm, ph%nm), ukq(ph%nm, ph%nm),ph3
   logical :: ltable(ph%nm, ph%nm, ph%nm)
   !
   integer, allocatable :: ist(:), bnd_idx(:)!, qq_idx(:,:)
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
  
   ! initialize array exp_ikr for all q's
   !allocate(expikr_set1%exps_storage(numq1,pht%na*pht%na*pht%na)) 
   !allocate(expikr_set2%exps_storage(numq1,pht%na*pht%na*pht%na)) 
   !idx = 0
   !do ia = 1, pht%na 
   !   do ja = 1, pht%na 
   !      do ka = 1, pht%na
   !         idx = idx+1
   !         tc => pht%phit(idx)%trip_ph
   !         do iq = 1, numq1
   !            allocate(expikr_set1%exps_storage(iq,idx)%exps(tc%nrt)) 
   !            !allocate(expikr_set2%exps_storage(iq,idx)%exps(tc%nrt)) 
   !            allocate(expikr_set2%exps_storage(iq,idx)%exps(tc%nrt))
   !         enddo
   !      enddo
   !   enddo
   !enddo
   !tc => null()
!$omp parallel do schedule(guided) default(shared) private(iq, xq, wq, uq, expikr_tmp, ia, ja, ka, idx, tc, iexp)
      do iq = 1, numq1
         xq = num2kpt(qg%kpt(iq), qg%ndim)
         call solve_phonon_modes(ph, xq, wq, uq)
         !idx = 0
         !do ia = 1, pht%na 
         !   do ja = 1, pht%na 
         !      do ka = 1, pht%na
         !         idx = idx+1
         !         tc => pht%phit(idx)%trip_ph
         !         do iexp = 1, tc%nrt
         !            expikr_tmp = pure_exp_ikr(xq, pht%rvect_set(:,1,tc%rvect(iexp)))
         !            expikr_set1%exps_storage(iq,idx)%exps(iexp) = expikr_tmp
         !            expikr_tmp = pure_exp_ikr(xq, pht%rvect_set(:,2,tc%rvect(iexp)))
         !            expikr_set2%exps_storage(iq,idx)%exps(iexp) = expikr_tmp
         !         enddo
         !      enddo
         !   enddo
         !enddo

         if(use_mem)  uqgrid(:,:, iq) = uq(:,:)

      enddo
!$omp end parallel do


   ! if adaptive smearing, calculate phonon group velocities
   if ( adapt_smear_phph ) then
      allocate(velocity(3, qg%numb, qg%nk))
      velocity = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(iq)
      do iq = 1, qg%nk
         !call solve_phonon_velocity_findiff(qg, iq, velocity(:,:,iq))
         call solve_phonon_velocity_findiff(iq, qg%enk, qg%num_neighbors, &
            qg%neighbors(:,iq), qg%bweight, qg%bvec, velocity(:,:,iq))
      enddo 
!$omp end parallel do 
   endif
   
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
      
      if ( adapt_smear_phph ) then
         allocate( adsmear(tot_chnl) )
         adsmear = 0.0_dp
      endif

!$omp parallel default(shared) private(j, iq, ikq, xq, xkq, uk, uq, ukq, wq, &
!$omp& qidx, ltable, it, im, n, m, bands, vk, vkq, adsmear_tmp, ph3, expikr_set_tot, ia, ja, ka, idx, ix, ongrid)
!$omp do schedule(guided)
      do j = 1, nkq
         !allocate ( expikr_set_tot%exps_storage(1, pht%na*pht%na*pht%na) )
         iq = qq_pair(i)%iq(j)
         !ikq = qq_pair(i)%ikq(j)
         ikq= kpts_plus_invert( qg%kpt(ik), qg%kpt(iq), qg%ndim) + 1
         !
         xq = num2kpt(qg%kpt(iq), qg%ndim)
         xkq = num2kpt(qg%kpt(ikq), qg%ndim)

         ! compute expikr and store them
         !do idx = 1, pht%na * pht%na * pht%na 
         !   allocate(expikr_set_tot%exps_storage(1,idx)%exps(pht%phit(idx)%trip_ph%nrt))
         !   expikr_set_tot%exps_storage(1,idx)%exps = &
         !      conjg(expikr_set1%exps_storage(iq,idx)%exps(:) * expikr_set2%exps_storage(ikq, idx)%exps(:))
         !enddo

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
         !

         !! iq -> (k+q) - k
         !if( qidx .ne. qg%kpt(iq) ) then
         !   write(msg,'(2(3f12.6,2x))') num2kpt(qidx, qg%ndim), xq
         !   call errore('compute_scatter_ph3_psi2','failed: ' // trim(msg), 1)
         !endif
         
         !
         call match_table_all(qg%enk(:,ikq), qg%enk(:,ik), qg%enk(:,iq), wcut, e_thr, ltable)
         !sanity check
         if( qq_pair(i)%nchl(j) .ne. count(ltable) ) &
            call errore('compute_scatter_ph3_psi2','qq_pair(i)%nchl(j) .ne. count(ltable)',1)


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
                  vk = velocity(:,im,iq) ! try vq first
                  vkq = velocity(:,m,ikq)
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
               call solve_phi3mat(pht, xk, xq, xkq, qg%enk(:,ik), qg%enk(:,iq), qg%enk(:,ikq), &
                  uk, uq, ukq, n, im, m, ph3)
               
               !call solve_phi3mat(pht, expikr_set_tot, qg%enk(n,ik), qg%enk(im,iq), qg%enk(m,ikq), &
               !   uk(:,n), uq(:,im), ukq(:,m), ph3)
               ph3_psi2(it) = abs2(ph3)


               !ph3_psi2(it) = 0.05_dp/ryd2mev/ryd2mev
            endif
         enddo; enddo; enddo
         ! deallocate expikr
         !do idx = 1, pht%na * pht%na * pht%na 
         !   deallocate(expikr_set_tot%exps_storage(1,idx)%exps)
         !enddo
         !deallocate(expikr_set_tot%exps_storage)
      
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
   
   if( use_mem ) deallocate( uqgrid )
   
   deallocate( cryst_tau)

   call mp_barrier(inter_pool_comm)
   call stop_clock('comp_ph3psi2')
end subroutine compute_scatter_ph3_psi2
