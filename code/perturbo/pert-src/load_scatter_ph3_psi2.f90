subroutine load_scatter_ph3_psi2(nq_loc, qq_pair, num_qq_pair, ph3_scat, tsmear, q_start_ind)
   use hdf5_utils
   use pert_param, only: prefix, tmp_dir, adapt_smear_phph
   implicit none
   integer, intent(in) :: nq_loc, num_qq_pair
   type(channel_ph), intent(in) :: qq_pair(nq_loc)
   type(t_scatt_ph), intent(out) :: ph3_scat(num_qq_pair)
   type(t_scatt_smear), allocatable, intent(out), optional :: tsmear(:) 
   integer, intent(out), optional :: q_start_ind(nq_loc)
   !
   integer :: n, i, ik, nkq, j, tot_chnl, nchl, it, ib
   integer, allocatable :: ist(:), bnd_idx(:)
   real(dp), allocatable :: ph3_psi2(:)
   real(dp), allocatable :: adsmear(:)
   !
   logical :: has_file
   integer(HID_T) :: file_id
   character(len=120) :: fname, dset_name
   character(len=6), external :: int_to_char
   integer :: i_scat, i_scat_chl
   
   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_ph3_psi2_p" 
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   ! check if file exist
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('load_scatter_ph3_psi2', "Missing file: "//trim(fname), 1)
   ! open file
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')

   it = 0
   i_scat = 0
   i_scat_chl = 0
   n = 0
   if (present( q_start_ind)) q_start_ind = 0 
   do i = 1, nq_loc
      if (present(q_start_ind)) then
         q_start_ind(i) = n + 1
      endif
      ik = qq_pair(i)%ik
      nkq = qq_pair(i)%nkq
      if (nkq < 1) cycle
      
      allocate( ist( nkq ) )
      ist(1) = 0 
      do j = 1, nkq-1
         ist(j+1) = ist(j) + qq_pair(i)%nchl(j)
      enddo

      tot_chnl = sum( qq_pair(i)%nchl(:) )
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
      do j = 1, nkq
         i_scat = i_scat + 1
         n = n + 1
         !
         ph3_scat(n)%ik = ik
         ph3_scat(n)%iq = qq_pair(i)%iq(j)
         !ph3_scat(n)%ikq = qq_pair(i)%ikq(j)
         !
         nchl = qq_pair(i)%nchl(j)
         ph3_scat(n)%nchl = nchl
         allocate( ph3_scat(n)%bands_index(nchl), ph3_scat(n)%eph_g2(nchl) )
         !allocate( ph3_scat(n)%iqf(2,64))
         !ph3_scat(n)%iqf(:,:) = 0

         if (adapt_smear_phph ) allocate( tsmear(n)%adsmear(nchl) )
         it = ist(j)
         do ib = 1, nchl
            it = it + 1
            i_scat_chl = i_scat_chl + 1
            !write(*,*) ik, it, i_scat, i_scat_chl
            ph3_scat(n)%bands_index(ib) = bnd_idx(it)
            ph3_scat(n)%eph_g2(ib) = ph3_psi2(it)
            if (adapt_smear_phph ) then
               tsmear(n)%adsmear(ib) = adsmear(it)
               !ph3_scat(n)%adsmear(ib) = adsmear(it)
            endif
         enddo
      enddo
      deallocate(ist, bnd_idx, ph3_psi2)
      if ( adapt_smear_phph ) deallocate( adsmear )
      !
   enddo
   !
   call hdf_close_file(file_id)
   call mp_barrier( inter_pool_comm )
end subroutine load_scatter_ph3_psi2
