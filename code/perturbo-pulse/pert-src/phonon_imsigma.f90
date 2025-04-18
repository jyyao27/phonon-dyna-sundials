subroutine phonon_imsigma()
   use pert_const, only: dp
   use boltz_utils, only: kpt2num
   use boltz_scatter, only: boltz_ph3scat_setup, boltz_ph3scat_setup_fast
   use boltz_scatter_integral, only: ph3_scat_int, ph3_scat_int_interp2, ph3_scat_int_interp_half
   use boltz_dynamics_mod, only: init_disp_bose
   use boltz_grid, only: grid, init_boltz_qgrid, init_boltz_qgrid_vel
   use boltz_grid_neighbors, only: set_grid_neighbors
   use qe_mpi_mod, only: ionode, stdout, inter_pool_comm, mp_sum, my_pool_id, world_comm, mp_barrier
   use pert_output, only: output_imsigma, output_imsigma_yaml
   use pert_data,  only:  epr_fid, kc_dim, qc_dim
   use vector_list, only: vlist, load_vector_list, rand_vector_list
   use pert_param, only: fklist, fqlist, ntemper, temper, efermi, &
      cauchy_scale, sampling, nsamples, adapt_smear_phph, delta_smear_ph, divide_grid
   !
   use phonon_dispersion, only: lattice_ifc, lattice_ifct, &
                 init_lattice_ifc, init_lattice_ifct
   use phonon_anharm2, only: calc_selfph
   implicit none
   integer :: ntpr !ntpr defines phonon distribution
   integer :: ivec, it, iq, ib
   integer :: ik
   !load k- and q- list
   type(vlist) :: kl, ql
   type(grid) :: qg
   !real(dp), allocatable:: wq(:,:), imsigma(:,:,:,:)
   real(dp), allocatable:: wq(:,:), ph3col(:,:), imsigma(:,:, :,:)
   real(dp), allocatable :: qdist(:,:)
   real(dp), allocatable :: velocity(:,:,:) 
   integer, allocatable :: ql_ind(:)
   !
   type(lattice_ifc)   :: phon
   type(lattice_ifct)  :: phont
   logical :: boltz_scat_flag
   logical :: isfine
   integer :: ierr, myrank, nprocs
   
   isfine = divide_grid
   !isfine = .false.
   !check if temperatures have been readin 
   if(ntemper .eq. 0) call errore('phonon_imsigma',&
   'ftemper is not specified, missing temperatures and chemical potential',1)
   ntpr = size(temper)
   !load k- and q-grid

   !call load_vector_list(fklist, kl)
   if( trim(fqlist) .ne. '' ) then
      call load_vector_list(fqlist, ql)
   else
      call rand_vector_list(ql, nsamples, sampling, cauchy_scale)
   endif
   if(ionode) then
      !write(stdout,'(5x,a,i11)') 'total number of k-points:', kl%nvec
      ! kelly - this is not related to k points
      write(stdout,'(5x,a,i11)') 'total number of q-points:', ql%nvec
      write(stdout,'(5x,a,i11)') 'total number of q-points:', ql%nvec
   endif
   !init
   call init_lattice_ifc(epr_fid, qc_dim, phon)
   call start_clock('init_lattice_ifct')
   call init_lattice_ifct(epr_fid, qc_dim, phont, .true.)
   call stop_clock('init_lattice_ifct')
  
   ! reinitialize irreducible qgrids
   call init_boltz_qgrid(qg, phon, .true.)
   boltz_scat_flag = .true.
   ! find indices of qlist
   if (boltz_scat_flag) then
      allocate( ql_ind(ql%nvec) )
      write(*,*) ql%nvec
      ql_ind = 0
      iq = 0
      do ivec = 1, ql%nvec
         iq = kpt2num(ql%vec(:,ivec), qg%ndim) + 1
         !if (qg%irrk(ivec)+1 .ne. iq) call errore('phonon_imsigma', 'qg index error', 1)
         ql_ind(ivec) = iq
      enddo
         

      ! setup q-grid and e-ph scattering channel
      if ( adapt_smear_phph ) then
         call set_grid_neighbors(qg)
         call init_boltz_qgrid_vel(qg, phon, .True.) 
         write(stdout,'(5x,a, f10.5)') 'Using adaptive phonon-phonon smearing with maxsmear = ', delta_smear_ph 
      endif

      !call boltz_ph3scat_setup(qg, phon, phont, qg%ndim)
      
      if (isfine) then 
         !call mp_barrier(inter_pool_comm)
         !!! TEMP CHANGE by kyao
         !call boltz_ph3scat_setup(qg, phon, phont, isfine)
         call boltz_ph3scat_setup_fast(qg, phon, phont, isfine)
      else
         call boltz_ph3scat_setup(qg, phon, phont)
      endif

      allocate(imsigma(phon%nm, 1, ntpr, ql%nvec), wq(phon%nm, ql%nvec))
      allocate( qdist(qg%numb, qg%nk) )
      allocate( ph3col(qg%numb, qg%nk) )
      if (qg%numb .ne. phon%nm) then
         call errore('phonon_imsigma', 'qg%numb is not equal to phon%nm', 1)
      endif

      do it = 1, ntpr
         iq = 0
         qdist = 0.0_dp
         ph3col = 0.0_dp
         imsigma = 0.0_dp; wq = 0.0_dp   
         
         call init_disp_bose(qg, qdist, temper(it))
         if (isfine) then 
            !call ph3_scat_int(qg, qdist, ph3col, .true.) 
            call ph3_scat_int_interp_half(qg, qdist, ph3col, .true.) 
            !call ph3_scat_int_interp2(qg, qdist, ph3col, .true.) 
            !call mp_barrier(inter_pool_comm)
         else
            call ph3_scat_int(qg, qdist, ph3col, .true.) 
         endif
         do ivec = 1, ql%nvec
            iq = ql_ind(ivec)
            !write(*,*) ph3col(:,iq)
            imsigma(:, 1, it, ivec) = ph3col(:, iq)
            wq(:, ivec) = qg%enk(:, iq)
         enddo
      enddo
   !call mp_sum(imsigma, inter_pool_comm)
      call mp_barrier(inter_pool_comm)

   else
      ! kelly - this is not related to k points
      allocate(imsigma(phon%nm, phon%nm, ntpr, ql%nvec), wq(phon%nm, ql%nvec))
      imsigma = 0.0_dp; wq = 0.0_dp   
      write(stdout, '(5x, a)') 'Start calculating self-energy.'
      call calc_selfph(phon, phont, ql%vec, qg, temper, wq, imsigma)
   endif

   !output results
   if(ionode) then
      !call output_imsigma(.true. , temper, efermi, wq, imsigma, 'ph')
      call output_imsigma(.false., temper, efermi, wq, imsigma, 'ph')
      call output_imsigma_yaml(temper, efermi, wq, imsigma, ql%nvec, ql%vec)
   endif
   !release memory
   deallocate(wq, imsigma)
   if (boltz_scat_flag) deallocate(qdist, ph3col) 
   
   if(ionode) write(stdout,'(5x, a)') '>Im[Sigma] done.'
end subroutine phonon_imsigma
