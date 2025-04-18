!transform the data structure of triplet to matrix, also drop the zero R's 
subroutine transform_fct()
   use pert_const, only: dp, ci, twopi
   use boltz_utils, only: kpt2num
   use hdf5_utils
   use qe_mpi_mod, only: ionode, stdout, inter_pool_comm, mp_sum, my_pool_id, world_comm, mp_barrier
   use pert_output, only: output_imsigma, output_imsigma_yaml
   use pert_data,  only: epr_fid, kc_dim, qc_dim
   use vector_list, only: vlist, load_vector_list, rand_vector_list
   use pert_param, only: prefix,fklist, fqlist, ntemper, temper, efermi, setup_rmap, &
      cauchy_scale, sampling, nsamples, adapt_smear_phph, delta_smear_ph
   use triplet_cell, only: trip_cell
   use phonon_dispersion, only: solve_phi3mat, solve_phonon_modes,lattice_ifc, lattice_ifct, &
                 init_lattice_ifc, init_lattice_ifct
   use third_order_fc
   implicit none
   type(fct_matrix) :: fct 
   integer :: ntpr !ntpr defines phonon distribution
   integer :: idx, irt
   integer :: a1,a2,a3,ir,ix1,ix2,ix3,im1,im2,im3,ir1,ir2,ir3
   complex(dp), allocatable :: ph3_trim(:,:,:) !phi3(\nu,\nu1,\nu2,q,q1) = phi(\nu q, \nu1 q_1,\nu2 -q-q_1)
   integer :: NRp
   real(dp), allocatable :: Rp(:,:) , rset(:,:)           !Rp grid
   integer, allocatable  :: rvect_index(:,:)
   !
   type(lattice_ifc)   :: phon
   type(lattice_ifct)  :: phont
   type(trip_cell), pointer :: tc
   real(dp) :: ifct(3,3,3)
   real(dp),allocatable :: r_phinorm(:), r_nonzero(:,:)
   integer,allocatable :: ir_nonzero(:), nonzero_index(:), im_nonzero(:), im_index(:) 
   integer :: icount, nat, nmod, nsize   
   
   integer(HID_T) :: group_id, fct_id 
      !local
   character(len=120) :: dset_name,fname

   !init
   call init_lattice_ifc(epr_fid, qc_dim, phon)
   call start_clock('init_lattice_ifct')

   call init_lattice_ifct(epr_fid, qc_dim, phont, .true.)
   call stop_clock('init_lattice_ifct')
   allocate(rvect_index(2,phont%nrvect))
   write(*,*) phont%nrvect, phont%max_nrt
   !output the total grid, it is a dual grid for both r1r2  
   open(101,file='rset.dat')
      do ir = 1,phont%nrvect 
         write(101,'(i10,6f10.5)') ir, phont%rvect_set(:,1,ir), phont%rvect_set(:,2,ir)
      enddo 
   close(101)
   !
   if(setup_rmap) then 
      write(*,'(A40)')'r is dumped to file. Run r-map-ir.py'
      return 
   endif  
   !reading the map from r to a integer: done by external python code   
   open(101,file='rset_index.dat',status='old')
      do ir = 1,phont%nrvect 
         read(101,*) rvect_index(:,ir)
      enddo 
   close(101)
   rvect_index = rvect_index + 1 !fortran starts from 1
   nmod = phon%nm; nat = phon%na
   !
   open(101,file='rset_irr.dat',status='old')
      read(101,*) nRp
      allocate(Rp(3,nRp))
      do ir = 1,nRp 
         read(101,*) Rp(:,ir)
      enddo 
   close(101)
   nsize = nRp*nmod
   write(*,'(A15,2i5)')'nRp, nsize = ', nRp, nsize !size of the irreducible R-grid 

   fct%nmod = nmod
   fct%nsize = nsize; fct%nr = nRp
   call allocate_fct(fct)

   !transform the phi3 from triplet to a matrix 
   fct%phi3_matrix = 0.d0 
   idx = 0
   do a3 = 1, nat
   do a2 = 1, nat
   do a1 = 1, nat
      idx = idx + 1
      tc => phont%phit(idx)%trip_ph
         do ir = 1, tc%nrt
            ifct = phont%phit(idx)%ifct(:,:,:,ir)      
            irt = tc%rvect(ir)   
            ir2 = rvect_index(1,irt)
            ir3 = rvect_index(2,irt)
            do ix1=1,3; do ix2 = 1,3; do ix3 = 1,3
               im1 = ix1 + (a1-1)*3
               im2 = ix2 + (a2-1)*3
               im3 = ix3 + (a3-1)*3
               fct%phi3_matrix( im3+(ir3-1)*phon%nm, im2+(ir2-1)*phon%nm, im1) = ifct(ix1,ix2,ix3)
            enddo;enddo;enddo 
         enddo
   enddo; enddo; enddo
   !
   allocate(r_phinorm(nsize))
   r_phinorm = 0.d0 
   do im1=1,nmod 
      do a1 = 1,nsize 
         r_phinorm(a1) = r_phinorm(a1) + maxval(abs(fct%phi3_matrix(:,a1,im1)))
      enddo  
   enddo 
   !drop all the (r,ka) that r_phinorm(r,ka) < 1e-14 
   icount = 0
   do a1 = 1,nsize 
      if(r_phinorm(a1)>1e-13) icount = icount + 1 
   enddo 
   allocate( nonzero_index(icount), ir_nonzero(icount), ph3_trim(icount,icount,nmod) )
   allocate( im_index(icount) )

   icount = 0; ph3_trim = 0.d0 
   do a1 = 1,nsize 
      if(r_phinorm(a1)>1e-13) then 
         icount = icount + 1 
         nonzero_index(icount) = a1
         im1 = mod(a1-1, nmod) + 1 
         im_index(icount) = im1 
         ir_nonzero(icount) = (a1-im1)/nmod + 1
      endif 
   enddo 
   write(*,'(A30,2i10)')'nsize, n-nonzero = ', nsize, icount
   
   do a1 = 1,icount   
      do a2 = 1,icount
         ph3_trim(a1,a2,:) = fct%phi3_matrix(nonzero_index(a1), nonzero_index(a2),:)
      enddo 
   enddo 

   !construct matrix 
   call deallocate_fct(fct)
   nsize = icount 
   fct%nsize = nsize 
   fct%nmod = nmod 
   fct%nr = nRp
   call allocate_fct(fct)

   fct%im_index = im_index 
   fct%rset = Rp
   fct%ir = ir_nonzero 
   fct%phi3_matrix = real(ph3_trim)
   !
   call body_order_decomp(fct)
   call write_fct(prefix, fct)

end subroutine

!phonon scattering rate, no mpi 
subroutine phonon_imsigma_fct()
   use omp_lib
   use pert_const, only: dp, pi, ryd2mev 
   use pert_utils, only : gauss, bose 
   use qe_mpi_mod, only: ionode, stdout, inter_pool_comm, mp_sum, my_pool_id, world_comm, mp_barrier
   use pert_data,  only:  epr_fid, kc_dim, qc_dim
   use vector_list, only: vlist, load_vector_list, rand_vector_list
   use pert_param, only: prefix, fklist, fqlist, ntemper, temper, phfreq_cutoff_ph, nsamples, delta_smear_ph
   use phonon_dispersion, only:  solve_phonon_modes,lattice_ifc, lattice_ifct, init_lattice_ifc, init_lattice_ifct
   use third_order_fc
   implicit none 
   type(lattice_ifc) :: ph 
   type(fct_matrix) :: fct 
   type(vlist) :: ql, ql_1
   real(dp) :: qpt(3), q1(3), q2(3)
   complex(dp), allocatable :: phi3_qq(:,:,:), u1(:,:), u2(:,:), uq(:,:)
   complex(dp), allocatable :: phiqr(:,:,:,:)
   real(dp), allocatable :: w1(:), w2(:), wqt(:), wq(:,:),bose_q(:), bose_q1(:),bose_q2(:)
   real(dp), allocatable :: phi3_2(:,:,:), imsgm(:,:,:)
   real(dp) :: ph_n, ph_n1, ph_n2, dt1, dt2, dt3, qwt, scat, x_imse 
   integer :: nmod   
   integer :: iq, iq_1, im1, im2, im, it, isvd, nth_omp, id_omp  
   logical :: active_channel 
   character(len=20) :: fname 
   
   if(ionode)write(*,*) 'init ph and load fct'
   call init_lattice_ifc(epr_fid, qc_dim, ph)
   call start_clock('init_lattice_ifct')
   call load_fct(prefix, fct)
   nmod = ph%nm

   !load in klist for gamma_q  
   call load_vector_list(fklist, ql)
   !load qlist for integeration q1 
   call load_vector_list(fqlist, ql_1)
   if(ionode) then 
      write(*,'(A10,i5)')'nq = ',ql%nvec
      write(*,'(A20,i5,f10.5)')'nq1,sum(weight) = ',ql_1%nvec, sum(ql_1%weight)
   endif 
   !
   allocate(wq(ph%nm, ql%nvec))
   allocate(imsgm(ph%nm, ntemper, ql%nvec))
   imsgm = 0.d0 

   
   write(*,'(A40)')'starting main loop for imse(ph) : '
   !$omp parallel default(shared) private(iq, qpt,wqt,uq,iq_1,qwt,q1,q2, &
   !$omp& active_channel,w1,u1,w2,u2,im,im1,im2,isvd,phi3_qq,phi3_2,dt1,dt2,dt3,&
   !$omp& ph_n1, ph_n2, scat, it, x_imse,id_omp,phiqr)  
   allocate(phi3_2(ph%nm,ph%nm,ph%nm))
   allocate(phi3_qq(ph%nm,ph%nm,ph%nm), phiqr(fct%nr,ph%nm,ph%nm,ph%nm))
   allocate(wqt(ph%nm),w1(ph%nm),w2(ph%nm))
   allocate(u1(ph%nm,ph%nm),u2(ph%nm,ph%nm),uq(ph%nm,ph%nm))
   id_omp = omp_get_thread_num()
   nth_omp = omp_get_num_threads()
   do iq = 1,ql%nvec
      if(mod(iq,nth_omp).ne.id_omp) cycle 
      qpt = ql%vec(:,iq)
      call solve_phonon_modes(ph, qpt, wqt, uq)
      wq(:,iq) = wqt 
      !Bose Einstein distribution at q 
      call phirr_to_phiqr(fct, qpt, phiqr)
      !write(*,'(E20.10)') maxval(abs(phiqr))

      do iq_1 = 1,ql_1%nvec     
         q1 = ql_1%vec(:,iq_1)
         qwt = ql_1%weight(iq_1)
         !momentum conservation 
         q2 = -qpt - q1 
         
         call solve_phonon_modes(ph, q1, w1, u1)
         call solve_phonon_modes(ph, q2, w2, u2)

         !check energy conservation 
         active_channel = .false. 
         do im = 1,ph%nm 
            if(active_channel) cycle 
            if(wqt(im)<phfreq_cutoff_ph) cycle 
            do im1 = 1,ph%nm
               if(w1(im1)<phfreq_cutoff_ph) cycle 
               do im2 = 1,ph%nm 
                  if(w2(im2)<phfreq_cutoff_ph) cycle 
                  if(abs(wqt(im) - w1(im1) - w2(im2))/delta_smear_ph < 3.d0) active_channel = .true.
                  if(abs(wqt(im) + w1(im1) - w2(im2))/delta_smear_ph < 3.d0) active_channel = .true. 
                  if(abs(wqt(im) - w1(im1) + w2(im2))/delta_smear_ph < 3.d0) active_channel = .true. 
               enddo 
            enddo 
         enddo 
         if(.not. active_channel) cycle 
         
         call phiqr_to_phiqq(fct, q1, phiqr, phi3_qq)
         !write(*,'(E20.10)') maxval(abs(phi3_qq))
         !transfer \phi to eigenmode 
         !\phi(\nu_2 q_2, \nu_1 q_1, \nu q)
         call phi3_transform(nmod, u2, u1, uq, phi3_qq)

         !V3 without frequency coeff 
         phi3_2 = abs(phi3_qq)**2

         !summation over \nu \nu_1,\nu_2
         do it = 1,ntemper   
            do im = 1,ph%nm 
               x_imse = 0.d0 
               if(wqt(im)<1e-6_dp) cycle 
               do im1 = 1,ph%nm
                  if(w1(im1)<phfreq_cutoff_ph) cycle 
                  do im2 = 1,ph%nm 
                     if(w2(im2)<phfreq_cutoff_ph) cycle 
                     !gaussians 
                     dt1 = gauss( delta_smear_ph, wqt(im) - w1(im1) - w2(im2)) 
                     dt2 = gauss( delta_smear_ph, wqt(im) + w1(im1) - w2(im2)) 
                     dt3 = gauss( delta_smear_ph, wqt(im) - w1(im1) + w2(im2)) 
                     !loop over temperature  
                  
                     ph_n1 = bose(temper(it), w1(im1))
                     ph_n2 = bose(temper(it), w2(im2))
                     scat = (ph_n1 + ph_n2 + 1.d0) * dt1 + (ph_n1 - ph_n2) * (dt2 - dt3) 
                     scat = scat * phi3_2(im2, im1, im)/(w2(im2)*w1(im1)*wqt(im)) * qwt * pi/ 16.d0 
                     x_imse = x_imse + scat 
                  enddo 
               enddo 
               !$omp atomic update
               imsgm(im,it,iq) = imsgm(im,it,iq) + x_imse 
            enddo
         enddo  
      enddo 
   enddo 
   deallocate(phi3_2)
   deallocate(phi3_qq,phiqr)
   deallocate(wqt,w1,w2)
   deallocate(u1,u2,uq)
   !$omp end parallel

   !output to file 
   !output to file 
   open(101, file='imsigma_ph.dat') 
   do it = 1, ntemper 
      write(101,'(A15, f10.5)')'# temp (meV) = ', temper(it)*ryd2mev
      do iq = 1, ql%nvec 
         do im = 1, ph%nm 
            write(101,'(2E20.10)')wq(im,iq)*ryd2mev, imsgm(im, it, iq)*ryd2mev
         enddo  
      enddo 
   enddo 
   close(101)
end subroutine

