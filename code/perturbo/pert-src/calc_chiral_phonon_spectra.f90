subroutine calc_chiral_phonon_spectra
   use pert_const, only: dp, ryd2mev
   use pert_data,  only: bg, epr_fid, qc_dim
   use pert_param, only: fqlist, prefix
   use vector_list,only: vlist, load_vector_list
   use qe_mpi_mod, only: ionode, mp_split_pools, mp_sum, inter_pool_comm, npool
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc, solve_phonon_modes, &
                        lattice_ifct
   implicit none
   type(vlist) :: ql
   integer :: qst, qend, iq
   character(len=90) :: fname
   real(dp), allocatable :: wq(:,:), xloc(:)!, Sz_ph(:,:)
   integer :: nmode
   complex(dp), allocatable :: Sz(:,:), Sz_ph(:,:)
   complex(dp), allocatable:: phmode(:,:)
   !
   type(lattice_ifc) :: phon

   
   call init_lattice_ifc(epr_fid, qc_dim, phon)

   call load_vector_list(fqlist, ql)
   call mp_split_pools(ql%nvec, qst, qend)
   if(npool > ql%nvec) &
      call errore('calc_band_structure','too many pools (npool > #.q-points)',1)

   nmode = phon%nm
   allocate(Sz(nmode, nmode))
   call phonon_polarization_operator(nmode, Sz)  

   allocate(wq(3*phon%na, ql%nvec), phmode(3*phon%na, 3*phon%na))
   allocate(Sz_ph(3*phon%na, ql%nvec))
   wq = 0.0_dp; phmode = cmplx(0.0_dp, 0.0_dp)
!$omp parallel do schedule(guided) default(shared) private(iq, phmode)
   do iq = qst, qend
      call solve_phonon_modes(phon, ql%vec(:,iq), wq(:,iq), phmode)
      call phonon_polarization(nmode, phmode, Sz, Sz_ph(:,iq))
   enddo
!$omp end parallel do 
   call mp_sum(wq, inter_pool_comm)
   call mp_sum(Sz_ph, inter_pool_comm)
   wq = wq * ryd2mev

   allocate(xloc(ql%nvec))
   call generate_path(ql, bg, xloc)
   fname=trim(prefix)//'.phdisp'
   if(ionode) call output_disp_pol(fname, 3*phon%na, ql%nvec, xloc, ql%vec, wq, Sz_ph)
   if(ionode) then
      do iq = 1, ql%nvec
         write(*,*) iq, Sz_ph(:,iq)
      enddo
   endif
end subroutine calc_chiral_phonon_spectra


subroutine phonon_polarization_operator(nmode, Sz)
   use pert_const, only: dp
   implicit none
   integer, intent(in) :: nmode
   complex(dp), intent(out) :: Sz(nmode,nmode)
   complex(dp) :: Ri(1,nmode), Li(1,nmode)
   complex(dp) :: Sz_tmp(nmode, nmode)
   integer :: na, i, x, y
   
   Sz_tmp(:,:) = cmplx(0.0_dp, 0.0_dp)
   na = nmode/3
   do i = 1, na
      x = 3*(i-1)+1
      y = 3*(i-1)+2
      !initialization
      Ri(1,:) = cmplx(0.0_dp, 0.0_dp)
      !basis set for right polarization Ri and RiT
      Ri(1,x) = cmplx(1.0_dp, 0.0_dp)
      Ri(1,y) = cmplx(0.0_dp, 1.0_dp)
      write(*,*) 'i =', i, 'Ri', Ri, 'transpose', transpose(Ri)

      !initialization
      Li(1,:) = cmplx(0.0_dp, 0.0_dp)
      !basis set for right polarization Ri and RiT
      Li(1,x) = cmplx(1.0_dp, 0.0_dp)
      Li(1,y) = cmplx(0.0_dp, -1.0_dp)
     ! write(*,*) 'i =', i, 'Li', real(Li)
     ! write(*,*) 'i =', i, 'Li',aimag(Li)
      Sz_tmp = Sz_tmp + matmul(conjg(transpose(Ri)),Ri) -  matmul(conjg(transpose(Li)),Li)
   enddo
      Sz = Sz_tmp/2
     ! write(*,*) 'real Sz =' , real(Sz)
     ! write(*,*) 'imag Sz =' , aimag(Sz)
   
end subroutine phonon_polarization_operator

subroutine phonon_polarization(nm, phmode, Sz, Sz_ph)
   use pert_const, only: dp
   use pert_data,  only: mass
   implicit none
   integer, intent(in) :: nm
   complex(dp), intent(in) :: phmode(nm,nm), Sz(nm,nm)
   !real(dp), intent(out) :: Sz_ph(nm)
   complex(dp) :: sz_ph_tmp(1,1), o(1,1), Sz_ph(nm)
   complex(dp) :: phvect(nm,nm), phm(1,nm)
   integer :: i, j, ia

   phvect = cmplx(0.0_dp,0.0_dp)
   phm = cmplx(0.0_dp,0.0_dp)

   do i = 1, nm
      ia = (i - 1) / 3 + 1
      phvect(i, 1:nm) = phmode(i, 1:nm) * sqrt( mass(ia) )
   enddo

   !orthonormal test -passed for comlex mode
   !sanity check
    do i = 1, nm; do j = 1, nm
       phm(1,1) = phvect(1,i)
       phm(1,2) = phvect(2,i)
       phm(1,3) = phvect(3,i)
       phm(1,4) = phvect(4,i)
       phm(1,5) = phvect(5,i)
       phm(1,6) = phvect(6,i)
       !o = dot_product(phvect(i, :), phvect(j, :))
       !phm(1, :) = phvect(i, :)
       o = matmul(phm, conjg(transpose(phm)))
       write(*,*) 'orth test', i, j, 'o =', o
    enddo; enddo


   do i = 1, nm
      ia = (i - 1) / 3 + 1
      !phm(1, 1:nm) = phvect(i, 1:nm)
       phm(1,1) = phvect(1,i)
       phm(1,2) = phvect(2,i)
       phm(1,3) = phvect(3,i)
       phm(1,4) = phvect(4,i)
       phm(1,5) = phvect(5,i)
       phm(1,6) = phvect(6,i)
      write(*,*) 'band index', i, 'eigenvector real', real(phm(1,:))
      write(*,*) 'band index', i, 'eigenvector imag', aimag(phm(1,:))
      sz_ph_tmp = cmplx(0.0_dp, 0.0_dp)
      sz_ph_tmp = matmul(phm,matmul(Sz,conjg(transpose(phm))))
     ! write(*,*) 'band index', i, 'polarization', real(sz_ph_tmp)
     ! write(*,*) 'band index', i, 'polarization', aimag(sz_ph_tmp)
      Sz_ph(i) = sz_ph_tmp(1,1)
   enddo

end subroutine

!output band structure or phonon dispersion
subroutine output_disp_pol(fname, nb, nvec, xloc, vectors, bands, pol)
   use pert_const, only: dp
   use pert_utils, only: find_free_unit
   implicit none
   character(len=80), intent(in) :: fname
   integer, intent(in) :: nb, nvec
   real(dp), intent(in) :: xloc(nvec), vectors(3,nvec), bands(nb,nvec)
   complex(dp), intent(in) :: pol(nb, nvec)
   integer :: iunit, iv, ib

   iunit = find_free_unit()
   open(iunit, file=trim(fname), form='formatted',status='unknown')
   do ib = 1, nb
      do iv = 1, nvec
         write(iunit,'(1x, f12.7, 2x, 3(1x,f10.5),2x,f16.10, 2x, f10.5, 2x, f10.5)') &
            xloc(iv), vectors(:,iv), bands(ib,iv), real(pol(ib, iv)), aimag(pol(ib,iv))
      enddo
      write(iunit,'(a)') " "
   enddo
end subroutine output_disp_pol
