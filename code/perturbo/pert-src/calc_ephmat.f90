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

subroutine calc_ephmat()
   use pert_const, only: dp
   use pert_utils, only: abs2
   use vector_list,only: vlist, load_vector_list
   use pert_data,  only: mass, epr_fid, kc_dim, qc_dim
   use pert_param, only: fqlist, fklist, prefix, band_min, band_max, phfreq_cutoff
   use qe_mpi_mod, only: ionode, mp_split_pools, mp_sum, inter_pool_comm, npool, stdout
   use band_structure, only: electron_wann, init_electron_wann, solve_eigenvalue_vector
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc, solve_phonon_modes
   use elphon_coupling_matrix, only: elph_mat_wann, init_elph_mat_wann, eph_fourier_el_para, &
      eph_fourier_elph, eph_transform
   implicit none
   type(vlist) :: kl, ql
   character(len=80) :: fname
   integer :: ik, iq, kst, kend, im, ib, jb, numb, nmod, i, j
   real(dp):: tot_mass, xk(3), xq(3), xkq(3)
   real(dp), allocatable :: enk(:), ekq(:), wqt(:), g2(:,:,:), dpot(:,:,:), gmod(:,:,:), wq(:,:), dp2(:), gm2(:)
   complex(dp), allocatable :: uk(:,:), ukq(:,:), mq(:,:)
   complex(dp), allocatable :: g_kerp(:,:,:,:,:), gkq(:,:,:)

   type(lattice_ifc) :: ph
   type(electron_wann) :: el
   type(elph_mat_wann) :: ep

   call load_vector_list(fqlist, ql)
   call load_vector_list(fklist, kl)
   !sanity check
   if(ql%nvec * kl%nvec > 10000000) &
      call errore('calc_ephmat','too many k- and q-points (nk*nq > 10^7)',1)
   if(kl%nvec < npool) &
      call errore('calc_ephmat','too many pools (npool > #.k-points)', 1)
   !init
   call init_lattice_ifc(epr_fid, qc_dim, ph)
   call init_electron_wann(epr_fid, kc_dim, el)
   call init_elph_mat_wann(epr_fid, kc_dim, qc_dim, ep)

   numb = band_max - band_min + 1
   tot_mass = sum(mass)
   nmod = 3 * ph%na
   
   allocate( enk(el%nb), uk(el%nb, el%nb), wq(nmod, ql%nvec))
   allocate( g_kerp(3, ep%max_nrp, ep%nb, ep%nb, ep%na) )
   !
   allocate(dpot(nmod, ql%nvec, kl%nvec), gmod(nmod, ql%nvec, kl%nvec))
   dpot = 0.0_dp;  gmod = 0.0_dp;  wq = 0.0_dp
   
   call mp_split_pools(kl%nvec, kst, kend)

   do ik = kst, kend
      xk = kl%vec(:,ik)
      !electronic wavefunction at ik
      call solve_eigenvalue_vector(el, xk, enk, uk)
      !fourier transformation on Re. bottleneck part, openmp parallel inside.
      call eph_fourier_el_para(ep, xk, g_kerp)
!$omp parallel default(shared) private(iq, xq, xkq, ekq, &
!$omp& ukq, wqt, mq, gkq, g2, im, dp2, gm2, jb, ib, i, j)
      allocate( ekq(el%nb), ukq(el%nb, el%nb), wqt(nmod), mq(nmod, nmod) )
      allocate( gkq(ep%nb, ep%nb, nmod), g2(ep%nb, ep%nb, nmod), dp2(nmod), gm2(nmod) )
!$omp do schedule(guided) 
      do iq = 1, ql%nvec
         xq = ql%vec(:,iq)
         xkq = xk + xq
         !get electronic wavefunction at ikq
         call solve_eigenvalue_vector(el, xkq, ekq, ukq)
         !get phonon frequcies and eigen-displacement at iq
         call solve_phonon_modes(ph, xq, wqt, mq)
         !get e-ph matrix elements in wannier gauge and cart. coord.
         call eph_fourier_elph(ep, xq, g_kerp, gkq)
         !transfor to phonon mode and bloch gauge.
         call eph_transform(ep, xq, mq, uk, ukq, gkq)
         if(ik .eq. kst) wq(:,iq) = wqt(:)
         ! compute |g|^2
         g2 = abs2(gkq)

         dp2 = 0.0_dp
         gm2 = 0.0_dp
         do im = 1, nmod
            do jb = band_min, band_max
            do ib = band_min, band_max
               dp2(im) = dp2(im) + g2(ib,jb,im)
               if(wqt(im)>phfreq_cutoff) gm2(im) = gm2(im) + g2(ib,jb,im)*0.5_dp/wqt(im)
            enddo; enddo
         enddo
         !check degenerate phonon modes
         im = 1
         do while( im < nmod )
            i = 0
            do j = im+1, nmod
               if( abs(wqt(j) - wqt(im)) > 1.0E-12_dp ) exit
               ! there are degenerate modes
               i = i + 1
            enddo
            if(i > 0) then
               dp2( im:(im+i) ) = sum(dp2( im:(im+i) )) / real(i+1, dp)
               gm2( im:(im+i) ) = sum(gm2( im:(im+i) )) / real(i+1, dp)
            endif
            !update to the next mode
            im = im + i + 1
         enddo
         
         do im = 1, nmod
            gmod(im, iq, ik) = sqrt( gm2(im) / numb)
            dpot(im, iq, ik) = sqrt( dp2(im) * tot_mass / numb )
         enddo
      enddo
!$omp end do
   deallocate(ekq, ukq, wqt, mq, gkq, g2, dp2, gm2)
!$omp end parallel
   enddo
   call mp_sum(dpot, inter_pool_comm)
   call mp_sum(gmod, inter_pool_comm)

   fname = trim(prefix)//'.ephmat'
   if(ionode) call output_ephmat(fname, kl, ql, nmod, wq, dpot, gmod)
   deallocate(enk, uk, wq, g_kerp, dpot, gmod)
end subroutine calc_ephmat

subroutine output_ephmat(fname, kl, ql, nm, phdisp, defpot, gabs)
   use pert_const, only: dp, bohr2ang, ryd2ev, ryd2mev
   use pert_utils, only: find_free_unit
   use pert_data,  only: bg
   use vector_list, only: vlist
   use yaml_utils, only: ymlout
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: nm
   type(vlist), intent(in) :: kl, ql
   real(dp), intent(in) :: phdisp(nm, ql%nvec), defpot(nm, ql%nvec, kl%nvec), gabs(nm, ql%nvec, kl%nvec)
   !local variables
   integer :: uout, ik, iq, im, nk, nq
   real(dp), allocatable :: kloc(:), qloc(:)
   character(len=6), external :: int_to_char
      
   nk = kl%nvec
   nq = ql%nvec
   allocate( kloc(nk), qloc(nq) )

   call generate_path(kl, bg, kloc)
   call generate_path(ql, bg, qloc)

   uout = find_free_unit()
   open(unit=uout, file=trim(fname), status='unknown', form='formatted')
   write(uout,'(a,4x,a,6x,a,11x,a)') "#  ik      xk     iq      xq   imod", &
      "omega(meV)", "deform. pot.(eV/A)", "|g|(meV)"
   do ik = 1, nk
      do iq = 1, nq
      do im = 1, nm
         write(uout, '(1x, 2(i4,1x,f9.5,1x), 1x,i3.3,2x,f12.6, 2(2x,E22.12))') &
            ik, kloc(ik), iq, qloc(iq), im, phdisp(im,iq)*ryd2mev, &
            defpot(im,iq,ik)*ryd2ev/bohr2ang,  gabs(im,iq,ik)*ryd2mev 
      enddo; enddo
      write(uout, '(a)') '  '
   enddo
   close(uout)

   ! Output to the YAML file
   write(ymlout, '(/,a)') 'ephmat:'

   ! Units
   write(ymlout, '(/,3x,a)') 'phonon energy units: meV'

   write(ymlout, '(/,3x,a)') 'deformation potential units: eV/A'

   write(ymlout, '(/,3x,a)') 'e-ph matrix elements units: meV'

   write(ymlout, '(/,3x,a,i4)') 'number of phonon modes:', nm

   ! k- and q-points
   write(ymlout, '(/,3x,a)') 'k-path coordinate units: arbitrary'
   write(ymlout, '(3x,a)') 'k-path coordinates:'
   do ik = 1, nk
      write(ymlout, '(6x,"-",1x,f12.7)') kloc(ik)
   enddo

   write(ymlout, '(/,3x,a)') 'q-path coordinate units: arbitrary'
   write(ymlout, '(3x,a)') 'q-path coordinates:'
   do iq = 1, nq
      write(ymlout, '(6x,"-",1x,f12.7)') qloc(iq)
   enddo

   write(ymlout, '(/,3x,a)') 'k-point coordinate units: crystal'
   write(ymlout, '(3x,a)') 'k-point coordinates:'
   do ik = 1, nk
      write(ymlout, '(6x,"-",1x, "[", 3(f10.5,",",2x), "]" )') kl%vec(:,ik)
   enddo

   write(ymlout, '(/,3x,a)') 'q-point coordinate units: crystal'
   write(ymlout, '(3x,a)') 'q-point coordinates:'
   do iq = 1, nq
      write(ymlout, '(6x,"-",1x, "[", 3(f10.5,",",2x), "]" )') ql%vec(:,iq)
   enddo

   ! Values per mode
   write(ymlout, '(/,3x,a)') 'phonon mode:'

   do im = 1, nm
      write(ymlout, '(/,6x,a,a)') trim(int_to_char(im)), ':'

      write(ymlout, '(/,9x,a)') 'phonon energy:'
      do iq = 1, nq
         write(ymlout,'(12x,"-", f16.10)') phdisp(im,iq)*ryd2mev
      enddo

      write(ymlout, '(/,9x,a)') 'deformation potential:'
      do ik = 1, nk; do iq = 1, nq
         write(ymlout,'(12x,"-", 1x, E22.12)') defpot(im,iq,ik)*ryd2ev/bohr2ang
      enddo; enddo

      write(ymlout, '(/,9x,a)') 'e-ph matrix elements:'
      do ik = 1, nk; do iq = 1, nq
         write(ymlout,'(12x,"-", 1x, E22.12)') gabs(im,iq,ik)*ryd2mev
      enddo; enddo

   enddo

   deallocate(kloc, qloc)
end subroutine output_ephmat
