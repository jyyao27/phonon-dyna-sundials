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

subroutine calc_band_structure()
   use pert_const, only: dp, ryd2ev
   use pert_data,  only: bg, epr_fid, kc_dim
   use pert_param, only: prefix, fklist
   use vector_list,only: vlist, load_vector_list
   use qe_mpi_mod, only: ionode, mp_split_pools, mp_sum, inter_pool_comm, npool
   use band_structure, only: electron_wann, init_electron_wann, solve_eigenvalue_vector
   implicit none
   type(vlist) :: kl
   integer :: ik, kst, kend
   character(len=80) :: fname
   real(dp), allocatable :: enk(:,:), xloc(:)
   !
   type(electron_wann) :: elec
   !
   call init_electron_wann(epr_fid, kc_dim, elec)
   !
   call load_vector_list(fklist, kl)
   call mp_split_pools(kl%nvec, kst, kend)
   if(npool > kl%nvec) &
      call errore('calc_bandstructure','too many pools (npool > #.k-points)',1)
    
   allocate(enk(elec%nb, kl%nvec))
   enk = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(ik)
   do ik = kst, kend
      call solve_eigenvalue_vector(elec, kl%vec(:,ik), enk(:,ik))
   enddo
!$omp end parallel do
   call mp_sum(enk, inter_pool_comm)
   enk = enk * ryd2ev

   allocate(xloc(kl%nvec))
   call generate_path(kl, bg, xloc)
   fname=trim(prefix)//'.bands'
   if(ionode) call output_bands(fname, elec%nb, kl%nvec, xloc, kl%vec, enk, 'eV')
end subroutine calc_band_structure

subroutine calc_phonon_spectra()
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
   real(dp), allocatable :: wq(:,:), xloc(:)
   !
   type(lattice_ifc) :: phon

   call init_lattice_ifc(epr_fid, qc_dim, phon)

   call load_vector_list(fqlist, ql)
   call mp_split_pools(ql%nvec, qst, qend)
   if(npool > ql%nvec) &
      call errore('calc_band_structure','too many pools (npool > #.q-points)',1)
    
   allocate(wq(3*phon%na, ql%nvec))
   wq = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(iq)
   do iq = qst, qend
      call solve_phonon_modes(phon, ql%vec(:,iq), wq(:,iq))
   enddo
!$omp end parallel do 
   call mp_sum(wq, inter_pool_comm)
   wq = wq * ryd2mev

   allocate(xloc(ql%nvec))
   call generate_path(ql, bg, xloc)
   fname=trim(prefix)//'.phdisp'
   if(ionode) call output_bands(fname, 3*phon%na, ql%nvec, xloc, ql%vec, wq, 'meV')
end subroutine calc_phonon_spectra

! compute phonon group velocity by finite difference
! using atomic units v = d omega /dq
! this is in cartesian coordinates!
subroutine calc_phonon_velocity()
   use pert_const, only: dp, ryd2mev
   use boltz_grid, only: grid, init_boltz_qgrid
   use boltz_grid_neighbors, only : set_grid_neighbors 
   use boltz_utils, only: kpt2num
   use pert_data,  only: bg, epr_fid, qc_dim
   use pert_param, only: fqlist, prefix
   use vector_list,only: vlist, load_vector_list
   use qe_mpi_mod, only: ionode, stdout, mp_split_pools, mp_sum, inter_pool_comm, npool
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc, solve_phonon_modes, &
      lattice_ifct, solve_phonon_velocity, solve_phonon_velocity_findiff
   implicit none
   type(vlist) :: ql
   type(grid) :: qg
   integer :: qst, qend, iq
   character(len=90) :: fname
   real(dp), allocatable :: xloc(:), vel(:,:,:)
   !
   type(lattice_ifc) :: phon

   call init_lattice_ifc(epr_fid, qc_dim, phon)
   if (ionode) then
      write(stdout,'(1x,a,/,5x,a)') "Computing phonon group velocities in lattice coordinates with finite difference"
   endif
   call load_vector_list(fqlist, ql)
   call mp_split_pools(ql%nvec, qst, qend)
   if(npool > ql%nvec) &
      call errore('calc_band_structure','too many pools (npool > #.q-points)',1)
    
   call init_boltz_qgrid(qg, phon, .true.)
   allocate(vel(3, qg%numb, ql%nvec ))
   vel = 0.0_dp
   call set_grid_neighbors(qg) 
!$omp parallel do schedule(guided) default(shared) private(iq)
   do iq = qst, qend
      !call solve_phonon_velocity_findiff(qg, kpt2num(ql%vec(:,iq), qg%ndim) + 1, vel(:,:,iq))
      call solve_phonon_velocity_findiff(iq, qg%enk, qg%num_neighbors, &
          qg%neighbors(:,iq), qg%bweight, qg%bvec, vel(:,:,iq))
   enddo
!$omp end parallel do 
   !call generate_path(ql, bg, xloc)
   fname=trim(prefix)//'.phvel'
   if(ionode) call output_velocity(fname, 3*phon%na, ql%nvec, ql%vec, vel)
   deallocate( vel)
end subroutine calc_phonon_velocity

subroutine generate_path(list, tmat, path)
   use pert_const, only: dp
   use vector_list, only: vlist
   implicit none
   type(vlist), intent(in) :: list
   real(dp), intent(in) :: tmat(3,3) !transformation matrix
   real(dp), intent(out) :: path(list%nvec)
   integer :: nvec, iv
   real(dp) :: dvec(3)
   real(dp), allocatable :: xvec(:,:)

   nvec = list%nvec
   allocate(xvec(3, nvec))

   xvec = list%vec
   call cryst_to_cart(nvec, xvec, tmat, 1) !kpt in cart. coor.
   path(1) = 0.0_dp
   do iv = 2, nvec
      dvec = xvec(:,iv) - xvec(:,iv-1)
      path(iv) = path(iv-1) + sqrt(dot_product(dvec, dvec))
   enddo
end subroutine generate_path

!output band structure or phonon dispersion
subroutine output_bands(fname, nb, nvec, xloc, vectors, bands, units)
   use pert_const, only: dp
   use pert_utils, only: find_free_unit
   use yaml_utils, only: ymlout
   use pert_param, only: calc_mode
   implicit none
   character(len=80), intent(in) :: fname
   integer, intent(in) :: nb, nvec
   real(dp), intent(in) :: xloc(nvec), vectors(3,nvec), bands(nb,nvec)
   character(len=*), intent(in) :: units
   integer :: iunit, iv, ib
   character(len=6), external :: int_to_char
   character(len=30) :: band_name, energy_name, index_name, point_prefix

   iunit = find_free_unit()
   open(iunit, file=trim(fname), form='formatted',status='unknown')
   do ib = 1, nb
      do iv = 1, nvec
         write(iunit,'(1x, f12.7, 2x, 3(1x,f10.5),2x,f16.10)') &
            xloc(iv), vectors(:,iv), bands(ib,iv)
      enddo
      write(iunit,'(a)') " "
   enddo

   close(iunit)

   if (calc_mode .eq. 'bands') then
      band_name = 'band'
      energy_name = 'band'
      index_name = 'band index'
      point_prefix = 'k'

   elseif (calc_mode .eq. 'phdisp') then
      band_name = 'mode'
      energy_name = 'phdisp'
      index_name = 'phonon mode'
      point_prefix = 'q'
   endif

   ! Output to the YAML file
   write(ymlout, '(/,a)') trim(calc_mode) // ":"
   write(ymlout, '(/,3x,a,i4)') 'number of ' // trim(band_name) // "s:", nb
   write(ymlout, '(/,3x,a)') trim(point_prefix) // '-path coordinate units: arbitrary'
   write(ymlout, '(3x,a)') trim(point_prefix) // '-path coordinates:'
   do iv = 1, nvec
      write(ymlout, '(6x,"-",1x,f12.7)') xloc(iv)
   enddo

   write(ymlout, '(/,3x,a)') trim(point_prefix) // '-point coordinate units: crystal'
   write(ymlout, '(3x,a)') trim(point_prefix) // '-point coordinates:'
   do iv = 1, nvec
      write(ymlout, '(6x,"-",1x, "[", 3(f10.5,",",2x), "]" )') vectors(:,iv)
   enddo

   write(ymlout, '(/,3x,a,1x,a)') trim(energy_name) // ' units:', trim(units)

   write(ymlout, '(3x,a)') trim(index_name) // ':'
   do ib = 1, nb
      write(ymlout, '(/,6x,a,a)') trim(int_to_char(ib)), ':'
      do iv = 1, nvec
         write(ymlout,'(9x,"-", f16.10)') bands(ib,iv)
      enddo
   enddo


end subroutine output_bands

!output group velocity
subroutine output_velocity(fname, nb, nvec, vectors, velocity)
   use pert_const, only: dp
   use pert_utils, only: find_free_unit
   implicit none
   character(len=80), intent(in) :: fname
   integer, intent(in) :: nb, nvec
   real(dp), intent(in) :: vectors(3,nvec), velocity(3,nb,nvec)
   !real(dp), intent(in) :: xloc(nvec), vectors(3,nvec), velocity(3,nb,nvec)
   integer :: iunit, iv, ib

   iunit = find_free_unit()
   open(iunit, file=trim(fname), form='formatted',status='unknown')
   do ib = 1, nb
      do iv = 1, nvec
         ! kelly 
         write(iunit,'(1x, 3(1x,f10.5),2x,3(1x,f10.5))') &
            vectors(:,iv), velocity(:,ib,iv)
      enddo
      write(iunit,'(a)') " "
   enddo
   close(iunit)
end subroutine output_velocity
