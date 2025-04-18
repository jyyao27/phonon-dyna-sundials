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
   use boltz_utils, only: kpt2num, num2kpt
   use pert_data,  only: bg, epr_fid, qc_dim
   use pert_param, only: fqlist, prefix, boltz_qdim, ph_vel_numerical!, adapt_smear_phph
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

   real(dp), allocatable :: pheig(:), enk(:,:)
   complex(dp), allocatable :: phmode(:,:)
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
    
   if (ph_vel_numerical) then
      call init_boltz_qgrid(qg, phon, .true.)
      call set_grid_neighbors(qg) 
   else
      write(stdout,'(1x,a)') "Initialize qg manually"
      qg%numb = phon%nm
      qg%bmin = 1
      qg%bmax = phon%nm
      qg%ndim(1:3) = boltz_qdim(1:3)
   endif
   allocate(vel(3, qg%numb, ql%nvec ))
   vel = 0.0_dp
   allocate(pheig(phon%nm), phmode(phon%nm, phon%nm), enk(phon%nm, ql%nvec))
   pheig = 0.0_dp
   enk = 0.0_dp
   phmode = (0.0_dp, 0.0_dp)
!$omp parallel do schedule(guided) default(shared) private(iq,pheig,phmode)
   do iq = qst, qend
      !call solve_phonon_velocity_findiff(qg, kpt2num(ql%vec(:,iq), qg%ndim) + 1, vel(:,:,iq))
      if (ph_vel_numerical) then
         call solve_phonon_velocity_findiff(iq, qg%enk, qg%num_neighbors, &
            qg%neighbors(:,iq), qg%bweight, qg%bvec, vel(:,:,iq))

      else
         ! ql%vec is in lattice coordinates
         call solve_phonon_modes(phon, ql%vec(:,iq), pheig, phmode)
         call solve_phonon_velocity(phon, ql%vec(:,iq), pheig, phmode, vel(:,:,iq))
         enk(:,iq) = pheig
      endif
   enddo
!$omp end parallel do 

if (.not. ph_vel_numerical) then
   call mp_sum(enk, inter_pool_comm)
endif
   call mp_sum(vel, inter_pool_comm)

   !call generate_path(ql, bg, xloc)
   fname=trim(prefix)//'.phvel'
   !syp - units conversion here,
   !    - after we get vel from the following subroutine
   !    - if we want to get velocity (km/s) VS energy (THz)
   !      - vel = vel * alat * bohr2km * ryd2hz / (2*pi) at least for finite difference method
   !      - enk = enk * ryd2thz_h for THz, enk = enk * ryd2hz_hbar for rad/picosecond
   if (ph_vel_numerical) then
      !syp - note here, actually ql is not used for calculation,
      !      so for consistency, we need input full_grid for calculation
      !      then qg%ik equal to ql, this output should be correct then.
      if(ionode) call output_velocity2(fname, 3*phon%na, ql%nvec, ql%vec, vel, qg%enk)
   else
      if(ionode) call output_velocity2(fname, 3*phon%na, ql%nvec, ql%vec, vel, enk)
   endif
   deallocate( vel, pheig, phmode, enk)
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
! the phonon velocity should be times alat before using this subroutine
subroutine output_velocity(kg)
   use pert_param, only: prefix, tmp_dir
   use pert_const, only: dp, ryd2thz_hbar, ryd2hz_hbar, ryd2thz_h, bohr2ang
   use pert_utils, only: find_free_unit
   use boltz_grid, only: grid
   use pert_data, only: alat
   use boltz_utils, only: num2kpt
   implicit none
   type(grid), intent(inout) :: kg !kgrid
   !real(dp), intent(in) :: xloc(nvec), vectors(3,nvec), velocity(3,nb,nvec)
   integer :: iunit, ik, ib
   real(dp) :: xq(3),vel

   iunit = find_free_unit()
   open(iunit, file=trim(prefix) // ".phvel", form='formatted',status='unknown')
   write(iunit,'(1x, a, a, a, a, a)') '      q-vec (crystal)    '  , &
      '              velocity (Km/s)     ', '        |V|   ',&
      '  fre. (THz)  ', ' irrep.'    
   do ib = 1, kg%numb
      do ik = 1, kg%nk
            xq = num2kpt(kg%kpt(ik), kg%ndim)
            vel = sqrt(sum(kg%vnk(:,ib,ik)**2))*alat*bohr2ang*ryd2hz_hbar*1.0E-13_dp
            write(iunit,'(1x, 3(f8.5,1x),2x,3(f10.5,1x),f10.5,2x,f10.6,3x,I4)') &
               xq, kg%vnk(:,ib,ik)*alat*bohr2ang*ryd2hz_hbar*1.0E-13_dp, vel, &
               kg%enk(ib,ik)*ryd2thz_h, kg%kpt2ir(ik)
      enddo
      write(iunit,'(a)') " "
   enddo
   close(iunit)
end subroutine output_velocity
!output group velocity
subroutine output_velocity2(fname, nb, nvec, vectors, velocity, enk)
   use pert_const, only: dp
   use pert_utils, only: find_free_unit
   implicit none
   character(len=80), intent(in) :: fname
   integer, intent(in) :: nb, nvec
   real(dp), intent(in) :: vectors(3,nvec), velocity(3,nb,nvec)
   real(dp), intent(in) :: enk(nb,nvec)
   !real(dp), intent(in) :: xloc(nvec), vectors(3,nvec), velocity(3,nb,nvec)
   integer :: iunit, iv, ib

   iunit = find_free_unit()
   open(iunit, file=trim(fname), form='formatted',status='unknown')
   do ib = 1, nb
      do iv = 1, nvec
         !sypeng   
         write(iunit,'(1x, 3(1x,f10.5),2x,3(1x,f10.5),2x,f10.5)') &
            vectors(:,iv), velocity(:,ib,iv), enk(ib,iv)
         !kelly
         !  write(iunit,'(1x, 3(1x,f10.5),2x,3(1x,f10.5))') &
         !     vectors(:,iv), velocity(:,ib,iv)
      enddo
      write(iunit,'(a)') " "
   enddo
   close(iunit)
end subroutine output_velocity2
