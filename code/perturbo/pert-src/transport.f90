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
!> Main subroutine to compute transport properties like
!! electrical conductivity and mobility, thermal conductivity
!! and Seebeck coefficient. Key ingredient is to solve for 
!! the occupation changes (defined mfd in the code), and compute
!! the transport distribution function (tdf). Check the PERTURBO paper
!! Comput. Phys. Commun. 107970 (2021) for definitions of mfd and tdf
subroutine transport()
   use pert_const, only: dp, ryd2ev
   use pert_utils, only: converged
   use pert_output,only: load_imsigma
   use pert_data,  only: volume, epr_fid, qc_dim, kc_dim, alat
   use boltz_grid, only: grid, init_boltz_grid, boltz_grid_load
   use boltz_grid_neighbors, only: set_grid_neighbors
   use qe_mpi_mod, only: ionode, stdout, mp_barrier, inter_pool_comm
   use pert_param, only: ntemper, temper, efermi, ftemper, hole, doping, debug, &
      band_min, band_max, boltz_de, boltz_emax, boltz_emin, spinor, boltz_kdim, &
      trans_thr, prefix, boltz_nstep, boltz_qdim, full_ite, boltz_bfield, trans_info
   use boltz_trans_output, only: output_mobility, output_rates, output_tdf
   use boltz_trans_mod, only: trans_cond_calc, trans_density_calc
   use boltz_scatter_integral, only: rates_scat_int, trans_scat_int
   use boltz_scatter, only: boltz_scatter_setup
   use HDF5_utils
   !
   use bfield_utils, only: deriv,crossprod 
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   implicit none
   type(grid) :: kg
   logical :: read_rate
   integer(HID_T) :: file_id, group_id, group_id2
   integer :: nestep, i, it, ik, ib, io_nstep, niter(ntemper), ncomp
   real(dp) :: emin, emax, ef, tmpr, enk
   real(dp), allocatable :: ene(:), rates(:,:,:), cond(:,:,:), dos(:), &
      mfd0(:,:,:), mfd1(:,:,:), mfd_tmp(:,:,:), imsgm(:,:,:), tdf(:,:), &
      cmfd(:,:,:)
   !
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec

   if(ntemper .eq. 0) call errore('tranpsort', &
      'ftemper is not specified, missing temperatures and chemical potential',1)
   ! Initialize electron info into elec object
   call init_electron_wann(epr_fid, kc_dim, elec)
   ! Using elec object, define properties of kgrid object kg
   ! such as velocity, energy, etc. 
   call init_boltz_grid(kg, elec, band_min, band_max)
   ! setup up energy windows for transport calculations
   emin = max(boltz_emin, minval(kg%enk))
   emax = min(boltz_emax, maxval(kg%enk))
   !write(stdout,'(5x,a,2(1x,f12.6),/)') 'Energy window (eV):', emin*ryd2ev, emax*ryd2ev
   if(emin > emax) call errore('transport','illegal energy window',1)
   !
   call mp_barrier(inter_pool_comm)
   write(stdout,'(5x,a)') '>finish init grid.'
   
   !setup energy grid: emin - boltz_de : emax + boltz_de
   nestep = int( floor((emax-emin)/boltz_de) ) + 3
   allocate( ene(nestep), dos(nestep) )
   do i = 1, nestep
      ene(i) = emin + (i-2)*boltz_de
   enddo
   call trans_density_calc(kg, ene, temper, efermi, hole, spinor, volume, doping, dos=dos)
  

   !read imsigma files: imaginary part of the selfenergy in meV.
   allocate( rates(kg%numb, kg%nk, ntemper), imsgm(kg%numb, kg%nk_irr, ntemper) )
   call load_imsigma(imsgm, read_rate)

   !Only setup scattering if imsigma file is absent for rta solution and always
   !for ITA solution
   if(trans_info%is_ita .or. (.not. trans_info%is_ita .and. .not. read_rate)) then
      !Setup phonon properties inside phon object
      call init_lattice_ifc(epr_fid, qc_dim, phon)
      !Setup scattering info by computing e-ph on the fly
      !Or reading it from tmp files if it exists
      call boltz_scatter_setup(kg, elec, phon, boltz_qdim)
   endif
 
   !Compute scattering rates
   if( read_rate ) then
      !Read rates from prefix.imsigma file if it exists
      write(stdout,'(6x,a)') "- Read scattering rate from file: " // trim(prefix) // ".imsigma"
      !
      do it = 1, ntemper
      do ik = 1, kg%nk
      do ib = 1, kg%numb
         !Scattering rates in Rydberg atomic unit
         rates(ib,ik,it) = 2.0_dp*imsgm(ib, kg%kpt2ir(ik), it)
      enddo; enddo; enddo
   else
      if(.not. trans_info%is_ita) then
         write(stdout,'(6x,a)') "- imsigma file not detected for RTA - computing scattering rates on the fly"
      else
         write(stdout,'(6x,a)') "- Compute scattering rates on the fly."
      endif
      do it = 1, ntemper
         !Compute scattering rates on the fly
         call rates_scat_int(kg, temper(it), efermi(it), rates(:,:,it))
      enddo
      !output rates
      if(debug .and. ionode) then
         do it = 1, ntemper
            !If debug=.true., output the scattering rates on the k grid
            call output_rates(kg, temper(it), efermi(it), rates(:,:,it), it>1)
         enddo
      endif
   endif
   deallocate( imsgm )
   if(any(rates(:,:,:) < 1.0E-9_dp)) write(stdout,'(5x, a)') &
      "Warn (transport): scattering rates less than 10^-9 a.u."
   !
   write(stdout,'(5x,a)') '>finish init rates'

   !! now all data we needed are ready. start to do the actual work.
   ! allocate work space first.
   io_nstep = max(1, boltz_nstep)
   niter(:) = 1

   !ncomp = 3 : Both electric and magnetic fields calculations; No T field
   !ncomp = 6 : Required if full_ite=true(Only for 0 magnetic field)
   ncomp = merge(6, 3, full_ite)  !if full_ite is true, both mfd for E and T field are computed

   if(.not. trans_info%is_mag) then
      !For a non-magnetic field calculation, need 6 elements of conductivity tensor
      allocate( mfd0(ncomp, kg%numb, kg%nk), cond(ncomp*2, io_nstep, ntemper), tdf(ncomp*2, nestep))
      if(boltz_nstep > 0) allocate(mfd1(ncomp,kg%numb,kg%nk), mfd_tmp(ncomp,kg%numb,kg%nk))
   else
      !For a magnetic field calculation, need all 9 elements of the conductivity tensor
      allocate( mfd0(ncomp, kg%numb, kg%nk), cond(ncomp*3, io_nstep, ntemper), tdf(ncomp*3, nestep))
      allocate(mfd1(ncomp,kg%numb,kg%nk), mfd_tmp(ncomp,kg%numb,kg%nk),cmfd(ncomp,kg%numb,kg%nk))
      if(boltz_nstep < 1) call errore('transport','boltz_nstep zero for magnetic field calculation',1)
      call set_grid_neighbors(kg) !call kgrid neighbor data for derivative computation
   endif

   !open prefix_tdf.h5 files
   !write basic information- num energy steps, energies, dos, hole
   if(ionode) then
      call hdf_open_file(file_id, trim(prefix)//"_tdf.h5", status='NEW')
   endif
   
   write(stdout,'(5x,a)') '>start transport'
   write(stdout,'(5x,a)') '>progress: (current T / total T)'
   cond = 0.0_dp
   do it = 1, ntemper
      ef = efermi(it);  tmpr = temper(it)
      !compute Fnk = tau*vnk
      do ik = 1, kg%nk
      do ib = 1, kg%numb
         ! F for uniform E-field:  \tau_nk * v_nk
         mfd0(1:3,ib,ik) = kg%vnk(:,ib,ik) / rates(ib, ik, it)
         ! F for T gradienet:  \tau_nk * v_nk * (enk - \mu)
         if(full_ite) mfd0(4:6, ib,ik) = mfd0(1:3, ib,ik) * (kg%enk(ib,ik) - ef)
      enddo; enddo
      !compute conductivity and mobility under RTA, 
      ! or starting point for iter solution
      call trans_cond_calc(kg, ene, tmpr, ef, spinor, volume, mfd0, cond(:,1,it), tdf, trans_info%is_mag)
      !add the omitted factor, now tdf is fully in rydberg atomic unit
      tdf = tdf * alat * alat
      !Write down the information at this iteration into the hdf5 file
      if(ionode) call write_conf_hdf5(file_id, it, 1, nestep, kg%nk, kg%numb, tdf, mfd0, trans_info%is_mag, group_id,&
                 group_id2,trans_info%nomag_rta)

      !for iterative solution
      if(boltz_nstep > 0)  mfd1(:,:,:) = mfd0(:,:,:)
      !
      !compute correction from iterative solution
      do i = 2, io_nstep
         niter(it) = i
         !computed the integration terms (without the factor of tau0)
         mfd_tmp = 0.0_dp
         !If iterative calculation, compute the backscattering term
         if(trans_info%is_ita) call trans_scat_int(kg, temper(it), efermi(it), mfd1, mfd_tmp)

         !Compute the magnetic field term in the BTE
         !Sum_i ((vxB)_{i}\frac{d}{dk_{i}}) mfd
         if(trans_info%is_mag) then
            do ik=1, kg%nk
            do ib=1, kg%numb
               cmfd(:,ib,ik) = deriv(kg,crossprod(alat*kg%vnk(:,ib,ik), &
                                    boltz_bfield(:,it)),mfd1(:,:,:), ib, ik)&
                                    /rates(ib,ik,it)
            enddo;enddo
         endif

         !compute the new mean free displacement(mfd):
         do ik=1, kg%nk
         do ib=1, kg%numb
            if(.not. trans_info%is_mag) then
               !No B field term if a non-magnetic field calculation is performed
               mfd1(:,ib,ik) = mfd0(:,ib,ik) + mfd_tmp(:,ib,ik)/rates(ib,ik,it)
            else
               mfd1(:,ib,ik) = mfd0(:,ib,ik) + mfd_tmp(:,ib,ik)/rates(ib,ik,it)&
                                   + cmfd(:,ib,ik)
            endif
         enddo; enddo

         !compute conductivity from updated mfd1
         call trans_cond_calc(kg, ene, tmpr, ef, spinor, volume, mfd1, cond(:,i,it), tdf, trans_info%is_mag)
         tdf = tdf * alat * alat   
         !save tdf to hdf5 file
         if(ionode) call write_conf_hdf5(file_id, it, i, nestep, kg%nk, kg%numb, tdf, mfd1, &
                         trans_info%is_mag, group_id, group_id2)
         !exit if converged (wthin relative error within trans_thr)
         if(converged(cond(:,i-1,it), cond(:,i,it), trans_thr)) then
            !Write the last iteration outside the iterations folder as well
            if(ionode) call write_conf_hdf5(file_id, it, i, nestep, kg%nk, kg%numb, tdf, mfd1,&
                            trans_info%is_mag, group_id, group_id2, .true.)
            exit
         endif

         if(i == io_nstep) write(stdout,'(1x, a)') 'Warning: max number of iterations reached without convergence!'
      enddo
      !output the converged tdf in text format
      if(ionode) call output_tdf(ene, tdf, tmpr, ef, trans_info%is_mag, it>1)
      !output progress
      write(stdout, '(5x, i4, a, i4)') it, " /", ntemper
      flush(stdout)
   enddo
   write(stdout,'(5x, a)') '>finish transport.'
   
   if(ionode) then
      call write_basic_hdf5(file_id, hole, niter, nestep, ene, dos)
      call hdf_close_file(file_id)
      !add the omitted factor alat**2 to cond.
      cond = cond * alat * alat
      !output conductivity and mobility
      call output_mobility(temper, efermi, doping, cond, niter)
      !compute other transport coefficient 
      call trans_postproc()
   endif
   deallocate(ene, rates, mfd0, cond, tdf, dos)
   if(allocated(mfd1))  deallocate(mfd1)
   if(allocated(cmfd))  deallocate(cmfd)
   if(allocated(mfd_tmp)) deallocate(mfd_tmp)
   !
   return

   contains

      subroutine write_conf_hdf5(fid, tid, iter, nstep, nk, nb, tdf, mfd1, is_mag, gid, gid2, last_iter)
         use pert_const, only: dp
         use pert_param, only: temper, efermi, boltz_bfield
         use HDF5_utils
         implicit none
         logical , intent(in) :: is_mag               !< .true. if a magnetic field calculation
         logical , intent(in), optional :: last_iter  !< is .true. if last iteration 
         integer(HID_T), intent(in) :: fid            !< File id for HDF5 file
         integer(HID_T), intent(inout) :: gid         !< group id of the configurations group
         integer(HID_T), intent(inout) :: gid2        !< group id of the iterations group
         integer, intent(in) :: tid                   !< index of the configuration
         integer, intent(in) :: iter                  !< index of iteration
         integer, intent(in) :: nstep                 !< number of energy steps
         integer, intent(in) :: nk                    !< number of k points
         integer, intent(in) :: nb                    !< number of bands
         real(dp), intent(in) :: tdf(:,:)             !< Array containing TDF
         real(dp), intent(in) :: mfd1(:,:,:)          !< Array containing MFD
         !
         character(len=120) :: dset_name, dset_name_fnk, group_name
         character(len=6), external :: int_to_char
         logical :: is_last_iter

         !If last iteration
         is_last_iter = .false.
         if(present(last_iter)) is_last_iter = last_iter

         group_name = 'configuration_' // trim(int_to_char(tid)) 

         !Create a group only at the first iteration
         if(iter .eq. 1) then
            !Initialize the group and keep it open till
            !all iterations are printed per configuration
            call hdf_create_group(fid, trim(group_name))
            call hdf_open_group(fid, trim(group_name), gid)

            !Write basic data for each configuration
            call hdf_write_dataset(gid, 'temperature', temper(tid))
            call hdf_write_dataset(gid, 'efermi', efermi(tid))
            if(is_mag) call hdf_write_dataset(gid, 'magnetic_field',boltz_bfield(:,tid)) 

            !Create another group containing iterations
            !Only create the iterations group if calc_mode .eq. ITA
            !Change this if after separating calc_modes
            if(.not. (is_last_iter .and. (iter .eq. 1))) then
               call hdf_create_group(gid, 'iterations')
               call hdf_open_group(gid, 'iterations', gid2)
            endif
         endif

         !If last iteration
         if(is_last_iter) then
            if(hdf_exists(gid,'iterations'))call hdf_close_group(gid2) !Close the iterations group

            dset_name = 'tdf' 
            call hdf_write_dataset(gid, trim(dset_name), tdf(:,1:nstep))

            dset_name = 'mfd' 
            call hdf_write_dataset(gid, trim(dset_name), mfd1(:,1:nb,1:nk))

            call hdf_close_group(gid) !Close the configurations group

         else
            dset_name = 'tdf_' //  trim(int_to_char(iter))
            call hdf_write_dataset(gid2, trim(dset_name), tdf(:,1:nstep))

            dset_name = 'mfd_' //  trim(int_to_char(iter))
            call hdf_write_dataset(gid2, trim(dset_name), mfd1(:,1:nb,1:nk))
         endif

      end subroutine write_conf_hdf5

      subroutine write_basic_hdf5(fid, is_hole, nit, nstep, energy, dos_out)
         use pert_param, only: ntemper, temper, efermi, boltz_bfield
         use pert_const, only: dp
         use HDF5_utils
         implicit none
         integer(HID_T), intent(in) :: fid   !< file id for HDF5 file
         logical, intent(in) :: is_hole      !< .true. if hole=.true.
         integer, intent(in) :: nit(:)       !< Number of iterations for all configurations
         integer, intent(in) :: nstep        !< Number of energy steps         
         real(dp), intent(in) :: energy(:)   !< Energy grid points
         real(dp), intent(in) :: dos_out(:)  !< Density of states at each energy step
         !
         integer :: hole_i

         call hdf_write_dataset(fid,'num_configurations', ntemper)

         call hdf_write_dataset(fid, 'energy_grid', energy(1:nstep))
         call hdf_write_attribute(fid, 'energy_grid', 'ne', nstep)

         call hdf_write_dataset(fid, 'temperature_array', temper(1:ntemper))
         call hdf_write_attribute(fid, 'temperature_array', 'nt', ntemper)

         call hdf_write_dataset(fid, 'efermi_array', efermi(1:ntemper))
         call hdf_write_dataset(fid, 'num_iter_array', nit(1:ntemper))

         if(allocated(boltz_bfield)) call hdf_write_dataset(fid, 'bfield_array', boltz_bfield(:,1:ntemper))
         call hdf_write_dataset(fid, 'dos', dos_out(1:nstep))

         hole_i = merge(1, 0, is_hole)
         call hdf_write_attribute(fid, 'dos', 'hole', hole_i)
      end subroutine write_basic_hdf5
      !
end subroutine transport
