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

Module mod_transport
   use pert_const, only: dp, ryd2ev, ryd2mev, bohr2ang, timeunit
   use pert_utils, only: converged, converged_vec
   use pert_output,only: load_imsigma
   use pert_data,  only: volume, epr_fid, qc_dim, kc_dim, alat
   use boltz_grid, only: grid, init_boltz_grid, boltz_grid_load, init_boltz_qgrid, &
      init_boltz_qgrid_vel_sym
   use boltz_grid_neighbors, only: set_grid_neighbors
   use qe_mpi_mod, only: ionode, stdout, mp_barrier, inter_pool_comm
   use pert_param, only: ntemper, temper, efermi, ftemper, hole, doping, debug, &
      band_min, band_max, boltz_de, boltz_de_ph, boltz_emax, boltz_emin, spinor, boltz_kdim, &
      trans_thr, prefix, boltz_nstep, boltz_qdim, full_ite, boltz_bfield, &
      trans_info, adapt_smear_phph, ph_vel_numerical, ph_e, drag, symm, converge_diag
   use boltz_trans_output, only: output_mobility, output_rates, output_tdf, output_ph_cond, convert_cond
   use boltz_trans_mod, only: trans_cond_calc_ph, trans_cond_calc, trans_density_calc, trans_impose_kelveinonsager
   use boltz_scatter_integral, only: rates_scat_int, trans_scat_int, rates_phe_scat_int, &
      drag_scat_4e, rates_ph3_scat_int, trans_ph3_scat_int, drag_scat_4ph 
   use boltz_scatter, only: boltz_scatter_setup, boltz_ph3scat_setup
   use HDF5_utils
   !
   use bfield_utils, only: deriv,crossprod 
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, lattice_ifct, init_lattice_ifc, init_lattice_ifct
   implicit none
   private

   public :: transport, transport_ph

   contains
      subroutine transport(kg_out, qg, it_out, e_niter, phon_out, elec_out, mfd_N, ph_tdf, ph_cond, mfd_F, lastornot)
         !- syp
         implicit none
         type(grid), intent(in), optional :: kg_out
         type(grid), intent(in), optional :: qg
         integer, intent(in), optional :: it_out
         integer, intent(inout), optional :: e_niter(:)
         type(lattice_ifc), intent(in), optional :: phon_out
         type(electron_wann), intent(in), optional :: elec_out
         real(dp), intent(in), optional :: mfd_N(:, :, :)!< mfd for e-BTE, (ncomp, qg%numb, qg%nk)
         real(dp), intent(in), optional :: ph_tdf(:,:)!< tdf for ph-BTE, (ncomp*2, nestep)
         real(dp), intent(in), optional :: ph_cond(:)!< conductivity for ph-BTE, (ncomp*2)
         real(dp), intent(out), optional :: mfd_F(:, :, :)!< mfd for e-BTE, (ncomp, kg%numb, kg%nk)
         logical, intent(in), optional :: lastornot
      
         type(grid) :: kg
         logical :: read_rate, lastornot_local
         integer(HID_T) :: file_id, group_id, group_id2
         integer :: nestep, i, it, ik, ib, io_nstep, niter(ntemper), ncomp
         integer :: it_start, it_end
         real(dp) :: emin, emax, ef, tmpr, enk
         real(dp), allocatable :: ene(:), rates(:,:,:), cond(:,:,:), dos(:), &
            mfd0(:,:,:), mfd1(:,:,:), mfd_tmp(:,:,:), imsgm(:,:,:), tdf(:,:), &
            cmfd(:,:,:), drag_4e(:,:,:)
         integer, allocatable :: converg_vec(:)
         real(dp) :: cond_unit
         !
         type(lattice_ifc)   :: phon
         type(electron_wann) :: elec
      
      
      
         lastornot_local = .false.
      
         if(ntemper .eq. 0) call errore('tranpsort', &
            'ftemper is not specified, missing temperatures and chemical potential',1)
      
         if(drag)then 
            if(.not. present(mfd_N) .or. .not. present(kg_out) .or. .not. present(phon_out) .or. .not. present(elec_out) &
               .or. .not. present(mfd_F) .or. .not. present(it_out) .or. .not. present(qg)) then
               call errore('transport_drag_e','mfd_N, mfd_F, qg, it_out, kg_out, phon_out, elec_out must be present for drag calculation',1)
            endif
            kg = kg_out
            elec = elec_out
            niter = e_niter
            it_end = it_out
            it_start = it_out
            lastornot_local = lastornot
            if(ionode .and. debug) write(stdout,'(10x,a)') '> kg grid and elec object are passed from phonon BTE.'

            !initialize
            mfd_F = 0.0_dp
            call mp_barrier(inter_pool_comm)
         else
            it_start = 1
            it_end = ntemper
            ! Initialize electron info into elec object
            call init_electron_wann(epr_fid, kc_dim, elec)
            ! Using elec object, define properties of kgrid object kg
            ! such as velocity, energy, etc. 
            call init_boltz_grid(kg, elec, band_min, band_max)
            lastornot_local = .true.
            write(stdout,'(10x,a)') '>finish init grid.'
         endif


         ! setup up energy windows for transport calculations
         emin = max(boltz_emin, minval(kg%enk))
         emax = min(boltz_emax, maxval(kg%enk))
        !if(.not. drag) write(stdout,'(5x,a,2(1x,f12.6),/)') 'Energy window (eV):', emin*ryd2ev, emax*ryd2ev
         if(emin > emax) call errore('transport','illegal energy window',1)
         !
         call mp_barrier(inter_pool_comm)
         
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
      
         if(drag)then 
            !for drag, we must need phon either ita or rta
            phon = phon_out
         else
            !Only setup scattering if imsigma file is absent for rta solution and always
            !for ITA solution
            if(trans_info%is_ita .or. (.not. trans_info%is_ita .and. .not. read_rate)) then
               !Setup phonon properties inside phon object
               call init_lattice_ifc(epr_fid, qc_dim, phon)
               !Setup scattering info by computing e-ph on the fly
               !Or reading it from tmp files if it exists
               call boltz_scatter_setup(kg, elec, phon, boltz_qdim)
            endif
         endif
      
         !Compute scattering rates
         if( read_rate ) then
            !Read rates from prefix.imsigma file if it exists
            write(stdout,'(12x,a)') "- Read scattering rate from file: " // trim(prefix) // ".imsigma"
            !
            do it = 1, ntemper
            do ik = 1, kg%nk
            do ib = 1, kg%numb
               !Scattering rates in Rydberg atomic unit
               rates(ib,ik,it) = 2.0_dp*imsgm(ib, kg%kpt2ir(ik), it)
            enddo; enddo; enddo
         else
            if(.not. trans_info%is_ita) then
               if(.not. drag) write(stdout,'(12x,a)') "- imsigma file not detected for RTA - computing scattering rates on the fly"
            else
               if(.not. drag) write(stdout,'(12x,a)') "- Compute scattering rates on the fly."
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
         if(.not. drag .and. any(rates(:,:,:) < 1.0E-9_dp)) write(stdout,'(5x, a)') &
            "Warn (transport): scattering rates less than 10^-9 a.u."
         !
         if(.not. drag) write(stdout,'(10x,a)') '>finish init rates'
      
         !! now all data we needed are ready. start to do the actual work.
         ! allocate work space first.
         io_nstep = max(1, boltz_nstep)
         niter(:) = 1
      
         !ncomp = 3 : Both electric and magnetic fields calculations; No T field
         !ncomp = 6 : Required if full_ite=true(Only for 0 magnetic field)
         if(drag) then
            ncomp = size(mfd_N,1)
         else
            ncomp = merge(6, 3, full_ite)  !if full_ite is true, both mfd for E and T field are computed
         endif
      
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

         if (converge_diag) then
            allocate(converg_vec(ncomp))
            if( ncomp == 3 ) then
               converg_vec = (/1, 3, 6/)
            elseif( ncomp == 6 ) then
               converg_vec = (/1, 3, 6, 7, 9, 12/)
            else
               call errore('transport','ncomp must be 3 or 6',1)
            endif
         else
            !use forall to traverse from 1 to ncomp
            if(trans_info%is_mag) then
               allocate(converg_vec(3*ncomp))
               forall(i=1:3*ncomp) converg_vec(i) = i
            else
               allocate(converg_vec(2*ncomp))
               forall(i=1:2*ncomp) converg_vec(i) = i
            endif
         endif
      
         !open prefix_tdf.h5 files
         !write basic information- num energy steps, energies, dos, hole
         if(ionode .and. lastornot_local) then
            if (drag .and. it_start > 1) then
               call hdf_open_file(file_id, trim(prefix)//"_tdf.h5", status='OLD', action='WRITE')!, action='READWRITE')
            else
               call hdf_open_file(file_id, trim(prefix)//"_tdf.h5", status='NEW')
            endif
         endif
      
        !write(stdout,'(5x,a)') '>start transport'
        !if(.not. drag) write(stdout,'(5x,a)') '>progress: (current T / total T)'
         cond = 0.0_dp
         do it = it_start, it_end
            ef = efermi(it); tmpr = temper(it) 
            !we do the drag term for electrons using the mfd_N from phonon parts
            if( drag ) then
               allocate( drag_4e(ncomp, kg%numb, kg%nk) )
               call drag_scat_4e(kg, qg, tmpr, ef, mfd_N, drag_4e)
            endif
      
            !compute Fnk = tau*vnk
            mfd0 = 0.0_dp
            do ik = 1, kg%nk
            do ib = 1, kg%numb
               ! F for uniform E-field:  \tau_nk * v_nk
               if (abs(rates(ib,ik,it)) > 1E-13_dp) then
                  mfd0(1:3,ib,ik) = kg%vnk(:,ib,ik) / rates(ib, ik, it)
                  if(drag) drag_4e(1:ncomp,ib,ik) = drag_4e(1:ncomp,ib,ik) / rates(ib, ik, it)
               endif
               ! F for T gradienet:  \tau_nk * v_nk * (enk - \mu)
               if(full_ite) mfd0(4:6, ib,ik) = mfd0(1:3, ib,ik) * (kg%enk(ib,ik) - ef)
            enddo; enddo
            !compute conductivity and mobility under RTA,
            ! or starting point for iter solution
            if (drag) then
               mfd1 = mfd0 + drag_4e
            else
               mfd1 = mfd0
            endif 

            if(drag) then
               if(present(ph_cond)) then
                  if(ionode .and. debug) write(stdout,'(12x,a)') 'Imposing Kelvin-Onsager relation'
                  call trans_impose_kelveinonsager(kg, ene, tmpr, ef, spinor, volume, -ph_cond(7:12), mfd1)
               endif
            endif
            call trans_cond_calc(kg, ene, tmpr, ef, spinor, volume, mfd1, cond(:,1,it), tdf, trans_info%is_mag)
            ! note: ph_cond(~T, ~E) = (thermal conductivity, ph_alpha / T)
            !add the omitted factor, now tdf is fully in rydberg atomic unit
            tdf = tdf * alat * alat
            call convert_cond(cond_unit)
            if(ionode .and. drag .and. debug) write(stdout,'(12x,a6,I3,a4,12(E13.6,1x),a8)') 'e-BTE:',1,'-RTA:',cond(:,1,it)*alat*alat*cond_unit,' 1/Ohm/m'
            !Write down the information at this iteration into the hdf5 file
            if(ionode .and. lastornot_local) call write_conf_hdf5(file_id, it, 1, nestep, kg%nk, kg%numb, tdf, mfd1, trans_info%is_mag, group_id,&
                    group_id2,trans_info%nomag_rta)
      
            !for iterative solution
            if(boltz_nstep > 0 .or. drag)  mfd1(:,:,:) = mfd0(:,:,:)
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
                     if(abs(rates(ib,ik,it)) > 1E-13_dp) then
                        mfd1(:,ib,ik) = mfd0(:,ib,ik) + mfd_tmp(:,ib,ik)/rates(ib,ik,it)
                     else
                        mfd1(:,ib,ik) = mfd0(:,ib,ik)
                     endif
                  else
                     if(abs(rates(ib,ik,it)) > 1E-13_dp) then
                        mfd1(:,ib,ik) = mfd0(:,ib,ik) + mfd_tmp(:,ib,ik)/rates(ib,ik,it)&
                                         + cmfd(:,ib,ik)
                     else
                        mfd1(:,ib,ik) = mfd0(:,ib,ik) + cmfd(:,ib,ik)
                     endif
                  endif
                  if(drag) mfd1(:,ib,ik) = mfd1(:,ib,ik) + drag_4e(:,ib,ik)
               enddo; enddo
               if(drag) then
                  if(present(ph_cond)) then
                     if(ionode .and. debug) write(stdout,'(12x,a)') 'Imposing Kelvin-Onsager relation'
                     call trans_impose_kelveinonsager(kg, ene, tmpr, ef, spinor, volume, -ph_cond(7:12), mfd1)
                  endif
               endif
      
               !compute conductivity from updated mfd1
               call trans_cond_calc(kg, ene, tmpr, ef, spinor, volume, mfd1, cond(:,i,it), tdf, trans_info%is_mag)
               tdf = tdf * alat * alat
               call convert_cond(cond_unit)
               if(ionode .and. drag .and. debug) write(stdout,'(12x,a6,I3,a4,12(E13.6,1x),a8)') 'e-BTE:',i,'-ITA:',cond(:,1,it)*alat*alat*cond_unit,' 1/Ohm/m'
               !save tdf to hdf5 file
               if(ionode .and. lastornot_local) call write_conf_hdf5(file_id, it, i, nestep, kg%nk, kg%numb, tdf, mfd1, &
                               trans_info%is_mag, group_id, group_id2)
               if(drag .and. debug)then
                  if(ionode) write(stdout,'(15x,a,12(E13.6,1x))')    'e-previous-cond:', cond(:,i-1,it)
                  if(ionode) write(stdout,'(15x,a,12(E13.6,1x))')    'e-current-cond: ', cond(:,i  ,it)
                  if(ionode) write(stdout,'(15x,a,12(E13.6,1x), /)') 'e-ratio-error: ', abs(cond(:,i  ,it)-cond(:,i-1,it))/abs(cond(:,i-1,it))
               endif

               !exit if converged (wthin relative error within trans_thr)
               if(converged_vec(cond(converg_vec,i-1,it), cond(converg_vec,i,it), trans_thr)) then
                  !Write the last iteration outside the iterations folder as well
                  if(ionode .and. lastornot_local) call write_conf_hdf5(file_id, it, i, nestep, kg%nk, kg%numb, tdf, mfd1,&
                                  trans_info%is_mag, group_id, group_id2, .true.)
                  exit
               endif
      
               if(i == io_nstep) then
                  if(ionode .and. lastornot_local ) call write_conf_hdf5(file_id, it, i, nestep, kg%nk, kg%numb, tdf, mfd1,&
                                  trans_info%is_mag, group_id, group_id2, .true.)
                  if(.not. drag) write(stdout,'(1x, a)') 'Warning: max number of iterations reached without convergence!'
               endif
            enddo
            !output the converged tdf in text format
            if(ionode .and. lastornot_local ) call output_tdf(ene, tdf, tmpr, ef, trans_info%is_mag, it>1)
            !output progress
            if(.not. drag)write(stdout, '(5x, i4, a, i4)') it, " /", ntemper
            flush(stdout)

         enddo
      
         if (drag) then
            mfd_F(1:3,:,:) = mfd1(4:6,:,:)
            mfd_F(4:6,:,:) = mfd1(1:3,:,:)
            e_niter = niter
         else
               write(stdout,'(5x, a)') '>finish transport.'
         endif
      
      
         if(ionode) then
            if(lastornot_local)then
               if(it_start == 1) call write_basic_hdf5(file_id, hole, niter, nestep, ene, dos)
               call hdf_close_file(file_id)
            endif
            !add the omitted factor alat**2 to cond.
            cond = cond * alat * alat
            !output conductivity and mobility
            if (drag ) then
               !compute other transport coefficient 
               call output_mobility(temper, efermi=efermi, dens=doping, cond=cond, niter_in=niter, phonon_iter=it_out, phbte=.false.)
            else
               call output_mobility(temper, efermi=efermi, dens=doping, cond=cond, niter_in=niter, phonon_iter=it_out, phbte=.false.)
               call trans_postproc()
            endif
         endif
         deallocate(ene, rates, mfd0, cond, tdf, dos)
         if(allocated(mfd1))  deallocate(mfd1)
         if(allocated(cmfd))  deallocate(cmfd)
         if(allocated(mfd_tmp)) deallocate(mfd_tmp)
         if(allocated(converg_vec)) deallocate(converg_vec)
         !
         return
      end subroutine 

      subroutine transport_ph()
         implicit none
         logical :: read_rate
         integer(HID_T) :: file_id, group_id, group_id2
         integer :: nestep, i, it, ik, ib, io_nstep, niter(ntemper), ncomp, e_niter(ntemper)
         real(dp) :: emin, emax, ef, tmpr, enk
         real(dp), allocatable :: ene(:), rates(:,:,:), cond(:,:,:), &
            mfd0(:,:,:), mfd1(:,:,:), mfd_tmp(:,:,:), imsgm(:,:,:), tdf(:,:), &
            rates_tmp(:,:), mfd_F(:,:,:), &
            drag_4ph(:,:,:), mfd_N(:,:,:), mfd_tmp2(:,:,:)
         type(grid) :: kg
         type(lattice_ifc)   :: phon   !< harmonic phonon
         type(grid) :: qg              !< q-grid
         type(lattice_ifct)  :: phont  !< anharmonic phonon
         type(electron_wann) :: elec
         integer, allocatable :: converg_vec(:)
         real(dp) :: cond_unit
      
      
         if(ionode) write(stdout,'(5x, a)') "Warn (transport_ph): the symm=.true. must hold for phonon transport calculation"
         if(ntemper .eq. 0) call errore('tranpsort_phonon', &
            'ftemper is not specified, missing temperatures and chemical potential',1)
         if(drag .or. ph_e) then
            ! Initialize electron info into elec object
            call init_electron_wann(epr_fid, kc_dim, elec)
            ! Using elec object, define properties of kgrid object kg
            ! such as velocity, energy, etc. 
            call init_boltz_grid(kg, elec, band_min, band_max)
            !Setup phonon properties inside phon object
            call init_lattice_ifc(epr_fid, qc_dim, phon)
            !Setup scattering info by computing e-ph on the fly
            !Or reading it from tmp files if it exists
            call boltz_scatter_setup(kg, elec, phon, boltz_qdim)
         endif
      
      
         !syp: if use phonon_velocity_finddiff, we should first 
         !     set up the qgrid neighbors, then calculate phonon
         !     properties elsewhere outside init_boltz_qgrid.
         call init_boltz_qgrid(qg, phon, .True.)
         if ( adapt_smear_phph .or. ph_vel_numerical ) then
            call set_grid_neighbors(qg)
         endif
         emin = minval(qg%enk)
         emax = maxval(qg%enk)
         if(emin > emax) call errore('transport-phonon','illegal phonon energy window',1)
         write(stdout,'(5x,a)') '>finish init qgrid.'
         
         !setup energy grid: emin - boltz_de : emax + boltz_de
         nestep = int( floor((emax-emin)/boltz_de_ph) ) + 3
         if(ionode) write(stdout,'(5x,a,i5,a,f12.6,a,f12.6)') &
            'Number of phonon energy steps: ', nestep, '  Emin: ', emin*ryd2ev, '  Emax: ', emax*ryd2ev
         allocate( ene(nestep) )
         do i = 1, nestep
            ene(i) = emin + (i-2)*boltz_de_ph
         enddo
        
         !syp - for phonon group velocity calculation
         call init_boltz_qgrid_vel_sym(qg, phon, ph_vel_numerical) 
         if (ionode .and. debug) call output_velocity(qg) 
      
         !read imsigma files: imaginary part of the selfenergy in meV.
         allocate( rates(qg%numb, qg%nk, ntemper), imsgm(qg%numb, qg%nk_irr, ntemper) )
         read_rate = .false.
         call load_imsigma(imsgm, read_rate)
      
         if(trans_info%is_ita .or. (.not. trans_info%is_ita .and. .not. read_rate)) then
            call init_lattice_ifct(epr_fid, qc_dim, phont, .true.)
            call boltz_ph3scat_setup(qg, phon, phont)
         endif
      
         !Compute scattering rates
         if( read_rate ) then
            !Read rates from prefix.imsigma file if it exists
            write(stdout,'(6x,a)') "- Read scattering rate from file: " // trim(prefix) // ".imsigma"
            !
            do it = 1, ntemper
            do ik = 1, qg%nk
            do ib = 1, qg%numb
               !Scattering rates in Rydberg atomic unit
               rates(ib,ik,it) = 2.0_dp*imsgm(ib, qg%kpt2ir(ik), it)
            enddo; enddo; enddo
         else
            if(.not. trans_info%is_ita) then
               write(stdout,'(6x,a)') "- imsigma file not detected for RTA - computing scattering rates on the fly"
            else
               write(stdout,'(6x,a)') "- Compute scattering rates on the fly."
            endif
            rates = 0.0_dp
            do it = 1, ntemper
               !Compute scattering rates on the fly
               call rates_ph3_scat_int(qg, temper(it), rates(:,:,it))
               if(ph_e) then
                  allocate( rates_tmp(qg%numb, qg%nk) ); rates_tmp = 0.0_dp
                  call rates_phe_scat_int(kg, temper(it), efermi(it), qg, rates_tmp) 
                  rates(:,:,it) = rates(:,:,it) + rates_tmp
                  if(allocated(rates_tmp)) deallocate(rates_tmp)
               endif
            enddo
            !output rates
            if(debug .and. ionode) then
               do it = 1, ntemper
                  call output_rates(qg, temper(it), efermi(it), rates(:,:,it), it>1)
               enddo
            endif
         endif
         deallocate( imsgm )
         if(any(rates(:,:,:) < 1.0E-9_dp)) write(stdout,'(5x, a)') "Warn (transport_ph): scattering rates less than 10^-9 a.u."
         !
         write(stdout,'(5x,a)') '>finish init rates'
      
         ! ----------------------------------------------------------------
         ! now all data we needed are ready. start to do the actual work. !
         ! ----------------------------------------------------------------
      
         ! allocate work space first.
         io_nstep = max(1, boltz_nstep)
         niter(:) = 1
         e_niter(:) = 1
      
         !ncomp = 3 :  T field
         !ncomp = 6 : Required if full_ite=true for both E and T field
         ncomp = merge(6, 3, full_ite)  !if full_ite is true, both mfd for E and T field are computed
         
         if( trans_info%is_mag) then
            call errore('transport_ph','magnetic field calculation is not implemented yet',1)
         endif
      
         if(.not. trans_info%is_mag) then
            !For a non-magnetic field calculation, need 6 elements of conductivity tensor
            allocate( cond(ncomp*2, io_nstep, ntemper) )
            allocate( mfd0(ncomp, qg%numb, qg%nk), tdf(ncomp*2, nestep))
            if(drag) allocate(mfd_F(ncomp,kg%numb,kg%nk), drag_4ph(ncomp,qg%numb,qg%nk))
            if(boltz_nstep > 0) allocate(mfd1(ncomp,qg%numb,qg%nk), mfd_tmp(ncomp,qg%numb,qg%nk), mfd_tmp2(ncomp,qg%numb,qg%nk))
         else
            call errore('transport_ph','magnetic field calculation is not implemented yet',1)
         endif

         if (converge_diag) then
            allocate(converg_vec(ncomp))
            if( ncomp == 3 ) then
               converg_vec = (/1, 3, 6/)
            elseif( ncomp == 6 ) then
               converg_vec = (/1, 3, 6/)
              !converg_vec = (/1, 3, 6, 7, 9, 12/)
            else
               call errore('transport','ncomp must be 3 or 6',1)
            endif
         else
            !use forall to traverse from 1 to ncomp
            if(trans_info%is_mag) then
               allocate(converg_vec(3*ncomp))
               forall(i=1:3*ncomp) converg_vec(i) = i
            else
               allocate(converg_vec(2*ncomp))
               forall(i=1:2*ncomp) converg_vec(i) = i
            endif
         endif
      
         !open prefix_tdf.h5 files
         !write basic information- num energy steps, energies, hole
         if(ionode) then
            call hdf_open_file(file_id, trim(prefix)//"_ph_tdf.h5", status='NEW')
         endif
         
         if (drag) then
            write(stdout,'(5x,a)') '>start coupled transport'
         elseif (ph_e) then
            write(stdout,'(5x,a)') '>start thermal transport with ph-e coupling'
         else
            write(stdout,'(5x,a)') '>start thermal transport'
         endif
         write(stdout,'(5x,a)') '>progress: (current T / total T)'
      
         cond = 0.0_dp; mfd0 = 0.0_dp
         do it = 1, ntemper
            ef = efermi(it);  tmpr = temper(it)
            !compute Fnk = tau*vnk
            do ik = 1, qg%nk
            do ib = 1, qg%numb
               ! F for uniform E-field:  \tau_nk * v_nk
               !syp - if rates != 0, then compute mfd0
               if(abs(rates(ib,ik,it)) > 1.0E-13_dp) mfd0(1:3,ib,ik) = qg%vnk(:,ib,ik)* qg%enk(ib,ik) / rates(ib, ik, it)
               if(full_ite) mfd0(4:6,ib,ik) = 0.d0 
            enddo; enddo
      
            !compute conductivity and mobility under RTA, 
            ! or starting point for iter solution
            call trans_cond_calc(qg, ene, tmpr, 0.d0, spinor, volume, mfd0, cond(:,1,it), tdf, trans_info%is_mag, .true.)
            call convert_cond(cond_unit, .true.)
            if(ionode .and. debug) write(stdout,'(7x,a7,I3,a4,12(E13.6,1x),a6)') 'ph-BTE:', 1, '-RTA:',cond(:,1,it)*alat*alat*cond_unit,' W/m/K' 
            !add the omitted factor, now tdf is fully in rydberg atomic unit
            tdf = tdf * alat * alat
            !Write down the information at this iteration into the hdf5 file
            if(ionode) call write_conf_hdf5(file_id, it, 1, nestep, qg%nk, qg%numb, tdf, mfd0, trans_info%is_mag, group_id,&
                       group_id2, trans_info%nomag_rta)
      
            !for iterative solution
            if(boltz_nstep > 0 .or. drag) mfd1(:,:,:) = mfd0(:,:,:)
      
            !compute correction from iterative solution
            do i = 2, io_nstep
      
               ! for RTA+coupled, we use io_nstep + rta to control the iterations over the coupled system
               if(drag) then
                  allocate(mfd_N(ncomp,qg%numb,qg%nk))
                  mfd_N(1:3,:,:) = mfd1(4:6,:,:)
                  mfd_N(4:6,:,:) = mfd1(1:3,:,:)
                  if(i>=4) then
                     call transport(kg, qg, it, e_niter, phon, elec, mfd_N, ph_tdf=tdf, ph_cond=cond(:,i-1,it), mfd_F=mfd_F, lastornot=.false.)
                  else
                     call transport(kg, qg, it, e_niter, phon, elec, mfd_N, ph_tdf=tdf, mfd_F=mfd_F, lastornot=.false.)
                  endif
                  call drag_scat_4ph(kg, qg, tmpr, ef, mfd_F, drag_4ph)
                  if(allocated(mfd_N)) deallocate(mfd_N)
               endif
      
               niter(it) = i
               !computed the integration terms (without the factor of tau0)
               mfd_tmp = 0.0_dp
               !If iterative calculation, compute the backscattering term
               if(trans_info%is_ita) then
                  call trans_ph3_scat_int(qg, tmpr, mfd1, mfd_tmp)
               endif
               call mp_barrier(inter_pool_comm)
      
               mfd_tmp2 = 0.0_dp
               !compute the new mean free displacement(mfd):
               do ik=1, qg%nk
               do ib=1, qg%numb
                  if(.not. trans_info%is_mag) then
                     !No B field term if a non-magnetic field calculation is performed
                     if(rates(ib,ik,it) > 1.0E-13_dp) then
                        mfd_tmp2(:,ib,ik) =  mfd_tmp(:,ib,ik)/rates(ib,ik,it)
                        mfd1(:,ib,ik) = mfd0(:,ib,ik) + mfd_tmp2(:,ib,ik)
                        if(drag) mfd1(:,ib,ik) = mfd1(:,ib,ik) + drag_4ph(:,ib,ik) / rates(ib,ik,it)
                     else
                        mfd1(:,ib,ik) = mfd0(:,ib,ik) 
                     endif
                  else
                     call errore('transport_ph','magnetic field calculation is not implemented yet',1)
                  endif
               enddo; enddo
               call trans_cond_calc(qg, ene, tmpr, 0.d0, spinor, volume, mfd1, cond(:,i,it), tdf, trans_info%is_mag, .true.)
               call convert_cond(cond_unit, .true.)
               if(ionode .and. debug) write(stdout,'(/,7x,a7,I3,a4,12(E13.6,1x),a6)') 'ph-BTE:',i,'-ITA:',cond(:,1,it)*alat*alat*cond_unit,' W/m/K'
      
               tdf = tdf * alat * alat   
               !save tdf to hdf5 file
               if(ionode) call write_conf_hdf5(file_id, it, i, nestep, qg%nk, qg%numb, tdf, mfd1, &
                               trans_info%is_mag, group_id, group_id2)
               !exit if converged (wthin relative error within trans_thr)
              !if(ionode) write(stdout,'(9x,a,12(E13.6,1x))')    'ph-previous-cond:', cond(:,i-1,it)
              !if(ionode) write(stdout,'(9x,a,12(E13.6,1x))')    'ph-current-cond: ', cond(:,i  ,it)
              !if(ionode) write(stdout,'(9x,a,12(E13.6,1x), /)') 'ph-ratio-error: ', abs(cond(:,i  ,it)-cond(:,i-1,it))/abs(cond(:,i-1,it))
               if(converged_vec(cond(converg_vec,i-1,it), cond(converg_vec,i,it), trans_thr)) then
                  !Write the last iteration outside the iterations folder as well
                  if(ionode) call write_conf_hdf5(file_id, it, i, nestep, qg%nk, qg%numb, tdf, mfd1,&
                                  trans_info%is_mag, group_id, group_id2, .true.)
                  exit
               endif
      
               if(i == io_nstep) then
                  if(ionode) call write_conf_hdf5(file_id, it, i, nestep, qg%nk, qg%numb, tdf, mfd1,&
                                  trans_info%is_mag, group_id, group_id2, .true.)
                  write(stdout,'(1x, a)') 'Warning: max number of iterations reached without convergence!'
               endif
      
            enddo
            ! for RTA+coupled, we use io_nstep + rta to control the iterations over the coupled system
            if(drag) then
               allocate(mfd_N(ncomp,qg%numb,qg%nk))
               mfd_N(1:3,:,:) = mfd1(4:6,:,:)
               mfd_N(4:6,:,:) = mfd1(1:3,:,:)
               !print mfd_N.shape
               call transport(kg, qg, it, e_niter, phon, elec, mfd_N, ph_tdf=tdf, ph_cond=cond(:,i,it), mfd_F=mfd_F, lastornot=.true.)
               if(allocated(mfd_N)) deallocate(mfd_N)
            endif
            !output progress
            write(stdout, '(5x, i4, a, i4)') it, " /", ntemper
            flush(stdout)
         enddo
         if (drag) then
            write(stdout,'(5x,a)') '>finish coupled transport'
         elseif (ph_e) then
            write(stdout,'(5x,a)') '>finish thermal transport with ph-e coupling'
         else
            write(stdout,'(5x,a)') '>finish thermal transport'
         endif
         
         if(ionode) then
            call write_basic_hdf5(file_id, hole, niter, nestep, ene)
            call hdf_close_file(file_id)
            !add the omitted factor alat**2 to cond.
            cond = cond * alat * alat
            call output_mobility(temper, cond=cond, niter_in = niter, phbte=.true.)
            call trans_postproc()
         endif
         if(allocated(ene)) deallocate(ene)
         if(allocated(rates)) deallocate(rates)
         if(allocated(mfd0)) deallocate(mfd0)
         if(allocated(cond)) deallocate(cond)
         if(allocated(tdf)) deallocate(tdf)
      
         if(allocated(mfd1))  deallocate(mfd1)
         if(allocated(mfd_tmp)) deallocate(mfd_tmp)
         if(allocated(mfd_tmp2)) deallocate(mfd_tmp2)
         if(allocated(mfd_F)) deallocate(mfd_F)
         if(allocated(drag_4ph)) deallocate(drag_4ph)
         if(allocated(converg_vec)) deallocate(converg_vec)
         !
         return
      end subroutine transport_ph
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
         real(dp), intent(in), optional :: dos_out(:)  !< Density of states at each energy step
         !
         integer :: hole_i

         call hdf_write_dataset(fid,'num_configurations', ntemper)

         call hdf_write_dataset(fid, 'energy_grid', energy(1:nstep))
         call hdf_write_attribute(fid, 'energy_grid', 'ne', nstep)

         call hdf_write_dataset(fid, 'temperature_array', temper(1:ntemper))
         call hdf_write_attribute(fid, 'temperature_array', 'nt', ntemper)

         call hdf_write_dataset(fid, 'num_iter_array', nit(1:ntemper))

         if (present(dos_out)) then
            call hdf_write_dataset(fid, 'efermi_array', efermi(1:ntemper))
            
            if(allocated(boltz_bfield)) call hdf_write_dataset(fid, 'bfield_array', boltz_bfield(:,1:ntemper))
            call hdf_write_dataset(fid, 'dos', dos_out(1:nstep))
            
            hole_i = merge(1, 0, is_hole)
            call hdf_write_attribute(fid, 'dos', 'hole', hole_i)
         endif
      end subroutine write_basic_hdf5
      !
end module mod_transport
