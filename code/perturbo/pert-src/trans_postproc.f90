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
!   compute more transport coefficients starting from the pre-calculated TDF
!
! Maintenance:
!===============================================================================

subroutine trans_postproc()
   use pert_const, only: dp, ryd2ev
   use pert_data,  only: volume, system_2d
   use pert_utils, only: mfermi_deriv, fermi
   use qe_mpi_mod, only: stdout, ionode
   use pert_param, only: ntemper, temper, efermi, ftemper, prefix, & 
                         calc_mode
   use boltz_trans_mod, only: extract_trans_coeff
   use boltz_trans_output, only: output_trans_coef, output_trans_coef_yaml
   use hdf5_utils
   implicit none
   logical :: has_file, is_hole, is_mag
   integer(HID_T) :: file_id, conf_id, iter_id
   integer :: ntmp, nestep, it, i, iiter, itmp, ncomp, tdf_size(2)
   integer :: ncomp2
   real(dp) :: de, tp, ef
   character(len=120) :: dset_name, fname, msg
   integer, allocatable  :: niter(:)
   real(dp), allocatable :: tmpr(:), efem(:), density(:), occup(:)
   real(dp), allocatable :: ene(:), tdf(:,:), dos(:), tdf_tmp(:,:), tdf_read(:,:)
   real(dp), allocatable :: tdf_iter(:,:,:), cond_iter(:,:,:)
   real(dp), allocatable :: cond(:,:), seebeck(:,:), therm_cond(:,:)
   real(dp), allocatable :: cond_seebeck(:), kk_coef(:), alpha(:),bfield(:,:)

   if(ntemper .eq. 0) call errore('tranpsort', &
      'ftemper is not specified, missing temperatures and chemical potential',1)

   !load data from hdf5 file
   fname = trim(prefix)//"_tdf.h5"
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('trans_postproc','missing '// trim(fname), 1)
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
   
   !read temperature
   call hdf_read_attribute(file_id, 'energy_grid', 'ne', nestep)
   call hdf_read_attribute(file_id, 'temperature_array', 'nt', ntmp)
   call hdf_read_attribute(file_id, 'dos', 'hole', itmp)
   is_hole = merge(.true., .false., itmp .eq. 1)
   allocate( tmpr(ntemper), efem(ntemper), niter(ntemper), ene(nestep), dos(nestep), &
             bfield(3,ntemper))
   bfield = 0.0_dp
   !
   write(msg,'(a)') "provided temperatures are different from those in " // trim(fname)
   if(ntmp .ne. ntemper) call errore('trans_postproc', trim(msg), 1)
   
   call hdf_read_dataset(file_id, 'temperature_array', tmpr)
   if( any(abs(tmpr-temper) > 1.0E-9_dp) ) call errore('trans_postproc', trim(msg), 1)
    
   call hdf_read_dataset(file_id, 'efermi_array', efem)
   write(msg,'(a)') "input chemical potentials are different from those in " // trim(fname)
   if( any(abs(efem - efermi) > 1.0E-6_dp) ) then
      write(stdout, '(5x, a)') "Warn (trans_postproc): "
      write(stdout, '(9x, a)') trim(msg)
   endif
   !
   deallocate( tmpr, efem )

   !If magnetic field is true, read its values
   is_mag = hdf_exists(file_id,'bfield_array')
   if(is_mag) call hdf_read_dataset(file_id, 'bfield_array', bfield)
   call hdf_read_dataset(file_id, 'dos', dos)
   call hdf_read_dataset(file_id, 'num_iter_array', niter)
   call hdf_read_dataset(file_id, 'energy_grid', ene)

   !Read tdf size (Maybe change it later?)
   call hdf_open_group(file_id, 'configuration_1', conf_id)
   call hdf_get_dims(conf_id, 'tdf', tdf_size)
   call hdf_close_group(conf_id)

   de = (ene(nestep)-ene(1)) / (nestep-1.0_dp)
   ncomp = tdf_size(1)

   !The variable ncomp2 is for the do loops in the next lines
   !Since we need to loop over 9 elements for non-zero magnetic field
   ncomp2 = merge(9,6,is_mag)

   if(ncomp.ne.6 .and. ncomp.ne.12 .and. ncomp.ne.9)&
   call errore('trans_postproc','wrong dimension in the tdf data',1)
   !
   allocate( tdf_tmp(ncomp2, nestep), tdf(ncomp, nestep), occup(nestep), density(ntemper) )
   allocate( cond_seebeck(ncomp2), kk_coef(ncomp2), alpha(ncomp2))
   allocate( cond(ncomp2, ntemper), seebeck(ncomp2, ntemper), therm_cond(ncomp2, ntemper) )
   cond = 0.0_dp; seebeck = 0.0_dp; therm_cond = 0.0_dp;

   allocate(cond_iter(ncomp2, maxval(niter), ntemper))
   allocate(tdf_read(ncomp, nestep))
   cond_iter = 0.0_dp

   do it = 1, ntemper
      tp = temper(it) !temperature
      ef = efermi(it) !chemical potential
      !load tdf in Rydberg atomic unit
      ! TDF_ij(E) = spinor/(N*V) sum_{nk} vnk*vnk*tau*delta(E-e_nk)
      call read_tdf_hdf5(file_id, it, niter(it), nestep, tdf, .true.)
      
      !compute carrier concentration: #. carrier / bohr^3
      occup(:) = fermi(tp, merge(ef-ene(:), ene(:)-ef, is_hole) )
      density(it) = dot_product(dos, occup) * de / volume
      
      !compute conductivity
      do i = 1, ncomp2
         tdf_tmp(i,:)  = tdf(i,:) * mfermi_deriv(tp, (ene(:) - ef))
      enddo
      !\sigma = q^2 * Int tdf(E) * (-df/dE) * dE
      cond(1:ncomp2, it) = sum(tdf_tmp, dim=2) * de

      ! for the YAML output, compute the onductibity for all iterations
      do iiter = 1, niter(it)

         call read_tdf_hdf5(file_id, it, iiter, nestep, tdf_read, iiter .eq. niter(it))

         do i = 1, ncomp2
            tdf_tmp(i,:)  = tdf_read(i,:) * mfermi_deriv(tp, (ene(:) - ef))
         enddo

         cond_iter(1:ncomp2, iiter, it) = sum(tdf_tmp, dim=2) * de

      enddo ! niter

      if(.not. is_mag) then
         ! alpha corresponding to T * \sigma * seebeck, using tdf of E-field
         !  = - q Int tdf(E) * (-df/dE) * (E-\mu) * dE   (Eq.7c in PRB 94, 085204, 2016)
         do i = 1, ncomp2
            tdf_tmp(i,:) = tdf_tmp(i,:) * (ene(:) - ef)
         enddo
         alpha(1:ncomp2) = - sum(tdf_tmp, dim=2) * de
      
         !
         ! \sigma * seebeck = - q/T Int tdf(E) * (-df/dE) * (E-\mu) * dE
         if(ncomp > 6) then
            !compute \sigma*seebeck, using tdf of T-field (note the minus sign below)
            do i = 1, 6
               tdf_tmp(i,:)  = tdf(i+6,:) * mfermi_deriv(tp, (ene(:) - ef))
            enddo
            cond_seebeck(1:6) = - sum(tdf_tmp, dim=2) * ( de / tp )
         else
            !!\sigma * seebeck = - q/T Int tdf(E) * (-df/dE) * (E-\mu) * dE
            !!cond_seebeck(1:6) = - sum(tdf_tmp, dim=2) * ( de / tp )
            cond_seebeck = alpha / tp
         endif

         !compute K for thermal conductivity, K = 1/T Int tdf(E) * (-df/dE) * (E-\mu)^2 * dE
         ! (using tdf of T-field (beta in Eq.7d PRB 94,085204,2016) if available.)
         do i = 1, ncomp2
            tdf_tmp(i,:) = tdf_tmp(i,:) * (ene(:) - ef)
         enddo
         kk_coef(1:ncomp2) = sum(tdf_tmp, dim=2) * ( de / tp )

         ! seebeck = sigma^-1 * [sigma*seebeck]
         ! 1, K_B/q is omitted, q is electron charge and K_b is boltzmann constant 
         ! 2, the vlue of seebeck here is dimensionless.
         !
         ! thermal conductivity \kappa = K - alpha * seebeck (or K - T*[sigma*seebeck]*seebeck)
         call extract_trans_coeff(tp, cond(:,it), cond_seebeck, kk_coef, system_2d, &
            seebeck(:,it), therm_cond(:,it), ncomp2, alpha)
      endif
   enddo
   call output_trans_coef(temper, efermi, density, cond, seebeck, therm_cond, is_mag)

   call output_trans_coef_yaml(temper, efermi, density, cond, cond_iter, niter, seebeck, therm_cond, bfield, is_mag)
   deallocate(tdf_read)

   deallocate( tdf_tmp, tdf, occup, density, cond, seebeck, therm_cond )

   call hdf_close_file(file_id)

   return

   contains

      subroutine read_tdf_hdf5(fid, tid, iter, nstep, tdf, last_iter)
         use pert_const, only: dp
         use HDF5_utils
         implicit none
         integer(HID_T), intent(in) :: fid
         logical, intent(in), optional :: last_iter 
         integer, intent(in) :: tid, iter, nstep
         real(dp), intent(out) :: tdf(:,:)
         !
         character(len=120) :: dset_name
         character(len=6), external :: int_to_char
         logical :: is_last_iter
         integer(HID_T) :: iter_id, conf_id

         is_last_iter = .false.
         if(present(last_iter)) is_last_iter = last_iter

         call hdf_open_group(fid, 'configuration_'//trim(int_to_char(tid)), conf_id)

         if(.not. is_last_iter) then
            call hdf_open_group(conf_id, 'iterations', iter_id)
            dset_name = 'tdf_' // trim(int_to_char(iter))
            call hdf_read_dataset(iter_id, trim(dset_name), tdf(:, 1:nstep))
            call hdf_close_group(iter_id)
         else
            call hdf_read_dataset(conf_id, 'tdf', tdf(:, 1:nstep))
         endif

         call hdf_close_group(conf_id)

      end subroutine read_tdf_hdf5
      !
end subroutine trans_postproc
