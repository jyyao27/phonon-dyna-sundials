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
!  main entrance for the perturbo code.
!
! Maintenance:
!===============================================================================

! Title for the Doxygen documentation
!> \mainpage PERTURBO code documentation
!! Perturbo is an open-source software package for first-principles calculations of charge transport and ultrafast carrier dynamics in solid state materials, including metals, semiconductors, insulators, and 2D materials. Visit [PERTURBO website](https://perturbo-code.github.io) for more details.

program perturbo
   use qe_mpi_mod, only: mp_startup, mp_barrier, mp_global_end, &
         inter_pool_comm, stdout, ionode, ionode_id
   use pert_param,  only: init_input_param, calc_mode, prefix, &
         dyna_phph
   use pert_param,  only: yaml_fname
   use pert_data,   only: epr_fid
   use yaml_utils,  only: open_output_yaml, close_output_yaml
   use yaml_utils_perturbo, only: output_param_to_yaml_perturbo, output_basic_data_to_yaml_perturbo
   use epr_hdf5_io, only: open_epr_file
   use hdf5_utils
   !
   use mod_transport, only: transport, transport_ph
   implicit none
   character(len=50) :: epr_filename
   logical :: has_epr_file, has_epwan_file
   
   call mp_startup()
   call environment_setup("PERTURBO", .true.)
   !readin contral parameter
   call init_input_param()
   write(stdout,'(5x,a)') "calc_mode = '" // trim(calc_mode) // "'"
   write(stdout,'(5x,a)') repeat("-",50)
   flush(stdout)

   !open HDF5 epr file and load data
   if(ionode) then
      call open_epr_file(prefix, epr_fid)
   end if

   call load_data_pert(epr_fid)

   if(ionode) then
      ! Open the YAML file
      call open_output_yaml(yaml_fname, 'perturbo')

      ! Output the input parameters to YAML
      call output_param_to_yaml_perturbo()

      ! Output the basic data from epr to YAML
      call output_basic_data_to_yaml_perturbo()
   end if

   select case (calc_mode)
      case ('bands')
      !compute wannier interpolated band structure
         call calc_band_structure()
      case ('phdisp')
      !compute phonon dispersion
         call calc_phonon_spectra()
      case ('phvel')
         call calc_phonon_velocity()
      case ('ephmat')
      ! compute |g| (both with and without 1/sqrt(w) factor)
         call calc_ephmat()
      case ('imsigma')
      !compute on-shell electron self-energy (only imaginary part): Im\Sigma(e_nk)
         call electron_imsigma()
      case ('transform-fct')
         call transform_fct()
      case ('phonlife')
      !compute on-shell phonon  self-energy (only imaginary part)
         call phonon_imsigma()
      case ('phonlife-fct')
         call phonon_imsigma_fct()
      case ('meanfp')
      !compute electron mean free path (require prefix.imsigma as input)
         call calc_electron_mfp()
      case ('setup')
      !setup a uniform k-grid for BTE dynamics or transport
         call boltz_setup()
      case ('dipole')
         call electron_dipole()
      case ('trans', 'trans-rta', 'trans-ita', 'trans-mag-rta', 'trans-mag-ita')
      !transport calculation, compute mobility using RTA, iterative(ITA)
      !magnetic RTA and ITA
      !Default value for trans is trans-rta
         call transport()
      case ('trans-ph-rta', 'trans-ph-ita')
         call transport_ph()
      case ('trans-pp')
      !compute more transport coefficients starting from the pre-calculated TDF
         if(ionode) call trans_postproc()
      case ('dynamics-run')
      !carrier dynamics, real-time dynamics by time-steping BTE. (experimental)
         call carrier_dynamics_run()
      case ('dynamics-pp')
      !postprocessing carrier dynamics simulation data. (experimental)
         call carrier_dynamics_postproc()
      case ('spectral-se')
         ! compute off-shell electron self-energy (only imaginary part): Im\Sigma(w)
         call calc_selfenergy()
      case ('spectral-re')
         ! compute real part of the on-shell electron self-enrgy: Re\Sigma(e_nk)
         call electron_resigma()
      case ('spectral-pp')
         ! for debugging, compute spectral function using retarded cumulant approach.
         call process_spectral()
      case ('spectral-cum')
         !compute electron spectral function using retarded cumulant approach
         call calc_spectral_cumulant()
      case ('spectral-trans')
         !transport calculation from spectral function using Kubo formula
         call calc_trans_cumulant()
      case ('spectral-trans-opcond')
         !compute optical conductivity
         call calc_trans_cumulant_opcond()
      case default
         write(stdout,'(1x,a,/,5x,a)') "Error: illegal calc_mode!"
   end select

   call mp_barrier(inter_pool_comm)
   if(ionode) call hdf_close_file(epr_fid)
   
   call environment_stop("PERTURBO")

   if (ionode) then
      ! Close the YAML file
      call close_output_yaml()
   endif

   call mp_global_end()
end program perturbo
