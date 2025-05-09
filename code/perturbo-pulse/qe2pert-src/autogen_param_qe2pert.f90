! This file was automatically generated by the f90tools.py Python script
! from the ./utils folder.
! To do any modifications, please modify the YAML files in the ./docs folder or the script directly.
! NOTE THAT the modifications will be erased when you run the Python script again.
! Date: April 09, 2025 15:55

module qe2pert_autogen_param
   use kinds, only: dp
   use io_files, only: prefix
   implicit none
   character(len=10), save :: asr !< Indicates the type of Acoustic Sum Rule imposed.
   logical, save :: debug !< Set to .true. to turn on the debug mode, in which the code stop after g(k,q) (does not compute g in wannier basis)
   integer, save :: dft_band_max !< Highest band index used in Wannier90. Be default, it will be reset to the highest band index in the DFT results.
   integer, save :: dft_band_min !< Lowest band index used in Wannier90.
   real(dp), save :: dis_win_min !< The 'dis_win_min' used in Wannier90, the lower boundary of the outer windows.
   character(len=256), save :: eig_corr !< File containing the electron eigenvalues on the (nk1, nk2, nk3) grid. The format of this file is the same as the file <code>prefix.eig</code> generated in Wannier90. if present, <code>qe2pert.x</code> will read the eigenvalues from this file, rather than Kohn-Sham eigenvalues from QE-nscf calculation. This is usually used when one wants to use modified eigenvalues (e.g., from GW).
   logical, save :: load_ephmat !< Set to .true. to load prefix_elph.h5 from the directory specified by the variable outdir.
   logical, save :: lwannier !< Set to .true. to rotate the wavefunctions using Wannier unitary matrix before computing e-ph matrix elements.
   integer :: nk1 !< Number of k points along x-axis used in the Wannierization.
   integer :: nk2 !< Number of k points along y-axis used in the Wannierization.
   integer :: nk3 !< Number of k points along z-axis used in the Wannierization.
   integer, save :: num_wann !< Number of Wannier functions.
   character(len=256) :: outdir !< Name of the directory where the QE nscf output directory prefix.save is located, and where the e-ph matrix elements prefix_elph.h5 will be stored.
   character(len=256), save :: phdir !< Name of the directory where the phonon "save" directory is located.
   real(dp), save :: polar_alpha !< Convergence parameter used in the Ewald sum when computing the polar correction in polar materials. The default value is 1.0.
   character(len=4), save :: spin_component !< Flag for LSDA spin polarized calculation.
   logical, save :: system_2d !< Set it to .true. if the system is 2D.
   logical, save :: tdep !< Flag for TDEP support.
   real(dp), save :: thickness_2d !< Thickness of the 2d system, used in the 2D polar e-ph correction. Only needed when system_2d=.true.
   character(len=80) :: yaml_fname !< Name of the YAML output file.

   character(len=10), save :: asr_beforeconv !< Indicates the type of Acoustic Sum Rule imposed. (before conversion of the unit)
   logical, save :: debug_beforeconv !< Set to .true. to turn on the debug mode, in which the code stop after g(k,q) (does not compute g in wannier basis) (before conversion of the unit)
   integer, save :: dft_band_max_beforeconv !< Highest band index used in Wannier90. Be default, it will be reset to the highest band index in the DFT results. (before conversion of the unit)
   integer, save :: dft_band_min_beforeconv !< Lowest band index used in Wannier90. (before conversion of the unit)
   real(dp), save :: dis_win_min_beforeconv !< The 'dis_win_min' used in Wannier90, the lower boundary of the outer windows. (before conversion of the unit)
   character(len=256), save :: eig_corr_beforeconv !< File containing the electron eigenvalues on the (nk1, nk2, nk3) grid. The format of this file is the same as the file <code>prefix.eig</code> generated in Wannier90. if present, <code>qe2pert.x</code> will read the eigenvalues from this file, rather than Kohn-Sham eigenvalues from QE-nscf calculation. This is usually used when one wants to use modified eigenvalues (e.g., from GW). (before conversion of the unit)
   logical, save :: load_ephmat_beforeconv !< Set to .true. to load prefix_elph.h5 from the directory specified by the variable outdir. (before conversion of the unit)
   logical, save :: lwannier_beforeconv !< Set to .true. to rotate the wavefunctions using Wannier unitary matrix before computing e-ph matrix elements. (before conversion of the unit)
   integer :: nk1_beforeconv !< Number of k points along x-axis used in the Wannierization. (before conversion of the unit)
   integer :: nk2_beforeconv !< Number of k points along y-axis used in the Wannierization. (before conversion of the unit)
   integer :: nk3_beforeconv !< Number of k points along z-axis used in the Wannierization. (before conversion of the unit)
   integer, save :: num_wann_beforeconv !< Number of Wannier functions. (before conversion of the unit)
   character(len=256) :: outdir_beforeconv !< Name of the directory where the QE nscf output directory prefix.save is located, and where the e-ph matrix elements prefix_elph.h5 will be stored. (before conversion of the unit)
   character(len=256), save :: phdir_beforeconv !< Name of the directory where the phonon "save" directory is located. (before conversion of the unit)
   real(dp), save :: polar_alpha_beforeconv !< Convergence parameter used in the Ewald sum when computing the polar correction in polar materials. The default value is 1.0. (before conversion of the unit)
   character(len=256) :: prefix_beforeconv !< Job name prefix. It should be the same as the prefix used in QE. (before conversion of the unit)
   character(len=4), save :: spin_component_beforeconv !< Flag for LSDA spin polarized calculation. (before conversion of the unit)
   logical, save :: system_2d_beforeconv !< Set it to .true. if the system is 2D. (before conversion of the unit)
   logical, save :: tdep_beforeconv !< Flag for TDEP support. (before conversion of the unit)
   real(dp), save :: thickness_2d_beforeconv !< Thickness of the 2d system, used in the 2D polar e-ph correction. Only needed when system_2d=.true. (before conversion of the unit)
   character(len=80) :: yaml_fname_beforeconv !< Name of the YAML output file. (before conversion of the unit)

   namelist / qe2pert / & 
      asr, & 
      debug, & 
      dft_band_max, & 
      dft_band_min, & 
      dis_win_min, & 
      eig_corr, & 
      load_ephmat, & 
      lwannier, & 
      nk1, & 
      nk2, & 
      nk3, & 
      num_wann, & 
      outdir, & 
      phdir, & 
      polar_alpha, & 
      prefix, & 
      spin_component, & 
      system_2d, & 
      tdep, & 
      thickness_2d, & 
      yaml_fname

contains
subroutine autogen_init_input_param()
   implicit none
   asr = 'crystal'
   debug = .false.
   dft_band_max = 10000
   dft_band_min = 1
   dis_win_min = -9999.0_dp
   eig_corr = ''
   load_ephmat = .false.
   lwannier = .true.
   num_wann = 1
   polar_alpha = 1.0_dp
   spin_component = 'none'
   system_2d = .false.
   tdep = .false.
   thickness_2d = 6.0_dp
   yaml_fname = 'qe2pert_output.yml'

end subroutine autogen_init_input_param

subroutine autogen_bcast_input_param()
   use io_global, only: ionode_id
   use mp_world, only: world_comm
   use mp, only: mp_bcast
   implicit none

   call mp_bcast(asr,ionode_id,world_comm)
   call mp_bcast(debug,ionode_id,world_comm)
   call mp_bcast(dft_band_max,ionode_id,world_comm)
   call mp_bcast(dft_band_min,ionode_id,world_comm)
   call mp_bcast(dis_win_min,ionode_id,world_comm)
   call mp_bcast(eig_corr,ionode_id,world_comm)
   call mp_bcast(load_ephmat,ionode_id,world_comm)
   call mp_bcast(lwannier,ionode_id,world_comm)
   call mp_bcast(nk1,ionode_id,world_comm)
   call mp_bcast(nk2,ionode_id,world_comm)
   call mp_bcast(nk3,ionode_id,world_comm)
   call mp_bcast(num_wann,ionode_id,world_comm)
   call mp_bcast(outdir,ionode_id,world_comm)
   call mp_bcast(phdir,ionode_id,world_comm)
   call mp_bcast(polar_alpha,ionode_id,world_comm)
   call mp_bcast(prefix,ionode_id,world_comm)
   call mp_bcast(spin_component,ionode_id,world_comm)
   call mp_bcast(system_2d,ionode_id,world_comm)
   call mp_bcast(tdep,ionode_id,world_comm)
   call mp_bcast(thickness_2d,ionode_id,world_comm)
   call mp_bcast(yaml_fname,ionode_id,world_comm)

end subroutine autogen_bcast_input_param

subroutine autogen_input_param_beforeconv()
   implicit none
   asr_beforeconv = asr
   debug_beforeconv = debug
   dft_band_max_beforeconv = dft_band_max
   dft_band_min_beforeconv = dft_band_min
   dis_win_min_beforeconv = dis_win_min
   eig_corr_beforeconv = eig_corr
   load_ephmat_beforeconv = load_ephmat
   lwannier_beforeconv = lwannier
   nk1_beforeconv = nk1
   nk2_beforeconv = nk2
   nk3_beforeconv = nk3
   num_wann_beforeconv = num_wann
   outdir_beforeconv = outdir
   phdir_beforeconv = phdir
   polar_alpha_beforeconv = polar_alpha
   prefix_beforeconv = prefix
   spin_component_beforeconv = spin_component
   system_2d_beforeconv = system_2d
   tdep_beforeconv = tdep
   thickness_2d_beforeconv = thickness_2d
   yaml_fname_beforeconv = yaml_fname

end subroutine autogen_input_param_beforeconv

end module qe2pert_autogen_param
