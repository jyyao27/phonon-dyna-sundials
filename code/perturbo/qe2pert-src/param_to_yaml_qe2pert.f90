! This file was automatically generated by the f90tools.py Python script
! from the ./utils folder.
! To do any modifications, please modify the YAML files in the ./docs folder or the script directly.
! NOTE THAT the modifications will be erased when you run the Python script again.
! Date: February 16, 2024 10:29

module qe2pert_autogen_output_yaml
   use yaml_utils, only: ymlout, python_bool
   use input_param
   implicit none
   private

   public :: auto_output_beforeconv_to_yaml
   public :: auto_output_afterconv_to_yaml

contains
subroutine auto_output_afterconv_to_yaml()
   implicit none
   write(ymlout,'(6x, a, 10x, a)') 'asr:', trim(asr)
   write(ymlout,'(6x, a, 10x, a)') 'debug:', python_bool(debug)
   write(ymlout,'(6x, a, 10x, i10)') 'dft_band_max:', dft_band_max
   write(ymlout,'(6x, a, 10x, i10)') 'dft_band_min:', dft_band_min
   write(ymlout,'(6x, a, 10x, es23.16)') 'dis_win_min:', dis_win_min
   write(ymlout,'(6x, a, 10x, a)') 'eig_corr:', trim(eig_corr)
   write(ymlout,'(6x, a, 10x, a)') 'load_ephmat:', python_bool(load_ephmat)
   write(ymlout,'(6x, a, 10x, a)') 'lwannier:', python_bool(lwannier)
   write(ymlout,'(6x, a, 10x, i10)') 'nk1:', kdim(1)
   write(ymlout,'(6x, a, 10x, i10)') 'nk2:', kdim(2)
   write(ymlout,'(6x, a, 10x, i10)') 'nk3:', kdim(3)
   write(ymlout,'(6x, a, 10x, i10)') 'num_wann:', num_wann
   write(ymlout,'(6x, a, 10x, a)') 'outdir:', trim(tmp_dir)
   write(ymlout,'(6x, a, 10x, a)') 'phdir:', trim(phdir)
   write(ymlout,'(6x, a, 10x, es23.16)') 'polar_alpha:', polar_alpha
   write(ymlout,'(6x, a, 10x, a)') 'prefix:', trim(prefix)
   write(ymlout,'(6x, a, 10x, a)') 'spin_component:', trim(spin_component)
   write(ymlout,'(6x, a, 10x, a)') 'system_2d:', python_bool(system_2d)
   write(ymlout,'(6x, a, 10x, a)') 'tdep:', python_bool(tdep)
   write(ymlout,'(6x, a, 10x, es23.16)') 'thickness_2d:', thickness_2d
   write(ymlout,'(6x, a, 10x, a)') 'yaml_fname:', trim(yaml_fname)

end subroutine auto_output_afterconv_to_yaml

! the "before conversion" variabls has not been broadcast to other nodes, so it can just be used in master (io) node.
subroutine auto_output_beforeconv_to_yaml()
   implicit none
   write(ymlout,'(6x, a, 10x, a)') 'asr:', trim(asr_beforeconv)
   write(ymlout,'(6x, a, 10x, a)') 'debug:', python_bool(debug_beforeconv)
   write(ymlout,'(6x, a, 10x, i10)') 'dft_band_max:', dft_band_max_beforeconv
   write(ymlout,'(6x, a, 10x, i10)') 'dft_band_min:', dft_band_min_beforeconv
   write(ymlout,'(6x, a, 10x, es23.16)') 'dis_win_min:', dis_win_min_beforeconv
   write(ymlout,'(6x, a, 10x, a)') 'eig_corr:', trim(eig_corr_beforeconv)
   write(ymlout,'(6x, a, 10x, a)') 'load_ephmat:', python_bool(load_ephmat_beforeconv)
   write(ymlout,'(6x, a, 10x, a)') 'lwannier:', python_bool(lwannier_beforeconv)
   write(ymlout,'(6x, a, 10x, i10)') 'nk1:', nk1_beforeconv
   write(ymlout,'(6x, a, 10x, i10)') 'nk2:', nk2_beforeconv
   write(ymlout,'(6x, a, 10x, i10)') 'nk3:', nk3_beforeconv
   write(ymlout,'(6x, a, 10x, i10)') 'num_wann:', num_wann_beforeconv
   write(ymlout,'(6x, a, 10x, a)') 'outdir:', trim(outdir_beforeconv)
   write(ymlout,'(6x, a, 10x, a)') 'phdir:', trim(phdir_beforeconv)
   write(ymlout,'(6x, a, 10x, es23.16)') 'polar_alpha:', polar_alpha_beforeconv
   write(ymlout,'(6x, a, 10x, a)') 'prefix:', trim(prefix_beforeconv)
   write(ymlout,'(6x, a, 10x, a)') 'spin_component:', trim(spin_component_beforeconv)
   write(ymlout,'(6x, a, 10x, a)') 'system_2d:', python_bool(system_2d_beforeconv)
   write(ymlout,'(6x, a, 10x, a)') 'tdep:', python_bool(tdep_beforeconv)
   write(ymlout,'(6x, a, 10x, es23.16)') 'thickness_2d:', thickness_2d_beforeconv
   write(ymlout,'(6x, a, 10x, a)') 'yaml_fname:', trim(yaml_fname_beforeconv)

end subroutine auto_output_beforeconv_to_yaml

end module qe2pert_autogen_output_yaml