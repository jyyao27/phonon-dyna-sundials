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

program qe2pert
   use mp_global, only: mp_startup, mp_global_end
   use qe_mpi_mod, only: ionode
   use input_param, only: read_input_param, load_ephmat, debug, yaml_fname, &
               tdep
   use electronic_data, only: init_electronic_data, deallocate_evc_bec
   use lattice_data, only: init_lattice_data
   use elph_matrix, only: compute_elph_matrix
   use perturbo_data, only: save_perturbo_data
   use yaml_utils,  only: open_output_yaml, close_output_yaml
   use yaml_utils_qe2pert, only: output_param_to_yaml_qe2pert

   implicit none

   call mp_startup()
   call environment_setup("QE2PERT", .false.)
      
   call read_input_param()
 
   
   ! first check if tdep flag is true and if tdep is compiled
   if (tdep) then
#ifndef __TDEP
   call errore('qe2pert', "Cannot save TDEP data. Recompile Perturbo with -D__TDEP.", 1)
#endif
   endif

   ! Open the YAML file
   if(ionode) then
      call open_output_yaml(yaml_fname, 'qe2pert')
      call output_param_to_yaml_qe2pert()
   end if 

   ! read prefix.save from pwscf calculation.
   call init_electronic_data()
   ! read lattice dynamical from PHonon calculations
   call init_lattice_data() 


   call save_perturbo_data()
      
   if( load_ephmat ) then
      call save_elph_mat_wann()
   else
      call compute_elph_matrix()
      !release large arrays that will not be used.
      call deallocate_evc_bec()
      !output data to hdf5 file
      call collect_elph_matrix()
      ! compute and save g(Re, Rp)
      if(.not. debug) call save_elph_mat_wann()
   endif

   call environment_stop("QE2PERT")
   
   ! Close the YAML file
   if(ionode) then
      call close_output_yaml()
   end if
   call mp_global_end()
end program qe2pert
