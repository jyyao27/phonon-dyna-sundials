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
!  read parameters from the input file
!
! Maintenance:
!===============================================================================

module input_param
   use kinds, only: dp
   use constants, only: angstrom_au, RYTOEV
   use io_files,  only: prefix, tmp_dir
   use io_global, only: stdout, stdin, ionode, ionode_id
   use mp_world, only: world_comm, root
   use mp, only: mp_bcast, mp_sum
   use qe2pert_autogen_param
   implicit none
   public

   !dimension of the regular k-grid.
   integer, save :: kdim(3)

   integer, save :: num_band

   public :: read_input_param
contains

subroutine read_input_param()
! read input from stdin
   implicit none
   ! local variable
   logical :: has_file
   CHARACTER(LEN=256) :: fname
   integer :: ios

   CHARACTER(LEN=256), EXTERNAL :: trimcheck

   call autogen_init_input_param()
   !
   ! master read
   if (ionode) then
      ! read input file
      call input_from_file()
      read(stdin, nml=qe2pert, err=100, iostat=ios)

      kdim(1) = nk1
      kdim(2) = nk2
      kdim(3) = nk3
   end if   

   ! keep the variables before unit conversion
   call autogen_input_param_beforeconv()

   ! broadcast
   call autogen_bcast_input_param()

   call mp_bcast(kdim, ionode_id, world_comm)
   tmp_dir = trimcheck(outdir)
   phdir = trimcheck(phdir)
   eig_corr = trim(adjustl(eig_corr))
   ! covert to rydberg atomic unit
   thickness_2d = thickness_2d * angstrom_au
   dis_win_min = dis_win_min / RYTOEV
   
   !sanity check
   if( polar_alpha < 1.0E-12 ) &
      call errore('input_param','polar_alpha too small !!!',1)
   
   if( system_2d .and. thickness_2d <= 0.0d0 ) then
      call errore('input_param','thickness_2d is negative (or too small) !!!',1)
   elseif (.not. system_2d) then
      thickness_2d = -6.0_dp  ! negative value means 3D system
   endif

   if( dft_band_min < 1 ) &
      call errore('input_param','dft_band_min shuld be larger than 0',1)
   
   if( load_ephmat ) then
      fname = trim(tmp_dir) // trim(prefix) // "_elph.h5"
      inquire(file=trim(fname), exist=has_file)
      if(.not. has_file) call errore('input_param', "Missing file: "//trim(fname), 1)
   endif
   
   if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'no')) then
      call errore('input_param','invalid Acoustic Sum Rule :' // asr, 1)
   endif

#if defined(__TDEP)
   if(tdep .and. system_2d) write(stdout,'(2x,a)') &
      "Warning: TDEP for 2D polar materials is not tested, use with caution!!"
#else
   if(tdep) call errore('input_param',"tdep=.true. requires Perturbo compiled with -D__TDEP",1)
#endif

   return
   100 call errore('init_input_parameters','reading input namelist',abs(ios))
end subroutine read_input_param

end module input_param
