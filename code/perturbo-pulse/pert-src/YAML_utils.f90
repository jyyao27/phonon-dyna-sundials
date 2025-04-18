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

module yaml_utils
   use kinds, only: dp
   use sys_utils, only: c_to_f_string, find_proc_file_value, bytes_to_gib

   integer, protected :: ymlout

   public :: open_output_yaml, print_clock_yaml, output_tensor_yaml
   public :: output_grid_yaml, python_bool

contains

subroutine open_output_yaml(fname, prog_name)
   implicit none
   character(len=*), intent(inout) :: fname
   character(len=*), intent(in) :: prog_name

   fname = adjustl(fname)
   open(newunit = ymlout, file = fname, action='write')
   write(ymlout, '(a)') '---'

   write(ymlout,'(a,a)') 'program: ', trim(prog_name)

   call date_and_time_yaml()

   ! For stdout, parallel info is called in the environment_setup subroutine
   ! However, environment setup happens before init_input_param,
   ! therefore, at that moment yaml_fname is not initialized yet. 
   ! Hence, we call parallel_info_yaml here.
   call parallel_info_yaml()

end subroutine open_output_yaml

subroutine close_output_yaml()
   use qe_mpi_mod, only: ionode
   implicit none

   call memory_info_yaml()

   if( ionode ) then
      call print_clock_yaml(' ')
   endif

   write(ymlout, '(a)') '...'
   close(ymlout)

end subroutine close_output_yaml

! Adapted from Modules/date_and_tim.f90
subroutine date_and_time_yaml()
   use qe_mpi_mod, only: ionode
   implicit none
   character(len=50)  :: cdate, ctime
   character(len=15), dimension(12) :: months
   data months /'January','February','March','April','May','June', &
        'July','August','September','October','November','December'/
   integer date_time(8)

   if( ionode ) then

      call date_and_time(values=date_time)

      write(cdate,'(i2,1x,a,1x,i4)') date_time(3), trim(months(date_time(2))),&
                                     date_time(1)
      write(ctime,'(i2,":",i2,":",i2)') date_time(5), date_time(6), date_time(7)

      write(ymlout,'(/,a)') 'start date and time:'
      write(ymlout,'(3x,a,1x,a)') 'date:', trim(cdate)
      write(ymlout,'(3x,a,1x,3a)') 'time:', '"', trim(ctime), '"'

   endif

end subroutine date_and_time_yaml


#if defined(_OPENACC)
subroutine acc_info_yaml()
   use openacc
#if defined(SHOW_NVML_INFO)
   use nvgpu, only: gpu_uuid, gpu_meminfo
   use iso_c_binding, only: c_char
#endif
   implicit none

   integer :: num_gpus
   integer, value :: orig_dev_num, dev_num
   integer(acc_device_kind), value :: dev_type

   character(len=80) :: str_name, str_driver
   integer(c_size_t), value :: device_mem, device_free_mem

#if defined(SHOW_NVML_INFO)
   ! If we have an NVIDIA GPU:
   character(kind=c_char, len=96) :: str_uuid
   integer(c_int) :: rc
   integer(c_long_long) :: free_mem, reserved_mem, total_mem, used_mem
#endif

   dev_type = acc_get_device_type()
   num_gpus = acc_get_num_devices(dev_type)
   orig_dev_num = acc_get_device_num(dev_type)

   write(ymlout, '(3x,a,i4)') 'number of OpenACC devices:', num_gpus
   write(ymlout, '(3x,a,i4)') 'OpenACC device type:', dev_type

   ! write (*,'(/5X,I0," ACC devices of type ",I0)') num_gpus, dev_type
   do dev_num = 0, num_gpus - 1
      call acc_set_device_num(dev_num, dev_type)
      call acc_get_property_string(dev_num, dev_type, acc_property_name, str_name)
      call acc_get_property_string(dev_num, dev_type, acc_property_driver, str_driver)

      device_mem = acc_get_property(dev_num, dev_type, acc_property_memory)
      device_free_mem = acc_get_property(dev_num, dev_type, acc_property_free_memory)

      str_name = c_to_f_string(str_name)
      str_driver = c_to_f_string(str_driver)

      write(ymlout, '(3x,"OpenACC device ",i2,": " ,a," (driver ",a, ") - ", F5.2," GiB memory,", F5.2," GiB free")') &
         dev_num, trim(str_name), trim(str_driver), bytes_to_gib(device_mem), bytes_to_gib(device_free_mem)

#if defined(SHOW_NVML_INFO)
      ! If we have an NVIDIA GPU:

      write(ymlout, '(3x,"OpenACC device ",i2," NVML info:")') dev_num
      rc = gpu_uuid(dev_num, str_uuid, 96)
      if (rc /= 0) str_uuid = "unknown"
      write (ymlout, '(6x,a,a)') 'uuid: ', c_to_f_string(str_uuid)
 
      rc = gpu_meminfo(dev_num, free_mem, reserved_mem, total_mem, used_mem)
      if (rc == 0) then
         write (ymlout, '(6x,a,f5.2,a)') 'mem free: ', bytes_to_gib(free_mem), ' GiB'
         write (ymlout, '(6x,a,f5.2,a)') 'mem total: ', bytes_to_gib(total_mem), ' GiB'
         write (ymlout, '(6x,a,f5.2,a)') 'mem reserved: ', bytes_to_gib(reserved_mem), ' GiB'
         write (ymlout, '(6x,a,f5.2,a)') 'mem used: ', bytes_to_gib(used_mem), ' GiB'
      endif
#endif
   enddo

   write(ymlout, '(3x,a,i2)') 'using OpenACC device: ', orig_dev_num

   ! Restore the original device number
   call acc_set_device_num(orig_dev_num, dev_type)
end subroutine acc_info_yaml
#endif


subroutine parallel_info_yaml()
   use mp_world,  only: nproc, nnode
   implicit none
   character(len=40) :: apis
#if defined(_OPENMP)
   integer, external :: omp_get_max_threads
#endif

   write(ymlout, '(/,a)') 'parallelization:'

   ! Possible results:
   !   MPI
   !   MPI+OpenMP
   !   MPI+OpenACC
   !   MPI+OpenMP+OpenACC
   !   OpenMP
   !   OpenMP+OpenACC
   !   OpenACC
   !   serial
#if defined(__MPI) || defined(_OPENMP) || defined(_OPENACC)

   apis = ""
#if defined(__MPI)
   apis = "+MPI"
#endif

#if defined(_OPENMP)
   apis = trim(apis) // "+OpenMP"
#endif

#if defined(_OPENACC)
   apis = trim(apis) // "+OpenACC"
#endif
   ! Chop off the leading "+" character.  We can assume the string has at
   ! least one character since the outer #if ensures that at least one of
   ! the inner #if guards will run.
   apis = apis(2:)

#else
   ! No parallelization
   apis = "serial"
#endif

   write(ymlout, '(3x,a)') 'type: ' // trim(apis)

#if defined(__MPI)
   write(ymlout, '(3x,a,i4)') 'number of MPI tasks:', nproc
#endif

#if defined(_OPENMP)
   write(ymlout, '(3x,a,i4)') 'number of OpenMP threads:', omp_get_max_threads()
#endif

#if defined(_OPENACC)
   call acc_info_yaml()  ! Output OpenACC info
#endif

#if (!defined(__GFORTRAN__) ||  ((__GNUC__>4) || ((__GNUC__==4) && (__GNUC_MINOR__>=8))) ) && defined(__MPI)
    write(ymlout, '(3x,a,i4)') 'number of nodes:', nnode
#endif

end subroutine parallel_info_yaml


subroutine memory_info_yaml()
#if defined(_OPENACC)
   use openacc
#endif
   implicit none

#if defined(_OPENACC)
   integer(acc_device_kind), value :: dev_type
   integer, value :: dev_num
   integer(c_size_t), value :: device_mem, device_free_mem
#endif

   ! These memory statistics are only accessible via /proc on Linux, so
   ! guard this code.  Other platforms can be added in time if desired.
#if defined(__linux__)
   write(ymlout, '(/,a)') 'memory:'

   write(ymlout, '(3x,a,a)') 'physical: ', find_proc_file_value('/proc/meminfo', 'MemTotal', 'unknown')
   write(ymlout, '(3x,a,a)') 'max rss: ', find_proc_file_value('/proc/self/status', 'VmHWM', 'unknown')
#endif

#if defined(_OPENACC)
   dev_type = acc_get_device_type()
   dev_num = acc_get_device_num(dev_type)
   device_mem = acc_get_property(dev_num, dev_type, acc_property_memory)
   device_free_mem = acc_get_property(dev_num, dev_type, acc_property_free_memory)

   write(ymlout, '(3x,a,f5.2," GiB")') 'gpu in-use: ', bytes_to_gib(device_mem - device_free_mem)
#endif

end subroutine memory_info_yaml


subroutine output_grid_yaml(num_tet, nk, nk_irr, g)
   implicit none
   integer, intent(in) :: num_tet, nk, nk_irr
   character(len=1), intent(in) :: g ! grid type 'k' or 'q'

   write(ymlout, '(/,a)') g//'-grid:'
   write(ymlout,'(3x,a,i11 )') 'number of  tetrahedra selected:', num_tet
   write(ymlout,'(3x,a,i11 )') 'number of reducible ' //g //'-points:', nk
   write(ymlout,'(3x,a,i11/)') 'number of irreducible ' // g //'-points:', nk_irr

end subroutine

subroutine output_scatter_info_yaml(npool, nq_col, collect_pools)
   use qe_mpi_mod, only: ionode
   implicit none
   integer, intent(in) :: npool, nq_col(npool), collect_pools(npool)
   !
   integer :: i
   character(len=6), external :: int_to_char

   if(ionode) then
      write(ymlout,'(/,a)') 'scatter parallelization:'

      write(ymlout,'(/,3x,a,i23)') 'total number of q-points:', sum(nq_col)
      write(ymlout,'(3x,a,i20)') 'total number of (k,q) pairs:', sum(INT8(collect_pools))

      write(ymlout,'(/,3x,a)') 'MPI task index:'

      do i = 1, npool

         write(ymlout, '(/,6x,a,a)') trim(int_to_char(i)), ':'

         write(ymlout,'(9x,a,i16)') 'number of q-points:', nq_col(i)
         write(ymlout,'(9x,a,i16)') 'number of (k,q) pairs:', collect_pools(i)

      enddo
   endif

end subroutine output_scatter_info_yaml

subroutine output_tensor_yaml(tensor_name, indent, tensor, cond_tensor_is_full)
   implicit none
   integer, intent(in)  :: indent
   character(len=*), intent(in) :: tensor_name
   real(dp), intent(in) :: tensor(:)
   logical, intent(in), optional  :: cond_tensor_is_full
   !
   logical :: cond_is_full
   character(len=64) :: fmtstr, fmtstr2, fmtstr3, fmtstr_key, fmtstr_key2
   character(len=6), external :: int_to_char

   cond_is_full = .false.

   if (present(cond_tensor_is_full)) then
      if( cond_tensor_is_full ) cond_is_full = .true.
   endif

   ! yaml keys formats
   write(fmtstr_key, '(5a)') "(/,", trim(int_to_char(indent)), "x,a,","':'",")"
   write(fmtstr_key2, '(5a)') "(/,", trim(int_to_char(indent+3)), "x,a,","':'",")"

   ! yaml number format
   write(fmtstr, '(3a)') "(", trim(int_to_char(indent+6)), "x,a,1x,E16.8)"

   ! yaml list format
   write(fmtstr3, *) "(", &
                        ! indent:
                        trim(int_to_char(indent+6)), "x,", &
                        ! list of 3 float
                        "'-'", ",1x,", "'['", ",2x,3(E16.8,", "','", ",2x),", "']'",&
                        ! end of Fortran format
                        ")"


   write(ymlout,fmtstr_key) trim(tensor_name)
   write(ymlout,fmtstr_key2) 'components'

   write(ymlout,fmtstr) 'xx:', tensor(1)
   write(ymlout,fmtstr) 'xy:', tensor(2)
   write(ymlout,fmtstr) 'yy:', tensor(3)
   write(ymlout,fmtstr) 'xz:', tensor(4)
   write(ymlout,fmtstr) 'yz:', tensor(5)
   write(ymlout,fmtstr) 'zz:', tensor(6)

   if( cond_is_full ) then
      write(ymlout,fmtstr) 'yx:', tensor(7)
      write(ymlout,fmtstr) 'zx:', tensor(8)
      write(ymlout,fmtstr) 'zy:', tensor(9)
   else
      write(ymlout,fmtstr) 'yx:', tensor(2)
      write(ymlout,fmtstr) 'zx:', tensor(4)
      write(ymlout,fmtstr) 'zy:', tensor(5)
   endif

   ! Output the tensor
   write(ymlout,'(a)') ''
   write(ymlout,fmtstr_key2) 'tensor'

   if( cond_is_full ) then
      write(ymlout, fmtstr3) tensor(1), tensor(2), tensor(4)
      write(ymlout, fmtstr3) tensor(7), tensor(3), tensor(5)
      write(ymlout, fmtstr3) tensor(8), tensor(9), tensor(6)
   else
      write(ymlout, fmtstr3) tensor(1), tensor(2), tensor(4)
      write(ymlout, fmtstr3) tensor(2), tensor(3), tensor(5)
      write(ymlout, fmtstr3) tensor(4), tensor(5), tensor(6)
   endif

end subroutine output_tensor_yaml

!convert Fortran logical variable (.true./.false.) to Python (True/False)
pure function python_bool(fortran_bool)
   implicit none
   logical, intent(in)  :: fortran_bool
   character(len=5)            :: python_bool

   if (fortran_bool) then
      python_bool = 'True'
   else
      python_bool = 'False'
   end if

end function python_bool

! Adapted from UtilXlib/clocks_handler.f90
subroutine print_clock_yaml( label )
   use util_param, only: stdout
   use mytime,     only: nclock, clock_label
   implicit none

   character(len=*) :: label
   character(len=12) :: label_
   integer          :: n

   write(ymlout, '(/,a)') 'timings:'

   write(ymlout, '(/,3x,a)') 'time units: s'

   if ( label == ' ' ) then
      !
      !
      do n = 1, nclock
         !
         call print_this_clock_yaml( n )
         !
      enddo
      !
   else
      !
      ! ... prevent trouble if label is longer than 12 characters
      !
      label_ = trim ( label )
      !
      do n = 1, nclock
         !
         if ( clock_label(n) == label_ ) then
            !
            call print_this_clock_yaml( n )
            !
            exit
            !
         endif
         !
      enddo
      !
   endif
   !
   return
   !
end subroutine print_clock_yaml


subroutine print_this_clock_yaml( n )
   use util_param, only : DP, stdout
   use mytime,     only : clock_label, cputime, walltime, mpi_per_thread, &
                          notrunning, called, t0cpu, t0wall, f_wall, f_tcpu
   !
   implicit none
   !
   integer  :: n
   real(dp) :: elapsed_cpu_time, elapsed_wall_time, nsec, msec
   integer  :: nday, nhour, nmin, nmax, mday, mhour, mmin
   !
   !
   if ( t0cpu(n) == notrunning ) then
      !
      ! ... clock stopped, print the stored value for the cpu time
      !
      elapsed_cpu_time = cputime(n)
      elapsed_wall_time= walltime(n)
      !
   else
      !
      ! ... clock not stopped, print the current value of the cpu time
      !
      elapsed_cpu_time   = cputime(n) + f_tcpu() - t0cpu(n)
      elapsed_wall_time  = walltime(n)+ f_wall() - t0wall(n)
      called(n)  = called(n) + 1
      !
   endif
   !
!#define PRINT_AVG_CPU_TIME_PER_THREAD
#if defined(PRINT_AVG_CPU_TIME_PER_THREAD)
   ! rescale the elapsed cpu time on a per-thread basis
   elapsed_cpu_time   = elapsed_cpu_time * mpi_per_thread
#endif
   !
   nmax = called(n)
   !
   ! ... In the parallel case there are several possible approaches
   ! ... The safest one is to leave each clock independent from the others
   ! ... Another possibility is to print the maximum across all processors
   ! ... This is done by uncommenting the following lines
   !
   ! call mp_max( elapsed_cpu_time, intra_image_comm )
   ! call mp_max( elapsed_wall_time, intra_image_comm )
   ! call mp_max( nmax, intra_image_comm )
   !
   ! ... In the last line we assume that the maximum cpu time
   ! ... is associated to the maximum number of calls
   ! ... NOTA BENE: by uncommenting the above lines you may run into
   ! ... serious trouble if clocks are not started on all nodes
   !

   write(ymlout, '(/,3x,a,":")') trim(clock_label(n))

   ! for clocks that have never been called
   if ( nmax == 0 ) then
      elapsed_cpu_time = 0.0_dp
      elapsed_wall_time = 0.0_dp
   endif

   write(ymlout, '(6x,a,f18.8)') 'cpu: ', elapsed_cpu_time
   write(ymlout, '(6x,a,f18.8)') 'wall:', elapsed_wall_time
   write(ymlout, '(6x,a,i10)') 'calls:', nmax

   return

end subroutine print_this_clock_yaml

end module yaml_utils
