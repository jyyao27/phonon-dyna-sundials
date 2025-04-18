!
!!set up enviroment, adapted from QE/Modules/environment_start
!-------------------------------------------------------------
! Copyright (C) 2002-2011 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

#if defined(_OPENACC)
subroutine print_acc_info()
   use openacc
   use sys_utils, only: c_to_f_string, bytes_to_gib
#if defined(SHOW_NVML_INFO)
   use nvgpu, only: gpu_uuid, gpu_meminfo
   use iso_c_binding, only: c_char
#endif
   implicit none

   integer :: num_gpus
   integer, value :: orig_dev_num, dev_num
   integer(acc_device_kind), value :: dev_type

   character(len=80) :: str_name, str_driver
   ! integer(acc_device_kind), value :: device_mem
   integer(c_size_t), value :: device_mem, device_free_mem

#if defined(SHOW_NVML_INFO)
   ! If we have an NVIDIA GPU:
   character(kind=c_char, len=96) :: str_uuid
   integer(c_int) :: rc
   integer(c_long_long) :: free_mem, reserved_mem, total_mem, used_mem
#endif

   dev_type = acc_get_device_type()
   num_gpus = acc_get_num_devices(dev_type)

   ! Get the device number so we can restore it later
   orig_dev_num = acc_get_device_num(dev_type)

   ! write (*,'(/5X,I0," ACC devices of type ",I0)') num_gpus, dev_type
   do dev_num = 0, num_gpus - 1
      call acc_set_device_num(dev_num, dev_type)
      call acc_get_property_string(dev_num, dev_type, acc_property_name, str_name)
      call acc_get_property_string(dev_num, dev_type, acc_property_driver, str_driver)
      device_mem = acc_get_property(dev_num, dev_type, acc_property_memory)
      device_free_mem = acc_get_property(dev_num, dev_type, acc_property_free_memory)

      str_name = c_to_f_string(str_name)
      str_driver = c_to_f_string(str_driver)
      write (*,'(5X,"ACC device ",I0,":  ",A," (driver ",A,") - ",F5.2," GiB memory, ",F5.2," GiB free")') &
            dev_num, trim(str_name), trim(str_driver), &
            bytes_to_gib(device_mem), bytes_to_gib(device_free_mem)

#if defined(SHOW_NVML_INFO)
      ! If we have an NVIDIA GPU:
      rc = gpu_uuid(dev_num, str_uuid, 96)
      if (rc /= 0) str_uuid = "unknown"
         write (*,'(5X," * NVML device UUID: ",A)') c_to_f_string(str_uuid)
 
         rc = gpu_meminfo(dev_num, free_mem, reserved_mem, total_mem, used_mem)
         if (rc == 0) then
            write (*,'(5X," * NVML device memory: ",F5.2," GiB free / ",F5.2," GiB total (",F5.2," GiB reserved, ",F5.2," GiB used)")') &
                  bytes_to_gib(free_mem), bytes_to_gib(total_mem), bytes_to_gib(reserved_mem), bytes_to_gib(used_mem)
         else
            write (*,'(5X," * NVML device emory: unknown")')
         endif
#endif
   enddo
   
   if (num_gpus > 1) then
      write (*,'(5X,A,I2)') 'Using OpenACC device: ', orig_dev_num
   end if

   ! Switch back to the original device number
   call acc_set_device_num(orig_dev_num, dev_type)
   
   write (*,*)
end subroutine print_acc_info
#endif


subroutine environment_setup(code, lclock)
   use io_files,  only: crash_file, nd_nmbr
   use mp_images, only: me_image, root_image, my_image_id
   use mp_pools,  only: npool
   use mp_world,  only: nproc
   USE io_global, ONLY: stdout, meta_ionode
   use fox_init_module, ONLY: fox_init
   !
#if defined(__HDF5)
   use qeh5_base_module, only: initialize_hdf5
#else
   use hdf5_utils, only: hdf_init
#endif

   implicit none
   CHARACTER(LEN=*), INTENT(IN) :: code
   LOGICAL, intent(in) :: lclock
  
   LOGICAL :: exst, debug = .false.
   INTEGER :: ios, crashunit
   CHARACTER(LEN=80) :: uname
   CHARACTER(LEN=9)  :: cdate, ctime
   
   INTEGER, EXTERNAL :: find_free_unit
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   
   ! ... The Intel compiler allocates a lot of stack space
   ! ... Stack limit is often small, thus causing SIGSEGV and crash
   ! ... One may use "ulimit -s unlimited" but it doesn't always work
   ! ... The following call does the same and always works
   ! 
#if defined(__INTEL_COMPILER)
   CALL remove_stack_limit()
#endif

   ! ... use ".FALSE." to disable all clocks except the total cpu time clock
   ! ... use ".TRUE."  to enable clocks
   CALL init_clocks( lclock )
   CALL start_clock( trim(code) )

   IF( meta_ionode ) THEN
      ! ...  search for file CRASH and delete it
      INQUIRE( FILE=TRIM(crash_file), EXIST=exst )
      IF( exst ) THEN
         crashunit = find_free_unit()
         OPEN( UNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD',IOSTAT=ios )
         IF (ios==0) THEN
            CLOSE( UNIT=crashunit, STATUS='DELETE', IOSTAT=ios )
         ELSE
            WRITE(stdout,'(5x,"Remark: CRASH file could not be deleted")')
         END IF
      END IF

   ELSE
      ! ... one processor per image (other than meta_ionode)
      ! ... or, for debugging purposes, all processors,
      ! ... open their own standard output file
!#define DEBUG
#if defined(DEBUG)
      debug = .true.
#endif
      IF (me_image == root_image .OR. debug ) THEN
         uname = 'out.' // trim(int_to_char( my_image_id )) // '_' // &
              trim(int_to_char( me_image))
         OPEN ( unit = stdout, file = TRIM(uname),status='unknown')
      ELSE
#if defined(_WIN32)
         OPEN ( unit = stdout, file='NUL:', status='unknown' )
#else
         OPEN ( unit = stdout, file='/dev/null', status='unknown' )
#endif
      END IF

   END IF

   !print some message
   call date_and_tim(cdate, ctime)
   WRITE( stdout, '(/3X,"Program ",A," starts on ",A9," at ",A9)' ) &
      TRIM(code), cdate, ctime

   ! ... for compatibility with PWSCF
#if defined(__MPI)
   nd_nmbr = TRIM ( int_to_char( me_image+1 ))
   CALL parallel_info(code, stdout)
#else
   nd_nmbr = ' '
   CALL serial_info(stdout)
#endif
   
   !open HDF5 interface
   !if __HDF5 is used in QE/make.inc, 
   !  then we initialized the HDF5 interface with qeh5_module subroutine
   !if __HDF5 is not used in QE, then we use our own HDF5 wrapper.
#if defined(__HDF5)
   call initialize_hdf5()
#else
   call hdf_init()
#endif

   CALL fox_init()
   !check if npool == nproc
   if (nproc /= npool) &
      CALL errore('qe2pert','npools must match the number of MPI processes',1)
end subroutine environment_setup


subroutine environment_stop( code )
   USE io_global, ONLY: stdout, meta_ionode
   use qe_mpi_mod, only: ionode
   !
#if defined(__HDF5)
   use qeh5_base_module, only: finalize_hdf5
#else
   use hdf5_utils, only: hdf_finalize
#endif

   implicit none
   CHARACTER(LEN=*), INTENT(IN) :: code
   !
   CHARACTER(LEN=9)  :: cdate, ctime
   CHARACTER(LEN=80) :: time_str

   CALL stop_clock(  TRIM(code) )
   !CALL print_clock( TRIM(code) )
   !print out all the clock
   CALL print_clock(' ')

   !close HDF5 interface
#if defined(__HDF5)
   call finalize_hdf5()
#else
   call hdf_finalize()
#endif

   call date_and_tim( cdate, ctime )
   time_str = 'Program was terminated on:  ' // ctime // ' ' // cdate

   if(meta_ionode) write(stdout,'(/3X, A60)') time_str
end subroutine environment_stop


  !==-----------------------------------------------------------------------==!
  SUBROUTINE parallel_info(code, stdout)
    use mp_world,  only: nproc, nnode
    implicit none
    CHARACTER(LEN=*), INTENT(IN) :: code
    integer, intent(in) :: stdout
#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_max_threads
#endif

    CHARACTER(LEN=40) :: apis

    apis = "MPI"
#if defined(_OPENMP)
    apis = trim(apis) // " & OpenMP"
#endif
#if defined(_OPENACC)
    apis = trim(apis) // " & OpenACC"
#endif

#if defined(_OPENMP)
    !
    WRITE( stdout, '(/5X,"Parallel version (' // trim(apis) // '), running on ", &
         &I7," processor cores")' ) nproc * omp_get_max_threads()
    !
    WRITE( stdout, '(5X,"Number of MPI processes:           ",I7)' ) nproc
    !
    WRITE( stdout, '(5X,"Threads/MPI process:               ",I7/)' ) &
         omp_get_max_threads()
#else
    WRITE( stdout, '(/5X,"Parallel version (' // trim(apis) // '), running on ",&
         &I5," processors")' ) nproc 
#endif
    !
#if !defined(__GFORTRAN__) ||  ((__GNUC__>4) || ((__GNUC__==4) && (__GNUC_MINOR__>=8)))
    WRITE( stdout, '(5X,"MPI processes distributed on ",&
         &I5," nodes",/)' ) nnode
#endif

#if defined(_OPENACC)
    call print_acc_info()
#endif
    !
  END SUBROUTINE parallel_info

  !==-----------------------------------------------------------------------==!
  SUBROUTINE serial_info ( stdout )
    implicit none
    integer, intent(in) :: stdout
    !
#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_max_threads
#endif

    CHARACTER(LEN=40) :: apis

    ! Set up apis to report OpenMP and/or OpenACC as appropriate
    apis = ""
#if defined(_OPENMP)
#if defined(_OPENACC)
    apis = "OpenMP & OpenACC"
#else
    apis = "OpenMP"
#endif
#elif defined(_OPENACC)
    apis = "OpenACC"
#endif

    !
#if defined(_OPENMP)
    WRITE( stdout, '(/5X,"Serial multi-threaded version (' // trim(apis) // '), running on ",&
         &I4," processor cores")' ) omp_get_max_threads()
    !
#elif defined(_OPENACC)
    WRITE( stdout, '(/5X,"Serial version (OpenACC)")' )
#else
    WRITE( stdout, '(/5X,"Serial version")' )
#endif

#if defined(_OPENACC)
    call print_acc_info()
#endif
    !
  END SUBROUTINE serial_info
