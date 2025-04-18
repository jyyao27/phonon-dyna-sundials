module nvgpu
   implicit none

   ! Interface to expose NVML functions that report GPU details
   interface
      ! Retrieve more detailed memory information from the NVIDIA GPU.  This
      ! includes the following:
      !  * "free" - Unallocated device memory (in bytes).
      !  * "reserved" - Device memory (in bytes) reserved for system use
      !    (driver or firmware).
      !  * "total" - Total physical device memory (in bytes).
      !  * "used" - Allocated device memory (in bytes).
      ! The OpenACC memory information reported is comprised from these
      ! details, so just using the OpenACC data is perfectly adequate.
      integer(c_int) function gpu_meminfo(gpu_index, free, reserved, total, used) bind(C,name="get_gpu_meminfo")
         use iso_c_binding, only: c_char, c_int, c_long_long
         implicit none
         integer(c_int), value :: gpu_index
         integer(c_long_long), intent(out) :: free, reserved, total, used
      end function gpu_meminfo

      ! Retrieve the GPU device's globally unique UUID.
      integer(c_int) function gpu_uuid(gpu_index, uuid_str, uuid_str_len) bind(C,name="get_gpu_uuid")
         use iso_c_binding, only: c_char, c_int
         implicit none
         integer(c_int), value :: gpu_index, uuid_str_len
         character(kind=c_char,len=*), intent(inout) :: uuid_str
      end function gpu_uuid
   end interface
end module

