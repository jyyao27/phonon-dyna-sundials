program badshape
   use iso_c_binding
   implicit none

   integer :: test_array(5, 2)

   integer :: dims(2)
   integer(kind=c_size_t) :: dims_csizet(2)
   integer(c_int16_t) :: dims_int16(2)
   integer(c_int32_t) :: dims_int32(2)
   integer(c_int64_t) :: dims_int64(2)

   integer(kind=c_size_t) :: dims_csizet2(2)

   dims = shape(test_array)

   ! Elements of varying sizes, with matching "kind=" options
   dims_csizet = shape(test_array, kind=c_size_t)
   dims_int16 = shape(test_array, kind=c_int16_t)
   dims_int32 = shape(test_array, kind=c_int32_t)
   dims_int64 = shape(test_array, kind=c_int64_t)

   ! 64-bit elements, but no "kind=" option
   dims_csizet2 = shape(test_array)

   write(*,*) "dims=", dims
   write(*,*) "dims_csizet=", dims_csizet
   write(*,*) "dims_int16=", dims_int16
   write(*,*) "dims_int32=", dims_int32
   write(*,*) "dims_int64=", dims_int64

   write(*,*) "dims_csizet2=", dims_csizet2

end program badshape

