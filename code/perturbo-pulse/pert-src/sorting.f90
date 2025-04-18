module sorting
   implicit none

   ! Interface to expose C Standard Library qsort() function
   interface
      ! Invoke the C Standard Library qsort() function to sort a 1D array.
      !
      ! array is a c_ptr to the start of the array in memory
      ! count is the total number of elements in the array
      ! elem_size is the size of each array-element in bytes, using c_sizeof()
      ! f_cmp is a c_funptr to the user-provided comparison function
      subroutine c_qsort(arr, elem_count, elem_size, f_cmp) bind(C,name="qsort")
         use iso_c_binding, only: c_ptr, c_size_t, c_funptr
         implicit none
         type(c_ptr), value :: arr
         integer(c_size_t), value :: elem_count
         integer(c_size_t), value :: elem_size
         type(c_funptr), value :: f_cmp
      end subroutine c_qsort
   end interface
end module
