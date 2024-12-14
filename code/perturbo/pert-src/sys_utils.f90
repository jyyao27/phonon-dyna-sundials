module sys_utils
   private
   public :: c_to_f_string, find_proc_file_value
contains


! Convert a C-string (with trailing NUL characters) into a Fortran string.
function c_to_f_string(cstr) result(fstr)
   use iso_c_binding, only: c_null_char
   implicit none

   character(len=*) :: cstr
   character(len=len(cstr)) :: fstr

   integer :: pos

   pos = index(cstr, c_null_char)
   if (pos > 0) then
      fstr = cstr(1:pos-1)
   else
      fstr = cstr
   endif
end function c_to_f_string


!!! Find the value of a key in a Linux process file, i.e. something under
!!! /proc/.  For example, to find the total physical memory available on
!!! the system, this can be run:
!!!     find_proc_file_value('/proc/meminfo', 'MemTotal', 'unknown')
function find_proc_file_value(filepath, key, default_value) result(val)
   use pert_utils, only: find_free_unit
   implicit none

   character(len=*) :: filepath, key, default_value
   character(len=80) :: line, val
   integer :: funit, ios, pos

   val = default_value
   funit = find_free_unit()
   open(unit = funit, file = filepath, action = 'read', iostat = ios)
   if (ios == 0) then
      ! File-open succeeded.  Scan through lines of the file, looking for
      ! the specified keyword.
      do
         read (funit, '(A)', iostat = ios) line
         if (ios /= 0) then
            exit
         endif

         if (index(line, key // ':') == 1) then
            val = line(len(key)+2:)
            exit
         endif
      enddo

      close(funit)

      ! Lines in /proc files uses spaces and tabs between the keyword and
      ! the associated value.  Trim these out.
      pos = 1
      do while (pos <= len(val) .and. (val(pos:pos) == ' ' .or. val(pos:pos) == achar(9)))
         pos = pos + 1
      end do
      val = val(pos:)
   endif

end function find_proc_file_value

end module sys_utils

