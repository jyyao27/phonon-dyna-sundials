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

module yaml_utils_qe2pert
   use yaml_utils, only: ymlout, python_bool

   public :: output_param_to_yaml_qe2pert

contains

subroutine output_param_to_yaml_qe2pert()
   use qe2pert_autogen_output_yaml, only: auto_output_beforeconv_to_yaml, &
       auto_output_afterconv_to_yaml
   implicit none

   write(ymlout, '(/,a)') 'input parameters:'
   write(ymlout, '(3x,a)') 'before conversion:'
   call auto_output_beforeconv_to_yaml()
   write(ymlout, '(3x,a)') 'after conversion:'
   call auto_output_afterconv_to_yaml()

end subroutine

end module yaml_utils_qe2pert
