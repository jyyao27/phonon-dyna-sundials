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

module yaml_utils_perturbo
   use yaml_utils, only: ymlout, python_bool

   public :: output_param_to_yaml_perturbo

contains

subroutine output_param_to_yaml_perturbo()
   use perturbo_autogen_output_yaml, only: auto_output_beforeconv_to_yaml, &
       auto_output_afterconv_to_yaml
   implicit none

   write(ymlout, '(/,a)') 'input parameters:'
   write(ymlout, '(3x,a)') 'before conversion (before reading from epr file):'
   call auto_output_beforeconv_to_yaml()
   write(ymlout, '(3x,a)') 'after conversion:'
   call auto_output_afterconv_to_yaml()

end subroutine output_param_to_yaml_perturbo

subroutine output_basic_data_to_yaml_perturbo()
   use pert_data, only: alat, at, bg, nat, tau, volume, nsym, symop, kc_dim, &
                        num_wann, polar_alpha, epsil, qc_dim, mass, zstar, &
                        system_2d, wannier_center, wannier_center_cryst
   use pert_param, only: spinor
   implicit none
   integer :: iv, jv
   character(len=6), external :: int_to_char
 
   write(ymlout, '(/, /,a)') 'basic data:'
   
   write(ymlout, '(3x, a, 7x, f16.10, /, 3x, a, 7x, a)') 'alat:', alat, &
                 'alat units:', 'bohr'

   write(ymlout, '(3x,a)') '# a vector is a row'
   write(ymlout, '(3x,a)') 'lattice vectors:'

   do iv = 1, 3
      write(ymlout,  '(9x,"-",1x, "[", 3(f10.5,",",2x), "]" )') at(:, iv)
   enddo
   
   write(ymlout, '(3x, a, 7x, a)') 'lattice vectors units:', 'alat'

   write(ymlout, '(3x, a)') 'reciprocal lattice vectors:'
   
   write(ymlout, '(3x,a)') '# a vector is a row'
   do iv = 1, 3
      write(ymlout, '(9x,"-",1x, "[", 3(f10.5,",",2x), "]" )') bg(:, iv)
   enddo

   write(ymlout, '(3x, a, 7x, a)') 'reciprocal lattice vectors units:', '2pi/alat'

   write(ymlout, '(3x, a, 7x, i10)') 'number of atoms in unit cell:', nat

   write(ymlout, '(3x, a)') 'atomic positions:'

   write(ymlout, '(3x,a)') '# a vector is a row'
   do iv = 1, nat
      write(ymlout, '(9x,"-",1x, "[", 3(f10.5,",",2x), "]" )') tau(:, iv)
   enddo

   write(ymlout, '(3x, a, 7x, a)') 'atomic positions units:', 'alat' 

   write(ymlout, '(3x, a, 7x, f16.10, /, 3x, a, 7x, a)') 'volume:', volume, &
                  'volume units:', 'bohr^3'
   
   write(ymlout, '(3x, a, 7x, i10)') 'number of symmetry operations:', nsym
   
   write(ymlout, '(3x, a)') 'symop:'

   do jv = 1, nsym
      write(ymlout, '(6x, a, a)') trim(int_to_char(jv)), ':'
      do iv = 1, 3
         write(ymlout, '(9x,"-",1x, "[", 3(i5,",",2x), "]" )') &
                        symop(:, iv, jv)
      enddo
   enddo

   write(ymlout, '(3x, a, /, 9x,"-",1x, "[", 3(i7,",",2x), "]")') &
                  'kc dimensions:', kc_dim(:)
   
   write(ymlout, '(3x, a, 7x, a)') 'spinor:', python_bool(spinor)
   
   write(ymlout, '(3x, a, 7x, f16.10)') 'polar_alpha:', polar_alpha
   
   write(ymlout, '(3x, a)') 'epsil:'
   write(ymlout, '(3x,a)') '# a vector is a row'
   do iv = 1, 3
      write(ymlout, '(9x,"-",1x, "[", 3(f16.8,",",2x), "]" )') epsil(:, iv)
   enddo

   write(ymlout, '(3x, a, /, 9x,"-",1x, "[", 3(i7,",",2x), "]")') &
                  'qc dimensions:', qc_dim(:)

   write(ymlout, '(3x, a)') 'mass:'
   do iv = 1, nat
      write(ymlout, '(6x, "-", 1x, f20.8)') mass(iv)
   enddo
   write(ymlout, '(3x, a, 7x, a)') 'mass units:', 'atomic'
   
   write(ymlout, '(3x, a)') 'zstar:'
   do jv = 1, nat
      write(ymlout, '(6x, a, a)') trim(int_to_char(jv)), ':'
      do iv = 1,3
      write(ymlout, '(9x,"-",1x, "[", 3(f10.5,",",2x), "]" )') &
                     zstar(iv, :, jv)
      enddo
   enddo
   
   write(ymlout, '(3x, a, 7x, a)') 'system_2d:', python_bool(system_2d)
   
   write(ymlout, '(3x, a, 7x, i10)') 'number of Wannier functions:', num_wann
   
   write(ymlout, '(3x, a)') 'wannier_center:'
   write(ymlout, '(3x,a)') '# a vector is a row'
   do iv = 1, num_wann
      write(ymlout, '(9x,"-",1x, "[", 3(f10.5,",",2x), "]" )') &
                    wannier_center(:, iv)
   enddo
   write(ymlout, '(3x, a)') 'wannier_center units: cartesian'

   write(ymlout, '(3x, a)') 'wannier_center_cryst:'
   write(ymlout, '(3x,a)') '# a vector is a row'
   do iv = 1, num_wann
      write(ymlout, '(9x,"-",1x, "[", 3(f10.5,",",2x), "]" )') &
                     wannier_center_cryst(:, iv)
   enddo
   write(ymlout, '(3x, a)') 'wannier_center_cryst units: crystal'
end subroutine output_basic_data_to_yaml_perturbo

end module yaml_utils_perturbo
