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
!   obtain inter-atomic force constants from TDEP
!
! Maintenance:
!===============================================================================

subroutine tdep_lattice_ifc(ph, pht, qdim, nat, at, cryst_tau)
   use kinds, only: dp
   use qe_mpi_mod, only: stdout
   use force_constant, only: lattice_ifc, set_ws_cell_ph
   use force_constant_thirdorder, only: lattice_ifct, set_tripcell_ph
   !
#if defined(__TDEP)
   use lattice_tdep, only: uc_tdep, fc_tdep, fct_tdep, init_ifc_tdep
   use lo_memtracker, only: lo_mem_helper
#endif
   !
   implicit none
   type(lattice_ifc), intent(inout) :: ph
   type(lattice_ifct), intent(inout) :: pht
   integer,  intent(in) :: qdim(3), nat
   ! at(:,i), lattice vector a_i (scaled by alat)
   ! cryst_tau: atomic position in crystal coordinate
   real(dp), intent(in) :: at(3, 3), cryst_tau(3, nat)

#if defined(__TDEP)
   type(lo_mem_helper) :: mem
   logical :: has_fct

   !setup wigner seitz cell
   call set_ws_cell_ph(ph, qdim, nat, at, cryst_tau)

   inquire(file='infile.forceconstant_thirdorder', exist=has_fct)

   if ( has_fct ) call set_tripcell_ph(pht, qdim, nat, at, cryst_tau)

   write(stdout,'(/5x,a)') "Interface to TDEP: "
   write(stdout,'( 5x,a)') "----------------------------------------------"
   write(stdout,'( 5x,a)') "read infile.ucposcar and infile.forceconstant."
   !read TDEP unit cell
   call uc_tdep%readfromfile('infile.ucposcar')
   !read TDEP forceconstant
   call mem%init()
   call fc_tdep%readfromfile(uc_tdep, 'infile.forceconstant', mem, 1)
   !read TDEP third order forceconstant
   if ( has_fct ) then
      write(stdout,'( 5x,a)') "read infile.forceconstant_thirdorder."
      call fct_tdep%readfromfile(uc_tdep, 'infile.forceconstant_thirdorder')
   endif

   !init ph with TDEP data
   if ( has_fct ) then
      call init_ifc_tdep(uc_tdep, fc_tdep, fct_tdep, ph, pht, at, cryst_tau)
   else
      call init_ifc_tdep(uc_tdep, fc_tdep, ph=ph, at=at, cryst_tau=cryst_tau)
   endif
   !
#else
   call errore('tdep_lattice_ifc',"Cannot load TDEP data. Recompile Perturbo with -D__TDEP.",1)
#endif

end subroutine tdep_lattice_ifc
