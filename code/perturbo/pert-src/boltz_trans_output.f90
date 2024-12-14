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

module boltz_trans_output
   use pert_const, only: dp, ryd2ev, bohr2ang, kelvin2ev, ryd2mev, unitcharge,&
                         timeunit, tesla2ryd
   use pert_data,  only: tpiba, bg, system_2d, alat, at, volume
   use pert_utils, only: find_free_unit, mfermi_deriv
   use pert_param, only: prefix, boltz_kdim
   use boltz_utils,only: num2kpt
   use boltz_grid, only: grid
   implicit none
   private
   public :: output_kgrid,   output_dos,   output_tdf, output_trans_coef
   public :: output_trans_coef_yaml
   public :: output_density, output_rates, output_mobility, output_ftemper
contains

!output kgrid
subroutine output_kgrid(kg)
   implicit none
   type(grid), intent(in) :: kg
   ! local variables
   integer :: ik, uout
   real(dp) :: xk(3), xkc(3)

   uout = find_free_unit()
   open(uout, file=trim(prefix)//'_fullgrid.kpt', status='unknown', form='formatted')
   write(uout, '(1x, i8, 2x, a)')  kg%nk, &
      'Col.1:index;  (2, 3, 4): kpt in crystal; (5, 6, 7): kpt in cartestian (2pi/a).'
   do ik = 1, kg%nk
      xk = num2kpt( kg%kpt(ik), boltz_kdim)
      xkc(:) = xk(:)
      call cryst_to_cart(1, xkc, bg, 1)
      write(uout,'(5x,i8,1x,3(f10.6,2x),2x,3(f10.6,2x))') ik, xk(1:3), xkc(1:3)
   enddo
   close(uout)
end subroutine output_kgrid

subroutine output_dos(energy, dos)
   use yaml_utils, only: ymlout
   implicit none
   real(dp), intent(in) :: energy(:), dos(:)
   ! local variables
   integer :: uout, i
   integer, external :: find_free_unit

   uout = find_free_unit()
   open(unit=uout, file=trim(prefix)//'.dos',status='unknown',form='formatted')
   write(uout,'(2x,"#   E (eV)     DOS (#.states/eV/u.c.)")')
   do i = 1, size(energy)
      write(uout,'(2x, E14.6, 3x, E20.10)') energy(i)*ryd2ev, dos(i)/ryd2ev
   enddo
   close(uout)

   ! Output density to the YAML file
   write(ymlout, '(/,a)') 'DOS:'

   write(ymlout, '(/,3x,a)') 'energy units: eV'

   write(ymlout, '(/,3x,a)') 'DOS units: num. of states/eV/u.c.'

   write(ymlout, '(/,3x,a)') 'energy:'

   do i = 1, size(energy)
      write(ymlout,'(6x, "-", 1x, E14.6)') energy(i)*ryd2ev
   enddo

   write(ymlout, '(/,3x,a)') 'DOS:'

   do i = 1, size(energy)
      write(ymlout,'(6x, "-", 1x, E20.10)') dos(i)/ryd2ev
   enddo

end subroutine output_dos

subroutine output_rates(kg, tempe, efermi, rates, append)
   implicit none
   type(grid), intent(in) :: kg
   logical, intent(in) :: append ! write in append mode
   real(dp), intent(in) :: tempe, efermi, rates(:,:) ! rates(kg%numb, kg%nk)
   ! local variables
   integer :: uout, ik, ib
   logical :: created

   uout = find_free_unit()
   inquire(file=trim(prefix)//'.rates', exist=created)
   if(append .and. created) then
      open(uout, file=trim(prefix)//'.rates',status='old',action='write',position='append')
   else
      open(uout, file=trim(prefix)//'.rates',status='unknown',form='formatted')
      write(uout,'(1x,"#    Scattering rates computed in transport mode        #")')
      write(uout,'(1x,"#    (\Gamma in meV and N.B. \Gamma = 2 ImSigma)        #")')
   endif
   write(uout, '(a)')
   write(uout,'(1x, a, f9.4, a, f10.6/)') &  
      '#  Temperature: ', tempe*ryd2ev/kelvin2ev, '  Chemical Potential: ', efermi*ryd2ev
   do ik = 1, kg%nk
      do ib = 1, kg%numb
         write(uout,'(1x, i10, 3x, i4.4, 3x, f12.6, 3x, E22.10)') &
            ik, ib, kg%enk(ib,ik)*ryd2ev, rates(ib,ik)*ryd2mev
      enddo
   enddo
   close(uout)
end subroutine output_rates

subroutine output_tdf(energy, tdf, tempe, efermi, calc_mag, append)
   implicit none
   real(dp), intent(in) :: energy(:), tdf(:,:), tempe, efermi   !tdf(6,:)
   logical, intent(in), optional :: append
   logical, intent(in) :: calc_mag
   ! local variables
   logical :: created
   integer :: uout, i, j
   
   uout = find_free_unit()
   inquire(file=trim(prefix)//'.tdf', exist=created)
   !
   if( present(append) ) created = append .and. created
   if(created) then
      open(uout, file=trim(prefix)//'.tdf',status='old',action='write',position='append')
   else
      open(uout, file=trim(prefix)//'.tdf',status='replace',form='formatted')
      if(.not. calc_mag) then
         write(uout,'(1x,"#  E(eV)    (-df/dE) (a.u.)    TDF(E)_(xx xy yy xz yz zz) (a.u.)   #")')
      else
         write(uout,'(1x,"#  E(eV)    (-df/dE) (a.u.)    TDF(E)_(xx xy yy xz yz zz yx zx zy) (a.u.)   #")')
   endif;endif
   write(uout,'(/1x, a, f9.4, a, f10.6/)') &
      '# Temperature: ', tempe*ryd2ev/kelvin2ev, '  Chemical Potential: ', efermi*ryd2ev
   do i = 1, size(energy)
      if(.not. calc_mag) then
         write(uout,'(f12.6, 3x, ES23.16, 2x, 6(E14.6,1x))') energy(i)*ryd2ev, &
            mfermi_deriv(tempe, (energy(i)-efermi)),  (tdf(j,i), j=1,6)
      else
         write(uout,'(f12.6, 3x, ES23.16, 2x, 9(E14.6,1x))') energy(i)*ryd2ev, &
            mfermi_deriv(tempe, (energy(i)-efermi)),  (tdf(j,i), j=1,9)
      endif
   enddo
   close(uout)
end subroutine output_tdf

subroutine output_density(temper, efermi, density)
   use yaml_utils, only: ymlout
   implicit none
   real(dp), intent(in) :: temper(:), efermi(:), density(:)
   character(len=7) :: label
   integer :: uout, it
   real(dp) :: dens_unit
   character(len=6), external :: int_to_char

   label = merge("(cm^-2)", "(cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   ! output carrer concentration
   uout = find_free_unit()
   open(uout, file=trim(prefix)//'.doping', status='unknown', form='formatted')
   write(uout,'(1x,a,a)') "#", repeat('=',79)
   write(uout,'(1x,a)') "# Temperature(K), Chemical Potential(eV), Carrier concentration" // label
   do it = 1, size(temper)
      write(uout,'(4x, f7.2, 7x, f15.10, 9x, E15.7)')  &
         temper(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, density(it)*dens_unit
   enddo
   close(uout)

   ! Output density to the YAML file
   write(ymlout, '(/,a)') 'carrier density:'

   write(ymlout, '(/,3x,a)') 'temperature units: K'

   write(ymlout, '(/,3x,a)') 'chemical potential units: eV'

   write(ymlout, '(/,3x,a,1x,a)') 'concentration units:', merge("cm-2", "cm-3", system_2d)

   write(ymlout, '(/,3x,a)') '# Cofiguration means one line in the temper file:'
   write(ymlout, '(3x,a)') '# temperature efermi concentration'
   write(ymlout, '(3x,a,i4)') 'number of configurations:', size(temper)

   write(ymlout, '(/,3x,a)') 'configuration index:'

   do it = 1, size(temper)
      write(ymlout, '(/,6x,a,a)') trim(int_to_char(it)), ':'

      write(ymlout,'(9x, a, f7.2)') 'temperature:', temper(it)*ryd2ev/kelvin2ev
      write(ymlout,'(9x, a, f15.10)') 'chemical potential:', efermi(it)*ryd2ev
      write(ymlout,'(9x, a, E15.7)') 'concentration:', density(it)*dens_unit

   enddo

end subroutine output_density


subroutine output_ftemper(temper, efermi, density, fname)
   implicit none
   real(dp), intent(in) :: temper(:), efermi(:), density(:)
   character(len=*), intent(in) :: fname
   integer :: uout, it
   real(dp) :: dens_unit
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)
   ! output carrer concentration
   uout = find_free_unit()
   open(uout, file=trim(fname), status='unknown', form='formatted')
   write(uout,'(i5, 5x)')  size(temper)
   do it = 1, size(temper)
      write(uout,'(1x, f7.2, 7x, f15.10, 9x, E15.7)')  &
         temper(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, density(it)*dens_unit
   enddo
   close(uout)
end subroutine


subroutine output_mobility(tempe, efermi, dens, cond, niter_in)
   implicit none
   ! cond(6 or 9, max_step, ntempe)
   real(dp), intent(in) :: tempe(:), efermi(:), dens(:), cond(:,:,:) !cond(6 or 9,:,:)
   integer, intent(in), optional :: niter_in(:)
   integer :: uout, it, ia, i, max_step, ntmp, last, ncomp_print, ncomp
   integer, allocatable :: niter(:)
   real(dp) :: cond_unit, mobi_unit, dens_unit
   character(len=220) :: ctmp, label, ctmp1, ctmp2, fmt1, fmt2, fmt3, fmt4,&
                         fmt5, fmt6
   character(len=3) :: sub(9) = (/'_xx', '_xy', '_yy', '_xz', '_yz', '_zz',&
                                  '_yx', '_zx', '_zy'/)
         
   label = merge("n_c (cm^-2)", "n_c (cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   ! unit conversion of cond. check below.
   ! \sigma = q^2*nspin/(N*V) sum_{nk} v_nk*v_nk*tau_nk* f *(1-f)/kT
   !  (N.B.: the factor alat^2 should have already included in cond at this point)
   !
   ! 1. q^2 is omitted in the calculation of cond(i), 
   ! 2. all the quantities are in Rydberg atomic unit.
   !--------------------------------------------------------------------------
   !cond_unit = q^2* (a0)^-3* (Ryd/(hbar*a0^-1))^2 * (hbar/Ryd) / Ryd
   !          = e^2 * (hbar*a0)^-1 = e^2 *(Ryd*t0*a0)^-1
   !convert cond to (ohm*m)^-1  (ohm = V/A = V*s/C; hbar = Ryd*t0)
   cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)
   !for 2D, we output conductance
   if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !convert from Rydberg atomic unit to cm^2/V/s
   !mobility mu = sigma / (n*q); (n:(a0)^-3); mu = sigma * a0^3 / q 
   !mobi_unit = q* (a0)^2/(Ryd*t0)
   mobi_unit = bohr2ang*bohr2ang/(timeunit*1.0E16_dp*ryd2ev)
   
   ntmp = size(tempe)
   max_step = size(cond, 2)

   ncomp = size(cond, 1)
   ! ncomp_print designates the number of cond elements to print
   ! ncomp_print is 6 for the symmetric cond tensor (non-magnetic case)
   ! ncomp_print is 9 for a cond tensor without symmetries (magnetic case)

   ! Non-magnetic case
   if( ncomp == 6 .or. ncomp == 12 ) then
      ncomp_print = 6
   ! Magnetic case
   else
      ncomp_print = 9
   endif

   !Create string formats based on whether to print 6 or 9 elements
   write(ctmp1,*) ncomp_print
   write(ctmp2,*) merge(134,90,(ncomp_print .eq. 9))
   write(fmt1,*) '('//trim(ctmp1)//'(5x,a8,2x))'
   write(fmt2,*) '(a8, 3x,a7, 3x,a11, 1x,a'//trim(ctmp2)//')'
   write(fmt3,*) '(f8.2, 1x, f9.5, 1x, E13.5, 1x, '//trim(ctmp1)//'(1x, E14.6))'
   write(fmt4,*) '(2x,a6,2x,a'//trim(ctmp2)//')'
   write(fmt5,*) '(2x,a1,i4,3x,'//trim(ctmp1)//'(1x, E14.6))'
   write(fmt6,*) '('//trim(ctmp1)//'(6x,a5,4x))'

   !
   allocate( niter( ntmp ) )
   niter = 1
   if(present( niter_in )) then
      if(size(niter) .ne. ntmp) call errore('output_mobility', 'niter dimension mismatch', 1)
      niter(:) = niter_in(:)
   else
      do it = 1, ntmp
         do i = 1, max_step
            if(sum(abs(cond(:,i,it))) < 1.0E-16_dp) exit
            niter(it) = i
         enddo
      enddo
   endif
   !
   ! array chcek
   if(size(dens).ne.ntmp .or. size(cond,3).ne.ntmp) &
      call errore('output_mobility','array dimension mismatch',1)
   uout = find_free_unit()
   open(uout, file=trim(prefix)//'.cond', status='unknown', form='formatted')

   write(uout,'(/10x,"#==========================================================#" )')

   if(system_2d) then
      write(uout,'( 10x,"#                   Conductivity (1/Ohm)                   #" )')
      write(uout,'( 10x,"#----------------------( 2D system )-----------------------#"/)')
   else
      write(uout,'( 10x,"#                  Conductivity (1/Ohm/m)                  #" )')
      write(uout,'( 10x,"#----------------------------------------------------------#"/)')
   endif

   write(ctmp, fmt1) ('sigma'//sub(ia), ia = 1, ncomp_print)
   write(uout, fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      last = niter(it)
      write(uout,fmt3) &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (cond(ia, last, it)*cond_unit, ia=1, ncomp_print)
      !for iterative approach
      if(max_step > 1) then
         write(uout,'( 10x,"#--------------------iterative process---------------------#" )')
         write(uout,fmt4) "#iter.", ctmp
         do i = 1, niter(it)
            write(uout,fmt5) '#', i, (cond(ia,i,it)*cond_unit, ia=1, ncomp_print)
         enddo
         write(uout,'( 10x,"#----------------------------------------------------------#"/)')
      endif
   enddo

   write(uout,'(/a)') 
   write(uout,'(/10x,"#==========================================================#" )')
   write(uout,'( 10x,"#                    Mobility (cm^2/V/s)                   #" )')
   write(uout,'( 10x,"#--------------------(for semiconductor)-------------------#"/)')

   write(ctmp,fmt6) ('mu'//sub(ia), ia = 1, ncomp_print)
   write(uout,fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      last = niter(it)
      write(uout,fmt3) &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (cond(ia, last, it)/dens(it)*mobi_unit, ia=1, ncomp_print)
   enddo
   deallocate(niter)
   close(uout)
end subroutine output_mobility

subroutine output_trans_coef(tempe, efermi, dens, cond, seebeck, t_cond, is_mag)
   implicit none
   real(dp), intent(in) :: tempe(:), efermi(:), dens(:), cond(:,:), seebeck(:,:), t_cond(:,:)
   logical, optional, intent(in) :: is_mag
   !
   character(len=220) :: ctmp, label, ctmp1, ctmp2, fmt1, fmt2, fmt3, fmt4,&
                         fmt5, fmt6
   character(len=3) :: sub(9) = (/'_xx', '_xy', '_yy', '_xz', '_yz', '_zz',&
                                  '_yx', '_zx', '_zy'/)
   integer  :: uout, it, ia, ntmp, ncomp
   real(dp) :: dens_unit, seebeck_unit, cond_unit, mobi_unit, t_cond_unit
   logical  :: calc_mag

   ntmp = size(tempe)
   calc_mag = .false.
   if(present(is_mag)) calc_mag = is_mag

   ! (k_b/q) = (- k_b/e) = - k_b*(K) / (e*K) = - (kelvin2ev * eV) / (e* K) = -kelvin2ev * (V/K)
   seebeck_unit = kelvin2ev * 1.0E3_dp  ! convert to mV/K
   
   label = merge("n_c (cm^-2)", "n_c (cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   !see comments in boltz_trans_output.f90/output_mobility for more detail.
   !cond_unit = q^2* (a0)^-3* (Ryd/(hbar*a0^-1))^2 * (hbar/Ryd) / Ryd
   !          = e^2 * (hbar*a0)^-1 = e^2 *(Ryd*t0*a0)^-1
   !convert cond to (ohm*m)^-1  (ohm = V/A = V*s/C; hbar = Ryd*t0)
   cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)
   if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !thermal conductivity (1/T * Ryd^2 *(cond / e^2) )
   ! \kappa = (1/T)* Ryd^2 * (Ryd*t0*a0)^-1 = (k_b/Ryd) * Ryd^2 * (Ryd*t0*a0)^-1
   !        = k_b / (t0 * a0) = (k_b*K) / (t0 * a0 * K) = kelvin2ev * (e*V/t0) / a0 / K
   ! convert to W/meter/K, W = J/s = C*V / s
   t_cond_unit = kelvin2ev * (unitcharge/timeunit) * 1.0E10_dp / bohr2ang
   if(system_2d) t_cond_unit = (t_cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !convert from Rydberg atomic unit to cm^2/V/s
   !mobility mu = sigma / (n*q); (n:(a0)^-3); mu = sigma * a0^3 / q 
   !mobi_unit = q* (a0)^2/(Ryd*t0)
   mobi_unit = bohr2ang*bohr2ang/(timeunit*1.0E16_dp*ryd2ev)

   ncomp = size(cond, 1)
   !Create string formats based on whether to print 6 or 9 elements
   write(ctmp1,*) ncomp
   write(ctmp2,*) merge(134,90,(ncomp .eq. 9))
   write(fmt1,*) '('//trim(ctmp1)//'(5x,a8,2x))'
   write(fmt2,*) '(a8, 3x,a7, 3x,a11, 1x,a'//trim(ctmp2)//')'
   write(fmt3,*) '(f8.2, 1x, f9.5, 1x, E13.5, 1x, '//trim(ctmp1)//'(1x, E14.6))'
   write(fmt4,*) '('//trim(ctmp1)//'(7x,a4,4x))'
   write(fmt5,*) '('//trim(ctmp1)//'(5x,a8,2x))'
   write(fmt6,*) '('//trim(ctmp1)//'(6x,a5,4x))'
   
   uout = find_free_unit()
   open(uout, file=trim(prefix)//'.trans_coef', status='unknown', form='formatted')

   write(uout,'(/10x,"#==========================================================#" )')
if(system_2d) then
   write(uout,'( 10x,"#                   Conductivity (1/Ohm)                   #" )')
   write(uout,'( 10x,"#----------------------( 2D system )-----------------------#"/)')
else
   write(uout,'( 10x,"#                  Conductivity (1/Ohm/m)                  #" )')
   write(uout,'( 10x,"#----------------------------------------------------------#"/)')
endif

   write(ctmp, fmt1) ('sigma'//sub(ia), ia = 1, ncomp)
   write(uout, fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      write(uout,fmt3) &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (cond(ia,it)*cond_unit, ia=1, ncomp)
   enddo

   write(uout,'(/a)') 
   write(uout,'(/10x,"#==========================================================#" )')
   write(uout,'( 10x,"#                    Mobility (cm^2/V/s)                   #" )')
   write(uout,'( 10x,"#--------------------(for semiconductor)-------------------#"/)')

   write(ctmp,fmt6) ('mu'//sub(ia), ia = 1, ncomp)
   write(uout,fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      write(uout,fmt3) &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (cond(ia,it)/dens(it)*mobi_unit, ia=1,ncomp)
   enddo
   if(.not. calc_mag) then
      write(uout,'(/a)') 
      write(uout,'(/10x,"#==========================================================#" )')
      write(uout,'( 10x,"#                Seebeck coefficient (mV/K)                #" )')
      write(uout,'( 10x,"#----------------------------------------------------------#"/)')
      !
      write(ctmp, fmt4) ('S'//sub(ia), ia = 1, ncomp)
      write(uout, fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
      do it = 1, ntmp
         write(uout,fmt3) &
            tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
            (seebeck(ia,it)*seebeck_unit, ia=1,ncomp)
      enddo
      
      write(uout,'(/a)') 
      write(uout,'(/10x,"#==========================================================#" )')
      if(system_2d) then
         write(uout,'( 10x,"#                Thermal conductivity ( W/K )              #" )')
         write(uout,'( 10x,"#----------------------( 2D system )-----------------------#" )')
      else
         write(uout,'( 10x,"#                Thermal conductivity (W/m/K)              #" )')
      endif
      write(uout,'( 10x,"#------------------(Electronic contribution)---------------#"/)')
      !
      write(ctmp, fmt5) ('kappa'//sub(ia), ia = 1, ncomp)
      write(uout, fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
      do it = 1, ntmp
         write(uout,fmt3) &
            tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
            (t_cond(ia,it)*t_cond_unit, ia=1,ncomp)
      enddo
   endif
end subroutine output_trans_coef


subroutine output_trans_coef_yaml(tempe, efermi, dens, cond, cond_iter, niter, seebeck, t_cond, bfield, is_mag)
   use yaml_utils, only: ymlout, output_tensor_yaml
   implicit none
   integer, intent(in)  :: niter(:)
   real(dp), intent(in) :: tempe(:), efermi(:), dens(:), cond(:,:), cond_iter(:,:,:)
   real(dp), intent(in) :: seebeck(:,:), t_cond(:,:), bfield(:,:)
   logical, optional, intent(in) :: is_mag
   !
   character(len=120) :: ctmp, label
   character(len=3) :: sub(6) = (/'_xx', '_xy', '_yy', '_xz', '_yz', '_zz'/)
   logical :: calc_mag
   integer  :: uout, it, ia, ntmp, iiter, max_step
   real(dp) :: dens_unit, seebeck_unit, cond_unit, mobi_unit, t_cond_unit
   character(len=6), external :: int_to_char

   calc_mag = .false.
   if(present(is_mag)) calc_mag = is_mag

   ntmp = size(tempe)

   ! (k_b/q) = (- k_b/e) = - k_b*(K) / (e*K) = - (kelvin2ev * eV) / (e* K) = -kelvin2ev * (V/K)
   seebeck_unit = kelvin2ev * 1.0E3_dp  ! convert to mV/K
   
   label = merge("n_c (cm^-2)", "n_c (cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   !see comments in boltz_trans_output.f90/output_mobility for more detail.
   !cond_unit = q^2* (a0)^-3* (Ryd/(hbar*a0^-1))^2 * (hbar/Ryd) / Ryd
   !          = e^2 * (hbar*a0)^-1 = e^2 *(Ryd*t0*a0)^-1
   !convert cond to (ohm*m)^-1  (ohm = V/A = V*s/C; hbar = Ryd*t0)
   cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)
   if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !thermal conductivity (1/T * Ryd^2 *(cond / e^2) )
   ! \kappa = (1/T)* Ryd^2 * (Ryd*t0*a0)^-1 = (k_b/Ryd) * Ryd^2 * (Ryd*t0*a0)^-1
   !        = k_b / (t0 * a0) = (k_b*K) / (t0 * a0 * K) = kelvin2ev * (e*V/t0) / a0 / K
   ! convert to W/meter/K, W = J/s = C*V / s
   t_cond_unit = kelvin2ev * (unitcharge/timeunit) * 1.0E10_dp / bohr2ang
   if(system_2d) t_cond_unit = (t_cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !convert from Rydberg atomic unit to cm^2/V/s
   !mobility mu = sigma / (n*q); (n:(a0)^-3); mu = sigma * a0^3 / q 
   !mobi_unit = q* (a0)^2/(Ryd*t0)
   mobi_unit = bohr2ang*bohr2ang/(timeunit*1.0E16_dp*ryd2ev)
   

   ! Output to YAML starts here
   write(ymlout, '(/,a)') 'trans:'

   write(ymlout, '(/,3x,a)') 'temperature units: K'

   write(ymlout, '(/,3x,a)') 'chemical potential units: eV'

   write(ymlout, '(/,3x,a,1x,a)') 'concentration units:', merge("cm-2", "cm-3", system_2d)

   if(calc_mag) write(ymlout, '(/,3x,a)') 'magnetic field units: T'

   if(system_2d) then
      write(ymlout, '(/,3x,a)') 'conductance units: 1/Ohm'
   else
      write(ymlout, '(/,3x,a)') 'conductivity units: 1/Ohm/m'
   endif

   write(ymlout, '(/,3x,a)') 'mobility units: cm2/V/s'

   if(.not. calc_mag) then
      write(ymlout, '(/,3x,a)') 'Seebeck coefficient units: mV/K'

      write(ymlout, '(/,3x,a)') '# Electronic contribution'
      write(ymlout, '(3x,a)') '# Not well tested, use with caution'
      write(ymlout, '(3x,a)') 'thermal conductivity units: W/m/K'
   endif

   write(ymlout, '(/,3x,a)') '# Configuration means one line in the temper file:'
   write(ymlout, '(3x,a)') '# temperature efermi concentration'
   write(ymlout, '(3x,a,i4)') 'number of configurations:', ntmp


   write(ymlout, '(/,3x,a)') 'configuration index:'

   do it = 1, ntmp

      write(ymlout, '(/,6x,a,a)') trim(int_to_char(it)), ':'

      write(ymlout,'(/,9x, a, f10.5)') 'temperature:', tempe(it)*ryd2ev/kelvin2ev

      write(ymlout,'(/,9x, a, f15.10)') 'chemical potential:', efermi(it)*ryd2ev

      write(ymlout,'(/,9x, a, E15.7)') 'concentration:', dens(it)*dens_unit

      if(calc_mag) write(ymlout,'(/,9x, a, 1x, "[",3(f10.5,",",2x),"]")') &
              'magnetic field:', (bfield(ia,it)/tesla2ryd,ia=1,3)

      ! output the final conductivity
      call output_tensor_yaml('conductivity', 9, cond(:, it) * cond_unit, calc_mag)

      ! mobility
      call output_tensor_yaml('mobility', 9, cond(:, it) / dens(it) *mobi_unit,&
      calc_mag)

      if(.not. calc_mag) then
         ! Seebeck coefficient
         call output_tensor_yaml('Seebeck coefficient', 9, seebeck(:, it) * seebeck_unit)
         ! Thermal conductivity
         write(ymlout,'(/,12x,a)') '# Thermal conductivity: not well tested, use with caution'
         call output_tensor_yaml('thermal conductivity', 9, t_cond(:, it) * t_cond_unit)
      endif

      ! output conductivity iterations
      max_step = size(cond_iter, 2)
      if(max_step > 1) then
         write(ymlout,'(/,9x, a, i5)') 'number of iterations:', niter(it)
         write(ymlout, '(/,9x,a)') 'iteration:'
         do iiter = 1, niter(it)
            write(ymlout, '(/,12x,a,a)') trim(int_to_char(iiter)), ':'

            call output_tensor_yaml('conductivity', 15, cond_iter(:, iiter, it) * cond_unit, calc_mag)

         enddo
      endif

   enddo ! ntmp

end subroutine output_trans_coef_yaml

end module boltz_trans_output
