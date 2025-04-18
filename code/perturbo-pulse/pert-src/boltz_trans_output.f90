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
                         timeunit, tesla2ryd, pi
   use pert_data,  only: tpiba, bg, system_2d, alat, at, volume
   use pert_utils, only: find_free_unit, mfermi_deriv
   use pert_param, only: prefix, boltz_kdim, ph_e, drag
   use boltz_utils,only: num2kpt
   use boltz_grid, only: grid
   implicit none
   private
   public :: output_kgrid,   output_dos,   output_tdf, output_trans_coef
   public :: output_trans_coef_yaml
   public :: output_density, output_rates, output_mobility, output_ftemper
   ! syp -
   public :: output_ph_cond, convert_cond
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

subroutine output_rates(kg, tempe, efermi, rates, append, tmpfile)
   implicit none
   type(grid), intent(in) :: kg
   logical, intent(in) :: append ! write in append mode
   real(dp), intent(in) :: tempe, efermi, rates(:,:) ! rates(kg%numb, kg%nk)
   character(len=*), intent(in), optional :: tmpfile
   ! local variables
   integer :: uout, ik, ib
   logical :: created
   character(len=100) :: filename_local

   !syp - added
   filename_local = trim(prefix)//'.rates'
   if(present(tmpfile)) filename_local = tmpfile

   uout = find_free_unit()
   inquire(file=filename_local, exist=created)
   if(append .and. created) then
      open(uout, file=filename_local,status='old',action='write',position='append')
   else
      open(uout, file=filename_local,status='unknown',form='formatted')
      write(uout,'(1x,"#    Scattering rates computed in transport mode        #")')
      write(uout,'(1x,"#    (\Gamma in meV and N.B. \Gamma = 2 ImSigma)        #")')
   endif
  !write(uout, '(a)')
   write(uout,'(1x, a, f9.4, a, f10.6/)') &  
      '#  Temperature: ', tempe*ryd2ev/kelvin2ev, '  Chemical Potential: ', efermi*ryd2ev
   write(uout,'(1x, a, 3x, a, 3x, a, 3x, a)') &
      '#  ik', '  ib', '  E (eV)', '  \Gamma (meV)'
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

subroutine output_ph_cond(tempe, cond, alpha, niter_in)
   implicit none
   ! cond(6 or 9, max_step, ntempe)
   real(dp), intent(in) :: tempe(:), cond(:,:,:), alpha(:,:,:)!cond(6 or 12,:,:)
   integer, intent(in), optional :: niter_in(:)
   integer :: uout, it, ia, i, max_step, ntmp, last, ncomp_print, ncomp, start_print
   integer, allocatable :: niter(:)
   real(dp) :: cond_unit, mobi_unit, dens_unit, alpha_unit
   character(len=220) :: ctmp, label, ctmp1, ctmp2, fmt1, fmt2, fmt3, fmt4,&
                         fmt5, fmt6
   character(len=3) :: sub(9) = (/'_xx', '_xy', '_yy', '_xz', '_yz', '_zz',&
                                  '_yx', '_zx', '_zy'/)
         
  !if (system_2d) call errore('output_ph_cond','2D system is not supported',1)
   label = merge("n_c (cm^-2)", "n_c (cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   ! unit conversion of phonon contributed conductivity. check below.
   ! \Kappa = 1/(N*V*T) sum_{\nu q} (v_{\nu q }/2pi)^2*(hbar*omega)^2 *tau_{\nu q} *n_{\nu q}*(1+n_{\nu q})/kT
   ! (\sigma = q^2*nspin/(N*V) sum_{nk} v_nk*v_nk*tau_nk* f *(1-f)/kT)
   !  (N.B.: the factor alat^2 should have already included in cond at this point)
   !
   ! 1. all the quantities are in Rydberg atomic unit.
   !--------------------------------------------------------------------------
   !cond_unit = k_b/(a0^3 * k_b*T)* (Ryd/(hbar*a0^-1))^2 *Ryd^2 * (hbar/Ryd)/Ryd
   !          = k_b/a0 * Ryd/hbar = k_b*K / (hbar/Ryd) /(a0*K)
   !          = (kelvin2eV / timeunit ) /a0/K
   !convert thermal conductivity to (W/m/K) 
   !kelvin2eV * unitcharge (C*V=J) / (s) = eV * unitcharge W /a0 /K 
   != kelvin2eV * unitcharge / timeunit / bohr2ang * 1.0E+10_dp (W/m/K)
   !cond_unit = kelvin2eV * unitcharge/ timeunit / bohr2ang * 1.0E+10_dp (W/m/K)
   !(W = J/s = V*A, eV = unitchage * J= 1.6*10E-19 J)
   cond_unit = (kelvin2ev * unitcharge * 1.0E10_dp)/ (timeunit * bohr2ang)! / (2.d0 * pi)**2 ! I have devide by 2pi^2 in the
   !phonon_velocity_findiff subroutines
   !for 2D, we output conductance
   if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !alpha = 1/(V*k_b*T) * \sum (hbar omega) e v^2 tau n (1+n), 
   !unit: 1/(a0^3 Ryd) Ryd e (a0/t0)^2 t0 = e/(a0*t0)
   ! e/(a0*t0) = unitcharge/(bohr2ang * 1E-10 * timeunit) C/(m*s). 
   !see Eq. (8) in Ref. PRB 102, 245202 (2020)
   !the 2*pi is for phonon velocity from find_diff subroutine
   alpha_unit = (unitcharge * 1.0E10_dp)/ (timeunit * bohr2ang) !/ (2.d0 * pi)**2 
   !cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)

   ntmp = size(tempe)
   max_step = size(cond, 2)

   ncomp = size(cond, 1)
   ! ncomp_print designates the number of cond elements to print
   ! ncomp_print is 6 for the symmetric cond tensor (E for non-magnetic case)
   ! ncomp_print is 12 for the symmetric cond tensor (E+T for non-magnetic case)

   ! Non-magnetic case
   if( ncomp == 6 .or. ncomp == 12 ) then
      ncomp_print = 6
      start_print = 1
   ! Magnetic case
   else
      ncomp_print = 9
      start_print = 1
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
      ! syp - I changed size(niter) to size(niter_in) in the following line
      if(size(niter_in) .ne. ntmp) call errore('output_ph_cond', 'niter dimension mismatch', 1)
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
   ! array check
   if(size(cond,3).ne.ntmp) &
      call errore('output_ph_cond','array dimension mismatch',1)
   uout = find_free_unit()
   open(uout, file=trim(prefix)//'.cond_ph', status='unknown', form='formatted')

   write(uout,'(/10x,"#==========================================================#" )')

  !if(system_2d) then
  !   write(uout,'( 10x,"#                   Conductivity (1/Ohm)                   #" )')
  !   write(uout,'( 10x,"#----------------------( 2D system )-----------------------#"/)')
  !else
      write(uout,'( 10x,"#              Thermal Conductivity (W/m/K)                #" )')
      write(uout,'( 10x,"#----------------------------------------------------------#"/)')
  !endif

   write(ctmp, fmt1) ('kappa'//sub(ia), ia = 1, ncomp_print)
   write(uout, fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      last = niter(it)
      write(uout,fmt3) &
         tempe(it)*ryd2ev/kelvin2ev, 0.d0, 0.d0, &
         (cond(ia, last, it)*cond_unit, ia=start_print, start_print + ncomp_print-1)
      !for iterative approach
      if(max_step > 1) then
         write(uout,'( 10x,"#--------------------iterative process---------------------#" )')
         write(uout,fmt4) "#iter.", ctmp
         do i = 1, niter(it)
            write(uout,fmt5) '#', i, (cond(ia,i,it)*cond_unit, ia=start_print, start_print + ncomp_print-1)
         enddo
         write(uout,'( 10x,"#----------------------------------------------------------#"/)')
      endif
   enddo

   write(uout,'(/a)') 
   write(uout,'(/10x,"#==========================================================#" )')
   write(uout,'( 10x,"#                     alpha^ph_E (C/m/s)                   #" )')
   write(uout,'( 10x,"#--------------------(for semiconductor)-------------------#"/)')

   write(uout,fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      last = niter(it)
      write(uout,fmt3) &
         tempe(it)*ryd2ev/kelvin2ev, 0.d0, 0.d0, &
         (alpha(ia, last, it)*alpha_unit, ia=start_print, start_print + ncomp_print-1)
   enddo
   deallocate(niter)
   close(uout)
end subroutine output_ph_cond

!convert conductivity (electron or thermal) to the SI unit for output in the middle of running for debugging.
!subroutine convert_cond(cond, icomp, iter, itemp, cond_value, thermalornot)
subroutine convert_cond(cond_unit, thermalornot)
   implicit none
  !real(dp), intent(in) :: cond(:,:,:) !cond(6 or 9,:,:)
  !integer, intent(in) :: icomp, iter, itemp
  !real(dp), intent(out) :: cond_value
   real(dp), intent(out) :: cond_unit
   logical, intent(in), optional :: thermalornot

  !integer :: ncomp, niter, ntemp
  !real(dp) :: cond_unit
   logical :: thermal

   thermal = .false.
   if(present(thermalornot)) thermal = thermalornot
         
   if(thermal) then 
      cond_unit = (kelvin2ev * unitcharge * 1.0E10_dp)/ (timeunit * bohr2ang) !/ (2.d0 * pi)**2 
      if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)
   else
      cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)
      if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)
   endif

  !ncomp = size(cond, 1)
  !niter = size(cond, 2)
  !ntemp = size(cond, 3)

  !if( icomp < 1 .or. icomp > ncomp ) call errore('out_cond','icomp out of range',1)
  !if( iter < 1 .or. iter > niter ) call errore('out_cond','iter out of range',1)
  !if( itemp < 1 .or. itemp > ntemp ) call errore('out_cond','itemp out of range',1)

  !cond_value = cond(icomp,iter,itemp)*cond_unit

end subroutine convert_cond

subroutine output_mobility(tempe, efermi, dens, cond, niter_in, phonon_iter, phbte)
   implicit none
   ! cond(6 or 9, max_step, ntempe)
   real(dp), intent(in) :: tempe(:), cond(:,:,:) !cond(6 or 9,:,:)
   real(dp), intent(in), optional :: efermi(:), dens(:)!, alpha_ph(:,:,:) !cond(6 or 9,:,:)
   integer, intent(in), optional :: niter_in(:)
   integer, intent(in), optional :: phonon_iter
   logical, intent(in), optional :: phbte
   integer :: uout, it, ia, i, max_step, ntmp, last, ncomp_print, ncomp
   integer, allocatable :: niter(:)
   real(dp) :: cond_unit, mobi_unit, dens_unit, alpha_unit
   logical :: phbte_local 
   character(len=220) :: ctmp, label, ctmp1, ctmp2, fmt1, fmt2, fmt3, fmt4,&
                         fmt5, fmt6, phonon_iter_str
   character(len=3) :: sub(9) = (/'_xx', '_xy', '_yy', '_xz', '_yz', '_zz',&
                                  '_yx', '_zx', '_zy'/)
        
   if((present(efermi) .or. present(dens) .or. present(phonon_iter)) &
      .and. (present(phbte))) then
      if(phbte) call errore('output_mobility', 'inconsistent input', 1)
   endif

   if(present(phbte)) then
      phbte_local = phbte
   else 
      phbte_local = .false.
   endif

   label = merge("n_c (cm^-2)", "n_c (cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   if(present(phonon_iter)) then
      write(phonon_iter_str, '(I4)') phonon_iter
      phonon_iter_str = trim(adjustl(phonon_iter_str))
   endif
   ! unit conversion of cond. check below.
   ! \sigma = q^2*nspin/(N*V) sum_{nk} v_nk*v_nk*tau_nk* f *(1-f)/kT
   !  (N.B.: the factor alat^2 should have already included in cond at this point)
   !
   ! 1. q^2 is omitted in the calculation of cond(i), 
   !        we could treat q^2 = 1 because we define unitchange = 1.60 * E_19
   !        if we strictly treat q^2 = 2, we need to define unitcharge = 1.1329*E-19
   ! 2. all the quantities are in Rydberg atomic unit.
   !--------------------------------------------------------------------------
   !cond_unit = q^2* (a0)^-3* (Ryd/(hbar*a0^-1))^2 * (hbar/Ryd) / Ryd
   !          = e^2 * (hbar*a0)^-1 = e^2 *(Ryd*t0*a0)^-1
   !convert cond to (ohm*m)^-1  (ohm = V/A = V*s/C; hbar = Ryd*t0)
   !cond_unit = e / ((Ryd/e)*t0*a0) = unitcharge/((ryd2eV/e)*timeunit*(bohr2ang*1.0E-10_dp))
   !                                ~ C / (V*s*m) = 1/(Ohm*m)
   if(phbte_local)then
      cond_unit = (kelvin2ev * unitcharge * 1.0E10_dp)/ (timeunit * bohr2ang) !/ (2.d0 * pi)**2 
      alpha_unit = (unitcharge * 1.0E10_dp)/ (timeunit * bohr2ang) !/ (2.d0 * pi)**2 
      if(system_2d) alpha_unit = (alpha_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)
   else
      cond_unit = unitcharge*1.0E10_dp/(ryd2ev*timeunit*bohr2ang)
      mobi_unit = bohr2ang*bohr2ang/(timeunit*1.0E16_dp*ryd2ev)
   endif
   !for 2D, we output conductance
   if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !convert from Rydberg atomic unit to cm^2/V/s
   !mobility mu = sigma / (n*q); (n:(a0)^-3); mu = sigma * a0^3 / q 
   !mobi_unit = q* (a0)^2/(Ryd*t0)
   
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
   if (.not. phbte_local) then
      if(size(dens).ne.ntmp .or. size(cond,3).ne.ntmp) &
         call errore('output_mobility','array dimension mismatch',1)
   endif
   uout = find_free_unit()
   if(phbte_local) then
      open(uout, file=trim(prefix)//'.cond_ph', status='unknown', form='formatted')
   else
      if(present(phonon_iter)) then
         open(uout, file=trim(prefix)//'.cond_'//trim(phonon_iter_str), status='unknown', form='formatted')
      else
         open(uout, file=trim(prefix)//'.cond', status='unknown', form='formatted')
      endif
   endif

   write(uout,'(/10x,"#==========================================================#" )')

   if(phbte_local) then
      if(system_2d) then
         write(uout,'( 10x,"#               Thermal Conductivity (W/K)                 #" )')
         write(uout,'( 10x,"#----------------------( 2D system )-----------------------#"/)')
      else
         write(uout,'( 10x,"#              Thermal Conductivity (W/m/K)                #" )')
         write(uout,'( 10x,"#----------------------------------------------------------#"/)')
      endif
   else
      if(system_2d) then
         write(uout,'( 10x,"#                   Conductivity (1/Ohm)                   #" )')
         write(uout,'( 10x,"#----------------------( 2D system )-----------------------#"/)')
      else
         write(uout,'( 10x,"#                  Conductivity (1/Ohm/m)                  #" )')
         write(uout,'( 10x,"#----------------------------------------------------------#"/)')
      endif
   endif

   if(phbte_local) then
      write(ctmp, fmt1) ('kappa'//sub(ia), ia = 1, ncomp_print)
      write(uout, fmt2) '#  T (K)', '       ', trim(label), ctmp
   else
      write(ctmp, fmt1) ('sigma'//sub(ia), ia = 1, ncomp_print)
      write(uout, fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
   endif
   do it = 1, ntmp
      last = niter(it)
      if(phbte_local)then
         write(uout,fmt3) &
            tempe(it)*ryd2ev/kelvin2ev, 0.d0, 0.d0, &
            (cond(ia, last, it)*cond_unit, ia=1, ncomp_print)
      else
         write(uout,fmt3) &
            tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
            (cond(ia, last, it)*cond_unit, ia=1, ncomp_print)
      endif
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

   if(phbte_local)then
      write(uout,'(/a)') 
      write(uout,'(/10x,"#==========================================================#" )')
      if(system_2d) then
         write(uout,'( 10x,"#                       alpha^ph (C/s)                     #" )')
         write(uout,'( 10x,"#------------------------( 2D system )---------------------#"/)')
      else
         write(uout,'( 10x,"#                      alpha^ph (C/m/s)                    #" )')
         write(uout,'( 10x,"#----------------------------------------------------------#"/)')
      endif
      
      write(uout,fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
      do it = 1, ntmp
         last = niter(it)
         write(uout,fmt3) &
            tempe(it)*ryd2ev/kelvin2ev, 0.d0, 0.d0, &
            (cond(ia, last, it)*alpha_unit*tempe(it), ia=ncomp_print+1, ncomp_print + ncomp_print)
      enddo
   else
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
   endif
   deallocate(niter)
   close(uout)
end subroutine output_mobility

subroutine output_trans_coef(tempe, efermi, dens, cond, seebeck, seebeck_e, t_cond, is_mag, alpha, seebeck_ph, alpha_ph, cond_ph)
   implicit none
   real(dp), intent(in) :: tempe(:)
   logical, optional, intent(in) :: is_mag
   real(dp), intent(in), optional :: efermi(:), dens(:), cond(:,:), seebeck(:,:), seebeck_e(:,:), t_cond(:,:), alpha(:)
   real(dp), intent(in), optional :: seebeck_ph(:,:), alpha_ph(:), cond_ph(:,:)
   !
   character(len=220) :: ctmp, label, ctmp1, ctmp2, fmt1, fmt2, fmt3, fmt4,&
                         fmt5, fmt6, fmt7
   character(len=3) :: sub(9) = (/'_xx', '_xy', '_yy', '_xz', '_yz', '_zz',&
                                  '_yx', '_zx', '_zy'/)
   integer  :: uout, it, ia, ntmp, ncomp
   real(dp) :: dens_unit, seebeck_unit, cond_unit, mobi_unit, t_cond_unit, alpha_unit
   logical  :: calc_mag, ebte, phbte

   ! intialize
   ebte = .false.
   phbte = .false.

   if(present(efermi) .or. present(cond)) ebte = .true.
   if(present(alpha_ph) .or. present(cond_ph)) phbte = .true.

   if(.not. ebte .and. .not. phbte) &
      call errore('output_trans_coef', 'either ebte or phbte should be present', 1)

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
   !(Why here 1/T=k_b/(k_b*T)=k_b/(Ryd), because we define T=T*kelvin2ev/Ryd2ev in pert_param.f90 file)
   !(So T is different with K in terms of unit conversion)
   !(Specifically, T*k_b = Ryd, k_b*K * kelvin2ev/Ryd2ev = Ryd)
   t_cond_unit = kelvin2ev * (unitcharge/timeunit) * 1.0E10_dp / bohr2ang
   if(system_2d) t_cond_unit = (t_cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   alpha_unit = (unitcharge * 1.0E10_dp)/ (timeunit * bohr2ang) !/ (2.d0 * pi)**2 

   !convert from Rydberg atomic unit to cm^2/V/s
   !mobility mu = sigma / (n*q); (n:(a0)^-3); mu = sigma * a0^3 / q 
   !mobi_unit = q* (a0)^2/(Ryd*t0)
   mobi_unit = bohr2ang*bohr2ang/(timeunit*1.0E16_dp*ryd2ev)
   
   if(ebte)then
      ncomp = size(cond, 1)
   else
      ncomp = size(cond_ph, 1)
   endif
   !Create string formats based on whether to print 6 or 9 elements
   write(ctmp1,*) ncomp
   write(ctmp2,*) merge(134,90,(ncomp .eq. 9))
   write(fmt1,*) '('//trim(ctmp1)//'(5x,a8,2x))'
   write(fmt2,*) '(a8, 3x,a7, 3x,a11, 1x,a'//trim(ctmp2)//')'
   write(fmt3,*) '(f8.2, 1x, f9.5, 1x, E13.5, 1x, '//trim(ctmp1)//'(1x, E14.6))'
   write(fmt7,*) '(f8.2, 1x, f9.5, 1x, E13.5, 1x, '//trim(ctmp1)//'(1x, E14.6),1x,a33)'
   write(fmt4,*) '('//trim(ctmp1)//'(7x,a4,4x))'
   write(fmt5,*) '('//trim(ctmp1)//'(5x,a8,2x))'
   write(fmt6,*) '('//trim(ctmp1)//'(6x,a5,4x))'
   
   uout = find_free_unit()
   open(uout, file=trim(prefix)//'.trans_coef', status='unknown', form='formatted')

   if (ebte) then
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
   endif

   if(.not. calc_mag) then
      if(ebte)then
         write(uout,'(/a)') 
         write(uout,'(/10x,"#==========================================================#" )')
         write(uout,'( 10x,"#                Seebeck coefficient (mV/K)                #" )')
         write(uout,'( 10x,"#----------------------------------------------------------#"/)')
         !
         write(ctmp, fmt4) ('S'//sub(ia), ia = 1, ncomp)
         write(uout, fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
         do it = 1, ntmp
            write(uout,fmt7) &
               tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
               (seebeck(ia,it)*seebeck_unit, ia=1,ncomp), '(total)'   
            write(uout,fmt7) &
               tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
               (seebeck_e(ia,it)*seebeck_unit, ia=1,ncomp),   '(electron Peltier thermalpower)'
            if(phbte) then
               if(drag .or. ph_e) then
                  write(uout,fmt7) &
                     tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
                     (seebeck_ph(ia,it)*seebeck_unit, ia=1,ncomp),'(phonon Peltier thermalpower)'
               else
                  write(uout,fmt7) &
                     tempe(it)*ryd2ev/kelvin2ev, 0.0_dp, 0.0_dp, &
                     (seebeck_ph(ia,it)*seebeck_unit, ia=1,ncomp),'(phonon Peltier thermalpower)'
               endif
            endif
         enddo
      endif
      
      write(uout,'(/a)') 
      write(uout,'(/10x,"#==========================================================#" )')
      if(system_2d) then
         write(uout,'( 10x,"#                Thermal conductivity ( W/K )              #" )')
         write(uout,'( 10x,"#----------------------( 2D system )-----------------------#" )')
      else
         write(uout,'( 10x,"#                Thermal conductivity (W/m/K)              #" )')
      endif
      !
      write(ctmp, fmt5) ('kappa'//sub(ia), ia = 1, ncomp)
      write(uout, fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
      do it = 1, ntmp
         if(ebte)then
            write(uout,fmt7) &
               tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
               (t_cond(ia,it)*t_cond_unit, ia=1,ncomp), '(electronic thermal conductivity)'
         endif
         if(phbte)then
            if(drag .or. ph_e) then
               write(uout,fmt7) &
                  tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
                  (cond_ph(ia,it)*t_cond_unit, ia=1,ncomp), '(phononic thermal conductivity)'
            else
               write(uout,fmt7) &
                  tempe(it)*ryd2ev/kelvin2ev, 0.0_dp, 0.0_dp, &
                  (cond_ph(ia,it)*t_cond_unit, ia=1,ncomp), '(phononic thermal conductivity)'
            end if
         endif
      enddo


      write(uout,'(/a)') 
      write(uout,'(/10x,"#==========================================================#" )')
      write(uout,'( 10x,"#                       alpha (C/m/s)                      #" )')
      write(uout,'( 10x,"#--------------------(for semiconductor)-------------------#"/)')
      
      write(uout,fmt2) '#  T (K)', 'E_f(eV)', trim(label), ctmp
      do it = 1, ntmp
         if (ebte) then
            write(uout,fmt7) &
               tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
               (alpha(ia)*alpha_unit, ia=1, ncomp), '(electronic)'
         endif
         if(phbte) then
            if(drag .or. ph_e) then
               write(uout,fmt7) &
                  tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
                  (alpha_ph(ia)*alpha_unit, ia=1, ncomp), '(phononic)'
            else
               write(uout,fmt7) &
                  tempe(it)*ryd2ev/kelvin2ev, 0.0_dp, 0.0_dp, &
                  (alpha_ph(ia)*alpha_unit, ia=1, ncomp), '(phononic)'
            endif
         endif
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
