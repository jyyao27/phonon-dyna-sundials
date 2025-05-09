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

module pert_output
   use pert_const, only: dp, ryd2mev, ryd2ev
   use pert_param, only: prefix
   use pert_utils, only: find_free_unit
   use qe_mpi_mod, only: ionode, stdout, mp_bcast, ionode_id, inter_pool_comm
   implicit none
   
   !adapted from wannnier90/io.f90
   type timing_data
      character(len=30) :: tag
      integer  :: ncalls !number of calls
      real(dp) :: ctime  !accumulated run time
      real(dp) :: stime  !cpu time when watch is started  `
   end type timing_data

   integer, parameter :: maxtags = 100
   type(timing_data), save, private :: clocks(maxtags)
   integer, save, private :: ntags = 0

   public :: progressbar_init, progressbar, set_progress_step
   public :: output_imsigma, output_imsigma_yaml, load_imsigma, output_scatrate
   public :: stopwatch, output_timings
contains

subroutine set_progress_step(total, step, nstep, step_max)
   implicit none
   integer, intent(in) :: total  !should larger than 0
   integer, intent(out) :: step, nstep
   integer, intent(in), optional :: step_max
   integer :: nst
   
   if(total < 1) call errore("set_progress_step", "total number of tasks < 1", 1)

   if(present(step_max) .and. step_max > 0) then
      step = min(total, step_max)
   else
      do nst = 20, 1, -1
         if(mod(total, nst) < total/nst) exit
      enddo
      step = total/nst
   endif
   nstep = total/step  + merge(1, 0, mod(total,step)>0)
end subroutine set_progress_step

!> Initializes the progress bar
subroutine progressbar_init(message)
   ! What is being calculated
   character(len=*), intent(in) :: message
   character(len=53) :: bar
   character(len=36) :: msg
   integer :: i
   
   bar=" |--------------------------------------------------|"
   msg="                                    "
   ! add the message
   do i=1,min(len_trim(message),36)
       msg(i:i)=message(i:i)
   enddo
   write(stdout,'(3x, a)')  bar
   write(stdout,'(5x, a)')  msg !, bar
   !write(stdout,'(1x, a)', advance='no') msg//' |'
end subroutine progressbar_init

subroutine progressbar(iter, maximun)
   integer, intent(in) :: iter, maximun
   !
   integer :: previous, current !, i
   character(len=53) :: bar
   bar=" |--------------------------------------------------|"
   
   !previous = nint( (iter-1)*50.0_dp/dble(maximun) )
   !current  = nint( (iter  )*50.0_dp/dble(maximun) )
   previous = nint( (iter-1)*25.0_dp/dble(maximun) )
   current  = nint( (iter  )*25.0_dp/dble(maximun) )
   if(current > previous) then
      !do i = previous+1, current
      !  write(stdout,'(a)', advance='no') '='
      !enddo
      !write(stdout,'(i5,3x,a1,i5)') iter, '/', maximun
      write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/maximun, '%'
      if(iter .eq. maximun) write(stdout,'(3x, a)')  bar
      flush(stdout)
   endif
   !if(iter .eq. maximun) write(stdout,'(a)') '|'
   if(current .eq. previous .and. iter .eq. maximun) then
      write(stdout,'(8x, f7.2, a1)') 100.0_dp, '%'
      write(stdout,'(3x, a)')  bar
   endif
end subroutine progressbar

subroutine output_scatrate(tmpr, ene, ef, scatrate, pn)
   use pert_param, only: boltz_de
   implicit none
   character(len=2), intent(in) :: pn
   !imsgm(nband, nmod, ntpr, nkpt)
   real(dp), intent(in):: tmpr(:), ef(:), ene(:), scatrate(:,:)
   !local
   character(len=80) :: fname
   integer :: ne, ntemp, ndim(2)
   integer :: uout, it, ie

   ndim = shape(scatrate)
   ne = ndim(1);  ntemp = ndim(2); 
   !array size check
   if(( size(tmpr)  .ne. ntemp) .or. (size(ene) .ne. ne) .or. (size(ef) .ne. ntemp) ) &
       call errore('output_scatrate','arguments dimension mismatch.',1)

   uout = find_free_unit()
   select case (pn)
      case ('el') 
        fname = trim(prefix)//'.el_scatrate'
      case ('ph')
        fname = trim(prefix)//'.ph_scatrate'
      case default
        fname = trim(prefix)//'.el_scatrate'
   end select
   open(unit=uout, file=trim(fname), status='unknown', form='formatted')
   !write head

   select case (pn)
      case ('el')
           write(uout,1001);  write(uout,1002);
      case ('ph')
           write(uout,1012)
      case default
           write(uout,1001);  write(uout,1002);
   end select

   write(uout,1003)
   write(uout,1004) ne, ntemp
   write(uout,1005)
   do it = 1, ntemp
      write(uout,'(a)') '#'
      write(uout,1007) tmpr(it)*ryd2mev,  ef(it)*ryd2ev, boltz_de*ryd2mev
      write(uout,1011)
      write(uout,1010)
      !
      do ie = 1, ne
         write(uout,1009) it, ie,  &
            ene(ie)*ryd2ev, scatrate(ie, it)*ryd2mev
      enddo
   enddo
   close(uout)
   return
1001 format('# Electron Self-Energy in the Migdal Approximation as a function of Energy#')
1002 format('#        ( only for bands within [band_min, band_max] )        #')
1003 format('#--------------------------------------------------------------#')
1004 format('# NO.E:',i7,3x, 'NO.T:',i4,3x)
1005 format('it', 2x, 'ie', 2x, 'E(ibnd)(eV)',4x,'Im(Sigma)(meV)/dE')
1007 format('# Temperature(T)=',f10.5,' meV;',2x,'Chem.Pot.(mu) (ev)=',f10.5,' eV','dE (meV)',f10.5)
1009 format(i3, 2x, i6, 2x, f20.16, 2x, es23.16)
1010 format('#------------------------------------------------------------')
1011 format('#============================================================')
1012 format('# Phonon (Imaginary) Self-Energy in the Migdal Approximation #')

end subroutine output_scatrate


subroutine output_imsigma(msolve, tmpr, ef, enk, imsgm, pn, suffix)
   implicit none
   logical, intent(in) :: msolve
   character(len=2), intent(in) :: pn
   !imsgm(nband, nmod, ntpr, nkpt)
   real(dp), intent(in):: tmpr(:), ef(:), enk(:,:), imsgm(:,:,:,:)
   character(len=*), intent(in), optional :: suffix
   !local
   character(len=80) :: fname
   integer :: nk, ntemp, numb, nmodes, ndim(4)
   integer :: uout, it, imode, ib, ik

   ndim = shape(imsgm)
   numb = ndim(1);  nmodes = ndim(2);  ntemp = ndim(3);  nk = ndim(4)
   !array size check
   if( size(tmpr)  .ne. ntemp .or. size(ef) .ne. ntemp .or. &
       size(enk,2) .ne. nk .or. size(enk,1) .ne. numb ) &
       call errore('output_imsigma','arguments dimension mismatch.',1)

   uout = find_free_unit()
   if(present(suffix)) then
      fname = trim(prefix)//'.'//trim(adjustl(suffix))
   else
      select case (pn)
         case ('el') 
           fname = trim(prefix)//'.imsigma'
         case ('ph')
           fname = trim(prefix)//'.phonlife'
         case default
           fname = trim(prefix)//'.imsigma'
      end select
   endif
   if(msolve) fname = trim(fname)//'_mode'
   open(unit=uout, file=trim(fname), status='unknown', form='formatted')
   !write head

   select case (pn)
      case ('el')
           write(uout,1001);  write(uout,1002);
      case ('ph')
           write(uout,1012)
      case default
           write(uout,1001);  write(uout,1002);
   end select

   write(uout,1003)
   write(uout,1004) nk, numb, ntemp, merge(nmodes, 1, msolve)
   do it = 1, ntemp
      write(uout,'(a)') '#'
      write(uout,1007) tmpr(it)*ryd2mev,  ef(it)*ryd2ev
      write(uout,1011)
      if(msolve) then
         write(uout,1005) 
      else
         write(uout,1006)
      endif
      write(uout,1010)
      !
      do ik = 1, nk
         do ib = 1, numb
            if(msolve) then
               do imode = 1, nmodes
                  write(uout,1008) it, ik, ib, &
                  enk(ib,ik)*ryd2ev, imode, imsgm(ib,imode,it,ik)*ryd2mev
               enddo
            else
               write(uout,1009) it, ik, ib, &
                  enk(ib,ik)*ryd2ev, sum(imsgm(ib,:,it,ik))*ryd2mev
            endif
         enddo
         if(msolve) write(uout,1010)
      enddo
   enddo
   close(uout)
   return
1001 format('# Electron (Imaginary) Self-Energy in the Migdal Approximation #')
1002 format('#        ( only for bands within [band_min, band_max] )        #')
1003 format('#--------------------------------------------------------------#')
1004 format('# NO.k:',i7,3x, 'NO.bands:',i4,3x, 'NO.T:',i4,3x, 'NO.modes:',i4)
1005 format('#',1x,'it',5x,'ik',3x,'ibnd',4x,'E(ibnd)(eV)',2x,'imode',4x,'Im(Sigma)(meV)')
1006 format('#',1x,'it',5x,'ik',3x,'ibnd',4x,'E(ibnd)(eV)',5x,'Im(Sigma)(meV)')
1007 format('# Temperature(T)=',f10.5,' meV;',2x,'Chem.Pot.(mu)=',f10.5,' eV')
1008 format(i3, 2x, i6, 2x, i4, 2x, f12.6, 2x, i4, 2x, es23.16)
1009 format(i3, 2x, i6, 2x, i4, 2x, f20.16, 2x, es23.16)
1010 format('#------------------------------------------------------------')
1011 format('#============================================================')
1012 format('# Phonon (Imaginary) Self-Energy in the Migdal Approximation #')

end subroutine output_imsigma


subroutine output_imsigma_yaml(tmpr, ef, enk, imsgm, nkvec, kvectors)
   use yaml_utils, only: ymlout
   implicit none
   !imsgm(nband, nmod, ntpr, nkpt)
   real(dp), intent(in) :: tmpr(:), ef(:), enk(:,:), imsgm(:,:,:,:)
   integer, intent(in)  :: nkvec 
   real(dp), intent(in) :: kvectors(3, nkvec) 
   !local
   character(len=80) :: fname
   integer :: nk, ntemp, numb, nmodes, ndim(4)
   integer :: uout, it, imode, ib, ik, iv
   character(len=6), external :: int_to_char

   ndim = shape(imsgm)
   numb = ndim(1);  nmodes = ndim(2);  ntemp = ndim(3);  nk = ndim(4)
   !array size check
   if( size(tmpr)  .ne. ntemp .or. size(ef) .ne. ntemp .or. &
       size(enk,2) .ne. nk .or. size(enk,1) .ne. numb ) &
       call errore('output_imsigma','arguments dimension mismatch.',1)

   !write head
   write(ymlout, '(/,a)') 'imsigma:'

   write(ymlout, '(/,3x,a)') 'temperature units: meV'

   write(ymlout, '(/,3x,a)') 'energy units: eV'

   write(ymlout, '(/,3x,a)') 'Im(Sigma) units: meV'

   write(ymlout, '(/,3x,a)') 'chemical potential units: eV'

   write(ymlout, '(/,3x,a)') 'k-point coordinate units: crystal'

   write(ymlout, '(/,3x,a,i8)') 'number of k-points:', nk

   write(ymlout, '(/,3x,a)') '# Band count starts from band_min'
   write(ymlout, '(3x,a,i4)') 'number of bands:', numb

   write(ymlout, '(/,3x,a,i4)') 'number of phonon modes:', nmodes

   write(ymlout, '(/,3x,a)') '# Cofiguration means one line in the temper file:'
   write(ymlout, '(3x,a)') '# temperature efermi concentration'
   write(ymlout, '(3x,a,i4)') 'number of configurations:', ntemp

   write(ymlout, '(/,3x,a)') 'k-point coordinates:'
   do iv = 1, nkvec
      write(ymlout, '(6x,"-",1x, "[", 3(f10.5,",",2x), "]" )') kvectors(:,iv)
   enddo

   write(ymlout, '(/,3x,a)') 'energy:'

   write(ymlout, '(/,6x,a)') 'band index:'

   do ib = 1, numb

      write(ymlout, '(/,9x,a,a)') trim(int_to_char(ib)), ':'

      do ik = 1, nk
         write(ymlout,'(12x,"-", f16.10)') enk(ib,ik)*ryd2ev
      enddo

   enddo

   write(ymlout, '(/,3x,a)') 'configuration index:'

   do it = 1, ntemp

      write(ymlout, '(/,6x,a,a)') trim(int_to_char(it)), ':'

      write(ymlout,'(/,9x, a, f10.5)') 'temperature:', tmpr(it)*ryd2mev

      write(ymlout,'(/,9x, a, f15.10)') 'chemical potential:', ef(it)*ryd2ev

      write(ymlout, '(/,9x,a)') 'band index:'

      do ib = 1, numb

         write(ymlout, '(/,12x,a,a)') trim(int_to_char(ib)), ':'

         write(ymlout, '(/,15x,a)') 'Im(Sigma):'

         write(ymlout, '(/,18x,a)') 'total:'
         do ik = 1, nk
            write(ymlout,'(21x,"-", 1x, es23.16)') sum(imsgm(ib,:,it,ik))*ryd2mev
         enddo

         write(ymlout, '(/,18x,a)') 'phonon mode:' 
         do imode = 1, nmodes
            
            write(ymlout, '(/,21x,a,a)') trim(int_to_char(imode)), ':'
            do ik = 1, nk
               write(ymlout,'(24x,"-", 1x, es23.16)') imsgm(ib,imode,it,ik)*ryd2mev
            enddo

         enddo

      enddo ! ib loop
   enddo ! it loop

end subroutine output_imsigma_yaml


subroutine output_gruneisen(wq, grun)
   implicit none
   real(dp), intent(in):: wq(:,:), grun(:,:)
   !local
   character(len=80) :: fname
   integer :: uout, ndim(2), nk, nmodes
   integer :: ik, imode

   ndim = shape(grun)
   nmodes = ndim(1); nk = ndim(2)
   !array size check
   if ( size(wq,1) .ne. nmodes .or. size(wq,2) .ne. nk ) &
        call errore('output_gruneisen','arguments dimension mismatch.',1)

   uout = find_free_unit()
   fname = trim(prefix)//'.gruneisen'
   open(unit=uout, file=trim(fname), status='unknown', form='formatted')

   write(uout,1001);  write(uout,1002);
   write(uout,1003) nk, nmodes
   write(uout,1004)

   do ik = 1, nk
      do imode = 1, nmodes
      write(uout,1005) ik, imode, wq(imode,ik)*ryd2mev, grun(imode,ik)
      enddo
   enddo
   write(uout,1006)
   close(uout)
   return
1001 format('# Gruneisen parameter for each phonon mode #')
1002 format('#--------------------------------------------------------------#')
1003 format('# NO.k:',i7,3x, 'NO.modes:',i4)
1004 format('#',1x,'ik',3x,'imode',4x,'E(imode)(meV)', 4x,'Gruneisen')
1005 format(i4, 3x, i3, 3x, E12.6, 5x, E12.6)
1006 format('#------------------------------------------------------------')
end subroutine output_gruneisen

subroutine load_imsigma(imsigma, success, enk, suffix)
   implicit none
   !imsigma(numb, numk, ntpr), enk(numb,numk)
   logical, intent(out)  :: success
   real(dp), intent(out) :: imsigma(:,:,:)
   real(dp), intent(out), optional :: enk(:,:)
   character(len=*), intent(in), optional :: suffix
   !
   logical :: has_file
   real(dp) :: wtmp, rtmp
   character(len=80) :: fname, ctmp0, ctmp1, ctmp2, ctmp3 
   integer :: iunit, i, ir_nk, tmp_nb, ntpr, ndim(3), it, ik, ib, itmp(3), ios
   
   success = .false.
   imsigma = 0.0_dp
   ndim = shape(imsigma)
   if( present(enk) ) then
      if( ndim(1).ne.size(enk,1) .or. ndim(2).ne.size(enk,2) )&
      call errore('load_imsigma','arguments dimension mismatch',1)
   endif
   
   if(present(suffix)) then
      fname = trim(prefix)//'.'//trim(adjustl(suffix))
   else
      fname = trim(prefix)//".imsigma"
   endif
   inquire(file=fname, exist=has_file)
   if(.not. has_file ) return

   if(ionode) then
      iunit = find_free_unit()
      !write(stdout,'(1x,a)') "Reading imsigma from file: "//trim(fname)
      open(iunit,file=fname,status='old',form='formatted',err=101,iostat=ios)
      !skip header
      do i = 1, 3
         read(iunit,*) ctmp1
      enddo
      read(iunit, *)  ctmp0, ctmp1, ir_nk,  ctmp2, tmp_nb,  ctmp3, ntpr
      if(ndim(1).ne.tmp_nb .or. ndim(2).ne.ir_nk .or. ndim(3).ne.ntpr) &
         call errore('load_imsigma','incompatible data in imsigma file',1)
      !load data
      do it = 1, ntpr
         !skip header
         do i = 1, 5
            read(iunit, *) ctmp1
         enddo
         !
         do ik = 1, ir_nk
         do ib = 1, tmp_nb
            read(iunit,*) itmp(3), itmp(1), itmp(2), wtmp, rtmp
            if(itmp(3) .ne. it) &
               call errore('load_imsigma','Temperature index mismatch.',1)
            if(present(enk) .and. it .eq. 1)  enk(ib,ik) = wtmp/ryd2ev
            ! prefix.imsigma file is the imaginary part of self-energy in unit of meV.
            ! convert to Ryd
            imsigma(ib,ik,it) = rtmp/ryd2mev
         enddo
         enddo
      enddo
      close(iunit)
   endif
   if( present(enk) ) call mp_bcast(enk, ionode_id, inter_pool_comm)
   call mp_bcast(imsigma, ionode_id, inter_pool_comm)
   success = .true.
   return
101   call errore('load_imsigma','opening file '//trim(fname),abs(ios))
end subroutine load_imsigma

!adapted from wannier90/io.f90
subroutine stopwatch(ctag, mode)
   implicit none
   character(len=*), intent(in) :: ctag
   character(len=*), intent(in) :: mode ! 'start' or 'stop'

   integer :: i
   real(dp) :: t
   
   call cpu_time(t)

   select case(mode)
   case('start')
      !tag exists
      do i = 1, ntags
         if(clocks(i)%tag .eq. ctag) then
            clocks(i)%stime = t
            clocks(i)%ncalls = clocks(i)%ncalls + 1
            return
         endif
      enddo
      !new tags
      ntags = ntags + 1
      if(ntags > maxtags) call errore('stopwatch','exceed max number of tags',1)
      clocks(ntags)%tag    = ctag
      clocks(ntags)%stime  = t
      clocks(ntags)%ctime  = 0.0_dp
      clocks(ntags)%ncalls = 1
   case('stop')
      do i = 1, ntags
         if(clocks(i)%tag .eq. ctag) then
            clocks(i)%ctime = clocks(i)%ctime + t - clocks(i)%stime
            return
         endif
      enddo
      !tag not found
      call errore('stopwatch',trim(ctag)//' not found',-1)
   case default
      call errore('stopwatch','Mode '//trim(mode)//' not recognised',-1)
   end select
end subroutine stopwatch

subroutine output_timings()
   implicit none
   integer :: i
   
   write(stdout,'(/1x,a)') 'Timing Summary:'
   write(stdout,'(1x,a)') repeat('-',53)
   write(stdout,1001) 
   do i = 1, ntags
      write(stdout,1002) clocks(i)%tag, clocks(i)%ncalls, clocks(i)%ctime
   enddo
1001 format(1x,4x,'Tag',23x,1x,'#.calls',2x, 5x,'Times(s)')
1002 format(1x,a30, 1x,i7, 2x,f13.3)
end subroutine output_timings

end module 
