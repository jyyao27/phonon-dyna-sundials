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
module pert_param
   use pert_const, only: dp, ryd2mev, ryd2ev, kelvin2eV, bohr2ang, timeunit,&
                         tesla2ryd 
   use perturbo_autogen_param
   implicit none
   save
   logical :: spinor ! if spinor is .true., then spin degeneracy will is added in DOS.
   type, public ::  mode_info
      logical :: is_mag, is_ita, nomag_rta, mag_rta, is_trans
   end type 

   type, public :: temperature
      integer :: nt
      real(dp), pointer :: kt(:) => null()
      real(dp), pointer :: ef(:) => null()
   end type

   integer :: ntemper
   real(dp), allocatable :: temper(:)
   real(dp), allocatable :: doping(:)
   real(dp), allocatable :: efermi(:)
   real(dp), allocatable :: boltz_bfield(:,:)
   type(temperature) :: tempers
   type(mode_info) :: trans_info, trans_ph_info

contains
subroutine init_input_param()
   use io_files, only: check_tempdir
   use pert_utils, only: find_free_unit
   use qe_mpi_mod, only: meta_ionode, meta_ionode_id, world_comm, mp_bcast, stdout, mp_barrier
   implicit none
   logical :: ltmp1, ltmp2
   integer :: iunit, i, ios, j, num_temper_el
   character(len=120) :: ctmp, ctmp2, msg
   !
   CHARACTER(LEN=256), EXTERNAL :: trimcheck

   call autogen_init_input_param()

   !readin parameters
   if(meta_ionode) then
      call input_from_file()


      read(5, perturbo, err=100, iostat=ios)
100   call errore('init_input_para','reading input namelist',abs(ios))
      
      tmp_dir = trimcheck(tmp_dir)
      !do some check
      !band_min and band_max should be no less than 1
      if(band_min > band_max .or. band_min < 1 .or. band_max < 1) then
         msg = "both band_min and band_max should > 0, and band_min > band_max"
         call errore('init_input_para', trim(msg), 1)
      endif
      
      fqlist     = adjustl(fqlist)
      fklist     = adjustl(fklist)
      ftemper   = adjustl(ftemper)
      calc_mode = adjustl(calc_mode)
      polar_split = adjustl(polar_split)

      !Initialize the object trans_info
      !Store info on whether mode is RTA,ITA, magnetic ITA or magnetic RTA
      if (calc_mode(1:8) == 'trans-ph') then
         ! for phonon transport
         call initialize_trans_info(.true.)
      else
         call initialize_trans_info()
      endif

      if(any(boltz_kdim(1:3) < 1)) &
         call errore('init_input_para','illegal boltz_kdim',1)

      if( all(boltz_qdim(1:3) .eq. 0) ) then
         !by default, boltz_qdim = boltz_kdim
         boltz_qdim(:) = boltz_kdim(:)
      elseif( any(boltz_kdim(1:3) < 1) ) then
         call errore('init_input_para','boltz_qdim should all be positive!',1)
      elseif( any( mod(boltz_kdim(:), boltz_qdim(:)) .ne. 0 ) ) then
         call errore('init_input_para','boltz_qdim is incommensurate with boltz_kdim',1)
      endif

      if(boltz_emin>boltz_emax .or. boltz_de<1.0E-3_dp ) call errore &
         ('init_input_para','illegal boltz_emax, boltz_emin or boltz_de.',1)
     
      call autogen_input_param_beforeconv()
 
      boltz_emin = boltz_emin/ryd2ev
      boltz_emax = boltz_emax/ryd2ev
      if(boltz_nstep < 0) boltz_nstep = 0
      if(boltz_nstep_min < 0) boltz_nstep_min = 100
      ! from mev to Ryd
      boltz_de = boltz_de/ryd2mev
      boltz_de_ph = boltz_de_ph/ryd2mev
      phfreq_cutoff = phfreq_cutoff/ryd2mev 
      phfreq_cutoff_ph = phfreq_cutoff_ph/ryd2mev 
      delta_smear = delta_smear/ryd2mev
      delta_smear_ph = delta_smear_ph/ryd2mev
      pulse = pulse / ryd2ev
      if(trans_thr < 1.0E-16) &
         call errore('init_input_param', 'trans_thr is too small or negative', 1)

      !open temperature file
      iunit = find_free_unit()
      if(trim(prefix)  .eq. '') prefix = 'pert'

      if ( trim(yaml_fname) .eq. 'pert_output.yml' ) then
         yaml_fname = trim(prefix) // '_' // trim(calc_mode) // '.yml'
      endif

      if(trim(ftemper) .eq. '') then
         ntemper = 0
      else
         open(iunit,file=trim(ftemper),form='formatted',err=101,iostat=ios)
101      call errore('init_input_para','opening file '//trim(ftemper),ios)
         read(iunit, '(a)', iostat=ios) ctmp
         if(ios .ne. 0) call errore('init_input_param', &
            'reading ntemper in file '//trim(ftemper), 1)

         call remove_comments(ctmp)
         num_temper_el = get_num_of_words(ctmp)

         ! The correct temper file format
         if(num_temper_el == 1) then

            read(ctmp, *, iostat=ios) ntemper

         ! Old format
         else if(num_temper_el == 2) then

            read(ctmp, *, iostat=ios) ntemper, ctmp2

            write(stdout,'(/,5x,a)') repeat("!", 50)
            write(stdout,'(5x,a)') "Warning: it seems like an old-format temper file was specified."
            write(stdout,'(5x,a,/)') "The <prefix>.temper file should be formatted as follows (example for two configurations):"
            write(stdout,'(a)') "2"
            write(stdout,'(a)') "200  1.0  1.0E+18"
            write(stdout,'(a,/)') "300  2.0  1.0E+19"
            write(stdout,'(5x, a)') "To compute the Fermi energy (chemical potential) given the carrier concentration,"
            write(stdout,'(5x, a)') "set find_efermi input parameter to .true."
            write(stdout,'(5x,a,/)') repeat("!", 50)

         ! Wrong format
         else
            call errore('init_input_param', &
                 'wrong format of ' // trim(ftemper) // ' file', 1)
         end if
         !
         if(ntemper < 1) call errore('init_input_param', &
            '#. temperatures < 1 in '// trim(ftemper), 1)

         !
         allocate(temper(ntemper), doping(ntemper), efermi(ntemper))
         if(trans_info%is_mag) allocate(boltz_bfield(3,ntemper))

         ! temper in K;  efermi in eV;   doping in cm^-3
         temper = 0.0_dp; efermi = 0.0_dp; doping = 0.0_dp;
         do i = 1, ntemper
            if(.not. trans_info%is_mag) then
               read(iunit,*,iostat=ios) temper(i), efermi(i), doping(i)
            else
               read(iunit,*,iostat=ios) temper(i), efermi(i), doping(i),&
               (boltz_bfield(j,i),j=1,3)
            endif
            if(ios .ne. 0) call errore('init_input_para', &
               'reading temper in file '//trim(ftemper), 1)
         enddo
         close(iunit)
         temper = temper*kelvin2eV/ryd2ev
         efermi = efermi / ryd2ev
         if(trans_info%is_mag) boltz_bfield = boltz_bfield*tesla2ryd
         ! do the conversion in a later stage since we don't know it's 2D or 3D system.
         !for 3D from #./cm^3 to #./bohr^3, for 2D from #./cm^2 to #./bohr^2
         !doping = doping*1.0E-24_dp*(bohr2ang)**3 
      endif

      if(trans_info%is_trans) then
         !Enforce full_ite = .false. if doing magnetic calculation
         if(trans_info%is_mag) full_ite = .false.
         !Enforce boltz_nstep .eq. 0 for non-magnetic RTA
         if(trans_info%nomag_rta .and. (boltz_nstep .ne. 0)) then
            write(stdout,'(5x,a)') 'Non-magnetic RTA - resetting boltz_nstep to 0'
            boltz_nstep = 0
         endif
         !Set boltz_nstep to 50 if it is zero for any case other than non-magnetic RTA
         if(.not.(trans_info%nomag_rta) .and. (boltz_nstep .eq. 0)) then
            write(stdout,'(5x,a)') 'boltz_nstep cannot be zero for ITA or magnetic calculation - resetting to 50'
            boltz_nstep = 50
         endif
      endif

      !for coupled method, full_ite should be .true. for coupled method
      if (drag) full_ite = .true.

      !for dynamics
      if(time_step < 0.0_dp) call errore('init_input_param','negative time step',1)
      !convert to Rydberg atomic unit
      boltz_init_e0 = boltz_init_e0 / ryd2ev
      boltz_init_smear = boltz_init_smear / ryd2mev
      time_step = time_step / (timeunit*1.0E15_dp)
      hs = hs / (timeunit*1.0E15_dp)
      !from e*V/cm to Rydberg atomic unit (e*E)
      ! convert to Rydberg atomic unit: eE, eV/cm -> E_Rydberg / Bohr
      boltz_efield(1:3)  = boltz_efield(1:3)*bohr2ang/ryd2ev*1.0E-8_dp
      boltz_init_dist = trim(adjustl(boltz_init_dist))
      solver = trim(adjustl(solver))
      if(trim(solver) .ne. 'euler' .and. solver .ne. 'rk4' &
         .and. solver .ne. 'sundials') &
         call errore('init_input_param',"solver should be 'euler' or 'rk4' or 'sundials'.", 1)

      sampling = adjustl(sampling)
      if(trim(sampling) .eq. '')  sampling = 'uniform'

      !for cumulant
      if(cum_inner_emin > 0.0_dp .or. cum_inner_emax < 0.0_dp) &
         call errore('init_input_param','(spectral_emin, spectral_emax) should enclose 0.0',1)
      if(cum_outer_emin > cum_inner_emin .or. cum_outer_emax < cum_inner_emax) &
         call errore('init_input_param','outer window should enclose inner window',1)
      spectral_emin = spectral_emin / ryd2ev
      spectral_emax = spectral_emax / ryd2ev
      cum_inner_emin =  cum_inner_emin / ryd2ev 
      cum_inner_emax =  cum_inner_emax / ryd2ev
      cum_outer_emin =  cum_outer_emin / ryd2ev
      cum_outer_emax =  cum_outer_emax / ryd2ev
      ! de in meV
      cum_de   = cum_de  / ryd2mev
      if(spectral_np < 1)  spectral_np = 1
      if(cum_outer_np < 1)  cum_outer_np = 1
   endif
   call autogen_bcast_input_param()

   call mp_bcast(trans_info%is_trans, meta_ionode_id, world_comm)
   call mp_bcast(trans_info%is_mag, meta_ionode_id, world_comm)
   call mp_bcast(trans_info%is_ita, meta_ionode_id, world_comm)
   call mp_bcast(trans_info%nomag_rta, meta_ionode_id, world_comm)
   call mp_bcast(trans_info%mag_rta, meta_ionode_id, world_comm)
   !
   call mp_bcast(ntemper, meta_ionode_id, world_comm)
   call mp_bcast(find_efermi, meta_ionode_id, world_comm)
   ! Exclude phonon modes
   tempers%nt = ntemper
   if(ntemper > 0) then
      if(.not. allocated(temper)) allocate(temper(ntemper))
      if(.not. allocated(efermi)) allocate(efermi(ntemper))
      if(.not. allocated(doping)) allocate(doping(ntemper))
      call mp_bcast(temper, meta_ionode_id, world_comm)
      call mp_bcast(doping, meta_ionode_id, world_comm)
      call mp_bcast(efermi, meta_ionode_id, world_comm)

      if(trans_info%is_mag) then
         if(.not. allocated(boltz_bfield)) allocate(boltz_bfield(3,ntemper))
         call mp_bcast(boltz_bfield, meta_ionode_id, world_comm)
      endif

      if(associated(tempers%kt)) deallocate(tempers%kt)
      if(associated(tempers%ef)) deallocate(tempers%ef)
      allocate( tempers%kt(ntemper), tempers%ef(ntemper) )
      tempers%kt(:) = temper(:)
      tempers%ef(:) = efermi(:)
   endif
   call check_tempdir(tmp_dir, ltmp1, ltmp2)

end subroutine init_input_param

!convert string to lower case
subroutine lower_case(string)
   character(len=*) :: string
   integer  :: i, ic, nlen

   nlen = len(string)
   do i = 1, nlen
      ic = ichar( string(i:i) )
      if( ic >= 65 .and. ic < 90 ) string(i:i) = achar(ic+32)
   end do
end subroutine lower_case

!Initialize trans_info
subroutine initialize_trans_info(phbte)
   implicit none
   logical, intent(in), optional :: phbte
   logical :: phbte_local

   phbte_local = .false.
   if(present(phbte)) phbte_local = phbte

   trans_info%is_trans  = .false.
   trans_info%is_mag    = .false.
   trans_info%is_ita    = .false.
   trans_info%nomag_rta = .false.
   trans_info%mag_rta   = .false.

   if (phbte_local) then
      trans_info%is_trans = (trim(calc_mode) .eq. 'trans-ph-mag-rta') .or. &
                            (trim(calc_mode) .eq. 'trans-ph-mag-ita') .or. &
                            (trim(calc_mode) .eq. 'trans-ph-ita') .or. &
                            (trim(calc_mode) .eq. 'trans-ph-rta') .or. &
                            (trim(calc_mode) .eq. 'trans-ph') 
      
      trans_info%is_mag = (trim(calc_mode) .eq. 'trans-ph-mag-rta') .or. &
                          (trim(calc_mode) .eq. 'trans-ph-mag-ita')
      
      trans_info%is_ita = (trim(calc_mode) .eq. 'trans-ph-ita') .or. &
                          (trim(calc_mode) .eq. 'trans-ph-mag-ita')
      
      trans_info%nomag_rta = (trim(calc_mode) .eq. 'trans-ph-rta') .or. &
                             (trim(calc_mode) .eq. 'trans-ph')
      
      trans_info%mag_rta = (trim(calc_mode) .eq. 'trans-ph-mag-rta')
   else
      trans_info%is_trans = (trim(calc_mode) .eq. 'trans-mag-rta') .or. &
                            (trim(calc_mode) .eq. 'trans-mag-ita') .or. &
                            (trim(calc_mode) .eq. 'trans-ita') .or. &
                            (trim(calc_mode) .eq. 'trans-rta') .or. &
                            (trim(calc_mode) .eq. 'trans') 
      
      trans_info%is_mag = (trim(calc_mode) .eq. 'trans-mag-rta') .or. &
                          (trim(calc_mode) .eq. 'trans-mag-ita')
      
      trans_info%is_ita = (trim(calc_mode) .eq. 'trans-ita') .or. &
                          (trim(calc_mode) .eq. 'trans-mag-ita')
      
      trans_info%nomag_rta = (trim(calc_mode) .eq. 'trans-rta') .or. &
                             (trim(calc_mode) .eq. 'trans')
      
      trans_info%mag_rta = (trim(calc_mode) .eq. 'trans-mag-rta')
   end if

end subroutine initialize_trans_info

!> Get number of words in a string separated by one or multiplace whitespaces
pure function get_num_of_words(string) result(num)
   implicit none
   character(len=*),intent(in)  :: string !< Input string
   integer                      :: num    !< Number of words
   !local
   integer   :: ip,pos

   pos = 1
   num = 0

   do
      ip = VERIFY(string(pos:),' ')  !-- Find next non-blank
      if( ip == 0 ) exit             !-- No word found
      num = num + 1                  !-- Found something
      pos = pos + ip - 1             !-- Move to start of the word
      ip = SCAN(string(pos:),' ')    !-- Find next blank
      if( ip == 0 ) exit             !-- No blank found
      pos = pos + ip - 1             !-- Move to the blank
   end do

end function get_num_of_words

!> Remove comments from sting.
!! Comment must start by an exclamation mark.
subroutine remove_comments(string)
    implicit none
    character(len=*), intent(inout) :: string !< Input string
    integer :: pos

    pos = index(string, '!')
    if (pos /= 0) then
        string = adjustl(trim(string(:pos-1)))
    end if

end subroutine remove_comments

end module pert_param
