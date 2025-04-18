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

module boltz_trans_mod
   use pert_const, only: dp
   use pert_param, only: calc_mode, drag
   use pert_output, only: stopwatch
   use boltz_grid, only: grid
   use pert_utils, only: find_free_unit, fermi, mfermi_deriv, mat_inv3, mat_inv2, &
                        mbose_deriv
   use qe_mpi_mod, only: mp_split_pools, mp_sum, mp_barrier, inter_pool_comm, ionode

   private
   public :: trans_dos_calc !compute density of states
   public :: trans_tdf_calc
   public :: trans_density_calc
   public :: trans_cond_calc
   public :: trans_cond_calc_ph
   public :: extract_trans_coeff
   public :: Tconduct
   public :: trans_impose_kelveinonsager
contains


! Straightforward implementation of the thermal conductivity as an integral
! over the whole Brillouin zone in terms of frequencies, velocities and F_n.
subroutine TConduct(nptk,nbands,T,vol,omega,velocity,F_n)
   use pert_const, only: kelvin2ev, ryd2ev, ryd2mev, pi, bohr2ang, ryd2hz_h, ryd2thz_h, ryd2thz_hbar
   use pert_data, only: alat
   implicit none

   integer, intent(in) :: nptk,nbands
   real(kind=dp),intent(in) :: T, vol
   real(kind=dp),intent(in) :: omega(Nbands,nptk),velocity(3,Nbands,nptk),F_n(3,Nbands,nptk)
   
   real(kind=dp) :: ThConductivity(3,3,nbands)
   real(kind=dp) :: ThConductivityMode(3,3,nbands,nptk)

   !real(kind=dp),parameter :: pi=3.141592653589793238_dp
   real(kind=dp),parameter :: kb=1.380648813e-23_dp ! J/K
   real(kind=dp),parameter :: hbar=1.05457172647e-22_dp ! J*THz

   real(kind=dp) :: omega_loc(nbands,nptk),velocity_loc(3,nbands,nptk),F_n_loc(3,nbands,nptk), T_loc
   real(kind=dp) :: vol_loc

   real(kind=dp) :: fBE,tmp(3,3)
   integer :: ii,jj,dir1,dir2

   T_loc = T/kelvin2eV*ryd2ev ! to Kelvin
   velocity_loc = velocity*alat*bohr2ang*ryd2hz_h*1.0E-13_dp ! to Km/s
   omega_loc = omega*ryd2thz_hbar ! to THz (rad/picosecond)
   F_n_loc = F_n*alat*bohr2ang*ryd2hz_h*1.0E-13_dp
   ! F = v * ometa/ scatter_rates = &
   ! (velocity*alat*bohr2ang*ryd2hz_h*1.0E-13_dp)*(omega*ryd2thz_hbar) / (scatter_rates*ryd2thz_hbar)
   vol_loc = vol * bohr2ang**3 * 1.0E-3_dp
   ! vol is in bohr^3, convert to nm^3, which is used in ShengBTE

   ThConductivity=0.0_dp
   ThConductivityMode=0.0_dp
   !$OMP PARALLEL DO default(none) collapse(2) schedule(static) &
   !$OMP & shared(nptk,nbands,omega_loc,velocity_loc,F_n_loc,ThConductivityMode,T_loc) &
   !$OMP & private(jj,ii,dir1,dir2,tmp,fBE) reduction(+:ThConductivity)
   do ii=2,nptk
      do jj=1,Nbands
         do dir1=1,3
            do dir2=1,3
               tmp(dir1,dir2)=velocity_loc(dir1,jj,ii)*F_n_loc(dir2,jj,ii)
            end do
         end do
         fBE=1.0_dp/(exp(hbar*omega_loc(jj,ii)/Kb/T_loc)-1.0_dp)
         ThConductivityMode(:,:,jj,ii)=fBE*(fBE+1.0_dp)*omega_loc(jj,ii)*tmp
         ThConductivity(:,:,jj)=ThConductivity(:,:,jj)+ThConductivityMode(:,:,jj,ii)
      end do
   end do
   !$OMP END PARALLEL DO
   write(*,*) 'volumn',vol
   write(*,*) 'omega_loc',omega_loc(6,:)
   ThConductivity=1.0e21_dp*hbar**2*ThConductivity/(kB*T_loc*T_loc*Vol_loc*nptk)
   ThConductivityMode=1.0e21_dp*hbar**2*ThConductivityMode/(kB*T_loc*T_loc*Vol_loc*nptk)
   if (ionode) write(*,*) 'Thermal Conductivity by ShengBTE (W/mK):', sum(ThConductivity,dim=3)
    
end subroutine TConduct

! compute dos using tetrahedra integration, spin is taken into account here
subroutine trans_dos_calc(kg, ene, dos)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: ene(:) 
   real(dp), intent(out) :: dos(:)
   !local variables
   integer :: nstep, nb, it, ic, ie, it_st, it_end
   real(dp) :: fnk(kg%numb,4), et(kg%numb,4)
   real(dp), allocatable :: dos_t(:)
   ! get the size of input array
   nstep = size(ene);   nb = kg%numb
   dos(1:nstep) = 0.0_dp;  fnk = 1.0_dp

   !mpi parallel: distribute tetras over pools
   call mp_split_pools(kg%num_tet, it_st, it_end)
!$omp parallel default(shared) private(it, ic, et, dos_t, ie)
   allocate( dos_t(nstep) );  dos_t = 0.0_dp
!$omp do schedule(guided)
   do it = it_st, it_end
      ! prepare data for the four vertex.
      do ic = 1, 4
         ! tet(2, ic, it): idx of kgrid. full klist
         et(:,ic) = kg%enk(:, kg%tetra(2,ic,it))
      enddo
      !compute the contribution from current tetra, and add it to dos_t
      !NOTE: dos_t is inout. initialization is required for dos_t
      call tetra_int(nb, et, fnk, nstep, ene, dos_t)
   enddo
!$omp end do nowait
   !collect results from each thread
   do ie = 1, nstep
!$omp atomic
      dos(ie) = dos(ie) + dos_t(ie)
   enddo
   deallocate(dos_t)
!$omp end parallel

   !DOS in the units of 'states/Ryd/unit.cell.'
   !tweight = 1/tot_tet_fbz = V_T/V_G, (spin degeneracy is not included here)
   dos = dos * kg%tweight
   ! collect results from different pools
   call mp_sum(dos, inter_pool_comm)
   !call mp_barrier(inter_pool_comm)
end subroutine trans_dos_calc

!compute transport distribution function (TDF)
! TDF_ij(E) = spinor/(N*V) sum_{nk} vnk * vnk * tau * delta(E-e_nk)
! the unit of velocity is Ryd * Bohr * alat, vnk = 1/hbar*dE/dk.
! However, here we does not include spinor and volume V. Be careful!
! syp @04252024 - I think this function applies to both e-BTE and ph-BTE
! but for phonons, the tetrahedron weights are always zero,
! so I use uniform grid integration for phonons
subroutine trans_tdf_calc(kg, ene, mfd, tdf, calc_mag)
   use boltz_utils, only: velprod
   implicit none
   type(grid), intent(in) :: kg
   logical, intent(in) :: calc_mag
   !mfd: mean free displacement = vnk*tau, mfd(3 or 6,numb,kg%nk)
   real(dp), intent(in) :: ene(:), mfd(:,:,:)
   real(dp), intent(out) :: tdf(:,:) !tdf(6 or 12,nstep)
   !local variables
   integer :: nb, nstep, it, ic, ik, ib, ia, ie, it_st, it_end, ncomp
   real(dp) :: et(kg%numb, 4)
   real(dp), allocatable :: tdf_t(:,:), fnk(:,:,:) !fnk(kg%numb, 4, 6 or 12)

   nb = kg%numb;  nstep = size(ene);   ncomp = size(tdf, 1)
   ! array check
   if((.not. calc_mag) .and. (size(mfd, 1) .ne. (ncomp/2)))&
      call errore('trans_tdf_calc', 'mismatch dimension in input arguments!',1)

   ! Currently, magnetic+thermal transport is not implemented, so we 
   ! expect mfd to have 3 as first dimension and ncomp as 9
   if((calc_mag) .and. (size(mfd, 1) .ne. (ncomp/3)))&
      call errore('trans_tdf_calc', 'mismatch dimension in input arguments!',1)

   if( size(mfd, 2) .ne. nb ) &
      call errore('trans_tdf_calc','array dimension mismatch', 1)

   tdf(:,:) = 0.0_dp
   !mpi parallel: distribute tetras over pools
   call mp_split_pools(kg%num_tet, it_st, it_end)
!$omp parallel default(shared) private(it, ic, et, ik, ib, fnk, ia, tdf_t, ie)
   allocate( tdf_t(nstep, ncomp), fnk(kg%numb, 4, ncomp) );  tdf_t = 0.0_dp
!$omp do schedule(guided)
   do it = it_st, it_end
      ! prepare data for the four vertex.
      do ic = 1, 4
         ! tet(2, ic, it): idx of kgrid. full klist
         ik = kg%tetra(2,ic,it)
         et(:,ic) = kg%enk(:,ik)
         ! get vnk_a * mfd_b - > TDF_a,b 
         do ib = 1, nb
            call velprod2( kg%vnk(1:3,ib,ik), mfd(:,ib,ik), fnk(ib,ic,:) )
         enddo
      enddo
      do ia = 1, ncomp
         call tetra_int(nb, et, fnk(:,:,ia), nstep, ene, tdf_t(:,ia))
      enddo
   enddo
!$omp end do nowait

   !collect results from each thread
   do ie = 1, nstep
      do ia = 1, ncomp
!$omp atomic
         tdf(ia, ie) = tdf(ia, ie) + tdf_t(ie, ia)
      enddo
   enddo
   deallocate(tdf_t, fnk)
!$omp end parallel

   !NOTE: volume and spinor is not included here.!!!
   !tweight = 1/tot_tet_fbz = V_T/V_G,
   tdf = tdf * kg%tweight
   ! collect results from different pools
   call mp_sum(tdf, inter_pool_comm)
   return
   !
   contains
      subroutine velprod2(v1, vv2, prod2)
         implicit none
         real(dp), intent(in) :: v1(3)
         real(dp), intent(in) :: vv2(:)
         real(dp), intent(out) :: prod2(:)
         !
         integer :: nsize, i

         nsize = size(vv2)
         if(size(prod2) .eq. 9) then !If magnetic field then
               prod2(1:6) = velprod(v1, vv2(1:3)) 
               prod2(7)   = v1(2)*vv2(1) !yx
               prod2(8)   = v1(3)*vv2(1) !zx
               prod2(9)   = v1(3)*vv2(2) !zy
         else
            do i = 0, nsize-1, 3
               prod2( (i*2+1):(i*2+6) ) = velprod(v1, vv2( (i+1):(i+3) ))
            enddo
         endif
      end subroutine velprod2
   ! 
end subroutine trans_tdf_calc


!compute carrier concentration: number of carriers per unit cell
!  or determine chemical potential at given concentration.
subroutine trans_density_calc(kg, ene, tmpr, ef, hole, ncol, vol, dens, find_ef, dos)
   implicit none
   type(grid), intent(in) :: kg
   !hole: true -> hole carrier; false -> electron carrier
   logical, intent(in) :: hole, ncol  !ncol: non-colliner
   real(dp), intent(in) :: ene(:), tmpr(:), vol !tmpr: temperature, volume
   !dens: density in Rydberg atomic unit: #.carrier / bohr^3
   real(dp), intent(inout) :: ef(:), dens(:)
   logical, intent(in), optional :: find_ef !if true, calc ef from input density
   real(dp), intent(out), optional :: dos(:)  ! density of states
   ! local variables
   logical :: less
   integer :: i, nestep, ntmpr
   real(dp) :: mid_dens, ist, iend, dist, mid, estep
   real(dp), allocatable :: dos_tmp(:), occup(:)

   call stopwatch('trans_density_calc','start')
   nestep = size(ene)
   allocate( dos_tmp(nestep), occup(nestep) )
   ntmpr = size(tmpr)
   estep = (ene(nestep)-ene(1)) / (nestep-1.0_dp) ! energy step
   ! some check
   if( (size(dens) .ne. ntmpr) .or. (size(ef) .ne. ntmpr)) &
      call errore('trans_density_calc','array dimension mismatch', 1)
   if( (ene(nestep)-ene(1)) < 1.0E-4_dp ) &
      call errore('trans_density_calc','energy range is too small', 1)
   ! compute density of states
   call trans_dos_calc(kg, ene, dos_tmp)
   !account for spin-degeneracy
   dos_tmp = dos_tmp * merge(1.0_dp, 2.0_dp, ncol)
   !  
   if( present(dos) ) then
      if( size(dos) .ne. nestep ) &
         call errore('trans_density_calc',"argument 'dos': mismatch dimension.", 1)
      !DOS in the units of 'states/Ryd/unit.cell.'
      dos(1:nestep) = dos_tmp(1:nestep)
   endif
   ! find ef with density close to dens
   if( present(find_ef) ) then
      if (find_ef) then
         do i = 1, ntmpr
            ! skip if dens is unavailable
            if( dens(i) < 1.0E-15_dp ) cycle
            ! do binary search
            ist = ene(1);  iend = ene(nestep);  dist = iend - ist;
            do while ( dist > 0.5E-4_dp )
               mid = (ist + iend)*0.5_dp
               occup = fermi(tmpr(i), merge(mid-ene, ene-mid, hole))
               mid_dens = dot_product(dos_tmp, occup) * estep / vol
               !if density matches, exit the loop
               if( abs(mid_dens-dens(i))/dens(i) < 1.0E-5_dp ) exit
               !which direction to search
               less = .not. hole .and. (dens(i) > mid_dens)
               less = less .or. (hole .and. (dens(i) < mid_dens) )
               if(less) then
                  ist = mid
               else
                  iend = mid
               endif
               dist = iend - ist
            enddo
            ef(i) = mid
         enddo
      endif
   endif
   ! compute carrier density: int_dos(E)*f(E)*dE 
   do i = 1, ntmpr
      occup = fermi(tmpr(i), merge(ef(i)-ene, ene-ef(i), hole))
      ! carrier density: #. carrier / bohr^3
      dens(i) = dot_product(dos_tmp, occup) * estep / vol
   enddo
   deallocate(dos_tmp, occup)
   !
   call stopwatch('trans_density_calc','stop')
end subroutine trans_density_calc

!calculate conductivity and mobility using tetrahedron integration. 
subroutine trans_cond_calc_ph(kg, tmpr, vol, mfd, cond, alpha, calc_mag)
   use boltz_utils, only: velprod
   implicit none
   type(grid), intent(in) :: kg
   logical,  intent(in) :: calc_mag
   real(dp), intent(in) :: tmpr, mfd(:,:,:), vol
   real(dp), intent(out) :: cond(:)
   real(dp), intent(out) :: alpha(:)
   ! local variables
   integer :: ncomp, it_st, it_end, ik, icomp, ib, nb
   real(dp), allocatable :: fnk(:,:)
   real(dp) :: cond1

   call stopwatch('trans_cond_calc_ph','start')

   !ncomp here is not the same as ncomp in transport.f90
   !ncomp = 6 : Just electric field
   !ncomp = 9 : Both electric and magnetic fields
   !ncomp = 12: Electric and T fields
   ncomp = size(mfd, 1) * merge(3,2,calc_mag)  
   nb = kg%numb

   cond = 0.0_dp; alpha = 0.0_dp

   !mpi parallel: distribute tetras over pools
   call mp_split_pools(kg%nk, it_st, it_end)
!$omp parallel default(shared) private(ik, ib, icomp, fnk)
   allocate( fnk(ncomp, kg%numb) ); fnk = 0.0_dp
!$omp do schedule(guided)
   ! calculate fnk for each components using integration over uniform grid not tetrahedra
   do ik = it_st, it_end
      do ib = 1, nb
         call velprod2( kg%vnk(1:3,ib,ik), mfd(:,ib,ik), fnk(:,ib) )
      enddo
      do icomp = 1, ncomp
         do ib = 1, nb
!$omp atomic
          cond(icomp) = cond(icomp) + &
              fnk(icomp,ib) * kg%enk(ib,ik) * ( mbose_deriv(tmpr, kg%enk(ib,ik)))
         enddo 
      enddo
   enddo
!$omp end do
   deallocate(fnk)
!$omp end parallel
   
   call mp_sum(cond, inter_pool_comm)
   !sypeng @ 09242024: I suspect it's wrong here because alpha and thermal conductivity should 
   ! be derived by different components of fnk.
   alpha = - cond / vol * kg%kweight
   cond = cond / (tmpr * vol) * kg%kweight

   call stopwatch('trans_cond_calc_ph','stop')
   return
   !
   contains
      subroutine velprod2(v1, vv2, prod2)
         implicit none
         real(dp), intent(in) :: v1(3)
         real(dp), intent(in) :: vv2(:)
         real(dp), intent(out) :: prod2(:)
         !
         integer :: nsize, i

         nsize = size(vv2)
         if(size(prod2) .eq. 9) then !If magnetic field then
               prod2(1:6) = velprod(v1, vv2(1:3)) 
               prod2(7)   = v1(2)*vv2(1) !yx
               prod2(8)   = v1(3)*vv2(1) !zx
               prod2(9)   = v1(3)*vv2(2) !zy
         else
            do i = 0, nsize-1, 3
               prod2( (i*2+1):(i*2+6) ) = velprod(v1, vv2( (i+1):(i+3) ))
            enddo
         endif
      end subroutine velprod2
end subroutine trans_cond_calc_ph

subroutine trans_impose_kelveinonsager(kg, ene, tmpr, ef, ncol, vol, ph_alpha_T, mfd)
   use pert_const, only: dp, ryd2ev, bohr2ang, kelvin2ev, ryd2mev, unitcharge,&
                         timeunit, tesla2ryd, pi
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: tmpr, ef, ene(:), vol
   logical,  intent(in) :: ncol
   real(dp), intent(in) :: ph_alpha_T(:)
   real(dp), intent(inout) :: mfd(:,:,:)
   ! local variables
   integer :: i, nestep, ncomp, ik, ib
   real(dp) :: estep
   real(dp) :: lambda(6)
   real(dp), allocatable :: cond_seebeck(:)
   real(dp), allocatable :: tdf_tmp(:,:)

   real(dp), allocatable :: mfd_diff(:,:,:)
   real(dp), allocatable :: mfd_drag(:,:,:)

   real(dp) :: seebeck_unit, cond_unit

   seebeck_unit = kelvin2ev * 1.0E3_dp
   cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)

   nestep = size(ene)

   !ncomp = 6: Electric and T fields
   ncomp = size(mfd, 1) 

  !if(ionode) write(*,*) 'ph_alpha_T', ph_alpha_T
   if(ncomp .ne. size(ph_alpha_T)) &
      call errore('trans_impose_kelveinonsager', 'mismatch dimension in input arguments!',1)
   allocate( tdf_tmp(ncomp, nestep) )
   allocate( mfd_diff(ncomp/2, kg%numb, kg%nk) )
   allocate( mfd_drag(ncomp/2, kg%numb, kg%nk) )
   allocate( cond_seebeck(ncomp) )
   estep = (ene(nestep)-ene(1)) / (nestep-1.0_dp)
   
   do ik = 1, kg%nk
      do ib = 1, kg%numb
         do i = 1, ncomp/2
            mfd_diff(i, ib, ik) = mfd(i, ib, ik) * (kg%enk(ib,ik) - ef)
         enddo
      enddo
   enddo

   mfd_drag = mfd(4:6, :, :) - mfd_diff 

   call trans_tdf_calc(kg, ene, mfd_drag, tdf_tmp, .false.)

   tdf_tmp = tdf_tmp * merge(1.0_dp, 2.0_dp, ncol) / (vol * tmpr)
   do i = 1, ncomp
      tdf_tmp(i,:) = tdf_tmp(i,:) * mfermi_deriv(tmpr, (ene-ef))
      cond_seebeck(i) = -sum(tdf_tmp(i,:)) * estep
   enddo

   !note that ph_alpha_T = alpha / T
   !     so alpha * T^-1 = ph_alpha_T
   do i = 1, ncomp
      lambda(i) = ph_alpha_T(i) / cond_seebeck(i)
     !if(ionode) write(*,'(1x,a,i2,a,f12.6,a,i2,a,f12.6,a,i2,a,f12.6)') 'cond_seebeck(',i,') = ', cond_seebeck(i) * seebeck_unit * cond_unit, ' ph_alpha_T(',i,') = ', ph_alpha_T(i)* seebeck_unit * cond_unit,' lambda(',i,') = ', lambda(i)
   enddo

  !if(ionode) write(*,*) 'lambda', lambda

   mfd(4:6, :, :) = mfd_diff + lambda(1) * mfd_drag

   deallocate(tdf_tmp, mfd_diff, mfd_drag, cond_seebeck)
end subroutine

!subroutine trans_impose_kelveinonsager(kg, ene, tmpr, ef, ncol, vol, ph_alpha_T, mfd)
!   use pert_const, only: dp, ryd2ev, bohr2ang, kelvin2ev, ryd2mev, unitcharge,&
!                         timeunit, tesla2ryd, pi
!   implicit none
!   type(grid), intent(in) :: kg
!   real(dp), intent(in) :: tmpr, ef, ene(:), vol
!   logical,  intent(in) :: ncol
!   real(dp), intent(in) :: ph_alpha_T(:)
!   real(dp), intent(inout) :: mfd(:,:,:)
!   ! local variables
!   integer :: i, nestep, ncomp, ik, ib
!   real(dp) :: estep
!   real(dp) :: lambda(6)
!   real(dp), allocatable :: cond_seebeck(:)
!   real(dp), allocatable :: tdf_tmp(:,:)
!
!   real(dp), allocatable :: mfd_diff(:,:,:)
!   real(dp), allocatable :: mfd_drag(:,:,:)
!
!   real(dp) :: seebeck_unit, cond_unit
!
!   seebeck_unit = kelvin2ev * 1.0E3_dp
!   cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)
!
!   nestep = size(ene)
!
!   !ncomp = 6: Electric and T fields
!   ncomp = size(mfd, 1) 
!
!   if(ionode) write(*,*) 'ph_alpha_T', ph_alpha_T
!   if(ncomp .ne. size(ph_alpha_T)) &
!      call errore('trans_impose_kelveinonsager', 'mismatch dimension in input arguments!',1)
!   allocate( tdf_tmp(ncomp, nestep) )
!   allocate( mfd_diff(ncomp/2, kg%numb, kg%nk) )
!   allocate( mfd_drag(ncomp/2, kg%numb, kg%nk) )
!   allocate( cond_seebeck(ncomp) )
!   estep = (ene(nestep)-ene(1)) / (nestep-1.0_dp)
!   
!   call trans_tdf_calc(kg, ene, mfd(4:6,:,:), tdf_tmp, .false.)
!   tdf_tmp = tdf_tmp * merge(1.0_dp, 2.0_dp, ncol) / (vol * tmpr)
!   do i = 1, ncomp
!      tdf_tmp(i,:) = tdf_tmp(i,:) * mfermi_deriv(tmpr, (ene-ef))
!      cond_seebeck(i) = -sum(tdf_tmp(i,:)) * estep
!   enddo
!   do i = 1, ncomp
!      if(ionode) write(*,'(1x,a,i2,a,f12.6,a,i2,a,f12.6,a,i2,a,f12.6)') '1-cond_seebeck(',i,') = ', cond_seebeck(i) * seebeck_unit * cond_unit, ' ph_alpha_T(',i,') = ', ph_alpha_T(i)* seebeck_unit * cond_unit,' lambda(',i,') = ', lambda(i)
!   enddo
!
!
!
!   do ik = 1, kg%nk
!      do ib = 1, kg%numb
!         do i = 1, ncomp/2
!            mfd_diff(i, ib, ik) = mfd(i, ib, ik) * (kg%enk(ib,ik) - ef)
!         enddo
!      enddo
!   enddo
!
!   call trans_tdf_calc(kg, ene, mfd_diff, tdf_tmp, .false.)
!   tdf_tmp = tdf_tmp * merge(1.0_dp, 2.0_dp, ncol) / (vol * tmpr)
!   do i = 1, ncomp
!      tdf_tmp(i,:) = tdf_tmp(i,:) * mfermi_deriv(tmpr, (ene-ef))
!      cond_seebeck(i) = -sum(tdf_tmp(i,:)) * estep
!   enddo
!   do i = 1, ncomp
!      if(ionode) write(*,'(1x,a,i2,a,f12.6,a,i2,a,f12.6,a,i2,a,f12.6)') '2-cond_seebeck(',i,') = ', cond_seebeck(i) * seebeck_unit * cond_unit, ' ph_alpha_T(',i,') = ', ph_alpha_T(i)* seebeck_unit * cond_unit,' lambda(',i,') = ', lambda(i)
!   enddo
!
!   mfd_drag = mfd(4:6, :, :) - mfd_diff 
!
!
!
!   call trans_tdf_calc(kg, ene, mfd_drag, tdf_tmp, .false.)
!   tdf_tmp = tdf_tmp * merge(1.0_dp, 2.0_dp, ncol) / (vol * tmpr)
!   do i = 1, ncomp
!      tdf_tmp(i,:) = tdf_tmp(i,:) * mfermi_deriv(tmpr, (ene-ef))
!      cond_seebeck(i) = -sum(tdf_tmp(i,:)) * estep
!   enddo
!
!   !note that ph_alpha_T = alpha / T
!   !     so alpha * T^-1 = ph_alpha_T
!   do i = 1, ncomp
!      lambda(i) = ph_alpha_T(i) / cond_seebeck(i)
!      if(ionode) write(*,'(1x,a,i2,a,f12.6,a,i2,a,f12.6,a,i2,a,f12.6)') 'cond_seebeck(',i,') = ', cond_seebeck(i) * seebeck_unit * cond_unit, ' ph_alpha_T(',i,') = ', ph_alpha_T(i)* seebeck_unit * cond_unit,' lambda(',i,') = ', lambda(i)
!   enddo
!
!   if(ionode) write(*,*) 'lambda', lambda
!
!   mfd(4:6, :, :) = mfd_diff + lambda(1) * mfd_drag
!
!
!   call trans_tdf_calc(kg, ene, mfd(4:6,:,:), tdf_tmp, .false.)
!   tdf_tmp = tdf_tmp * merge(1.0_dp, 2.0_dp, ncol) / (vol * tmpr)
!   do i = 1, ncomp
!      tdf_tmp(i,:) = tdf_tmp(i,:) * mfermi_deriv(tmpr, (ene-ef))
!      cond_seebeck(i) = -sum(tdf_tmp(i,:)) * estep
!   enddo
!   do i = 1, ncomp
!      if(ionode) write(*,'(1x,a,i2,a,f12.6,a,i2,a,f12.6,a,i2,a,f12.6)') '4-cond_seebeck(',i,') = ', cond_seebeck(i) * seebeck_unit * cond_unit, ' ph_alpha_T(',i,') = ', ph_alpha_T(i)* seebeck_unit * cond_unit,' lambda(',i,') = ', lambda(i)
!   enddo
!
!   deallocate(tdf_tmp, mfd_diff, mfd_drag, cond_seebeck)
!end subroutine


!calculate conductivity and mobility using tetrahedron integration. 
subroutine trans_cond_calc(kg, ene, tmpr, ef, ncol, vol, mfd, cond, tdf, calc_mag, calc_ph_thermal)
   implicit none
   type(grid), intent(in) :: kg
   logical,  intent(in) :: ncol, calc_mag
   real(dp), intent(in) :: tmpr, ef, ene(:), mfd(:,:,:), vol
   real(dp), intent(out) :: cond(:)
   real(dp), intent(out), optional :: tdf(:,:)
   logical,  intent(in), optional :: calc_ph_thermal
   ! local variables
   integer :: i, nestep, ncomp
   logical :: calc_ph_thermal_local
   real(dp) :: estep
   real(dp), allocatable :: tdf_tmp(:,:)

   call stopwatch('trans_cond_calc','start')

   if( present(calc_ph_thermal) ) then
      calc_ph_thermal_local = calc_ph_thermal
   else
      calc_ph_thermal_local = .false.
   endif

   nestep = size(ene)

   !ncomp here is not the same as ncomp in transport.f90
   !ncomp = 6 : Just electric field
   !ncomp = 9 : Both electric and magnetic fields
   !ncomp = 12: Electric and T fields
   ncomp = size(mfd, 1) * merge(3,2,calc_mag)  
   allocate( tdf_tmp(ncomp, nestep) )
   estep = (ene(nestep)-ene(1)) / (nestep-1.0_dp)
   ! compute TDF_ij(E) = spinor/(N*V) sum_{nk} vnk*vnk*tau*delta(E-e_nk)
   ! the unit of velocity is Ryd*Bohr*alat, vnk = 1/hbar*dE/dk.
   ! NOTE: tdf computed below does not include spinor and volume V !!!
   if (.not. calc_ph_thermal_local) then
      ! for electrons, using tetrahedra integration
      call trans_tdf_calc(kg, ene, mfd, tdf_tmp, calc_mag)
   else
      ! for phonons, using uniform grid integration
      call trans_tdf_calc(kg, ene, mfd, tdf_tmp, calc_mag)
   endif

   !account for spin-degeneracy and volume
   if (calc_ph_thermal_local) then
      ! syp @04252024 - for ph-BTE we do not need to account for spin degeneracy
      !               - also we have a temperature in the denominator
      tdf_tmp = tdf_tmp / (tmpr * vol)
   else
      tdf_tmp = tdf_tmp * merge(1.0_dp, 2.0_dp, ncol) / vol
   endif
   !
   if( present(tdf) ) then
      if( size(tdf, 1) .ne. ncomp .or. size(tdf, 2) .ne. nestep) &
         call errore('boltz_trans_mod',"argument 'tdf': mismatch dimension.", 1)
      !N.B.: tdf is in Rydberg atomic unit (omitted a factor of (alat)^2).
      tdf(:,:) = tdf_tmp(:,:)
   endif
   !
   !do the integration over energy to compute conductivity.
   do i = 1, ncomp
      !sigma = \sum_{nk} F_nk * vnk * -df/dE
      if (calc_ph_thermal_local) then
         ! for phonons, using bose-einstein distribution
         ! for phonon-thermal conductivity, there are two (hbar omega):
         !     one is absorbed into tdf already, the other is here.
         tdf_tmp(i,:) = tdf_tmp(i,:) * mbose_deriv(tmpr, ene) * ene
      else
         ! for electrons, using fermi-dirac distribution
         tdf_tmp(i,:) = tdf_tmp(i,:) * mfermi_deriv(tmpr, (ene-ef))
      endif

      !sigma = q^2*nspin/(N*V) sum_{nk} v_nk*v_nk*tau_nk* f *(1-f)/kT
      ! 1. q^2 is omitted in the calculation of cond(i), 
      ! 2. all the quantities are in Rydberg atomic unit.
      ! 3. (alat)^2 should be applied to cond(i) to get the actual sigma
      cond(i) = sum(tdf_tmp(i,:)) * estep
   enddo
   deallocate(tdf_tmp)
   !
   call stopwatch('trans_cond_calc','stop')
end subroutine trans_cond_calc

! compute seebeck and thermal conductivity
subroutine extract_trans_coeff(tmpr, cond, cond_seebeck, kk_coef, sys_2d, seebeck, t_cond, ncomp2, seebeck_e, alpha, seebeck_ph, alpha_ph)
   implicit none
   !Size is either 6 or 9 depending on whether
   !magnetic field is present
   integer,  intent(in) :: ncomp2
   real(dp), intent(in) :: tmpr, cond(:), cond_seebeck(:), kk_coef(:)
   logical, intent(in) :: sys_2d
   real(dp), intent(out) :: seebeck(ncomp2), t_cond(ncomp2)
   real(dp), intent(out), optional :: seebeck_e(ncomp2), seebeck_ph(ncomp2)
   real(dp), intent(in), optional :: alpha(ncomp2), alpha_ph(ncomp2)
   ! local variables
   real(dp) :: cs_tmp(3,3), tcs_tmp(3,3), c_tmp(3,3), s_tmp(3,3), c_aux(3,3)
   real(dp) :: s_tmp2(3,3), tcs_tmp2(3,3), s_tmp3(3,3)
   
   ! convert rank-2 tensor to matrix
   call tensor2matrix(cond, c_tmp) ! sigma
   call tensor2matrix(cond_seebeck, cs_tmp) ! sigma*seebeck

   !therm_cond = K - T * [sigma*seebeck] * seebeck
   ! if alpha is present, then replace "T * [sigma*seebeck]" with alpha
   if(present(alpha)) then
      call tensor2matrix(alpha, tcs_tmp) ! alpha_e
   else
      if (drag) then
         call errore('extract_trans_coeff', 'alpha need to be provided for coupled BTE calculation!',1)
      else
         tcs_tmp = tmpr * cs_tmp ! alpha = T * [sigma*seebeck]
      endif
   endif
   if(present(alpha_ph)) then
      call tensor2matrix(alpha_ph, tcs_tmp2) ! alpha_ph
   else
      tcs_tmp2 = 0 
   endif

   c_aux = 0.0_dp
   s_tmp = 0.0_dp
   if( sys_2d ) then
      c_aux(1:2, 1:2) = mat_inv2(c_tmp(1:2, 1:2)) ! sigma^-1
      ! seebeck = sigma^-1 * [sigma*seebeck]
      s_tmp(1:2, 1:2) = matmul(c_aux(1:2, 1:2), cs_tmp(1:2, 1:2)) ! sigma^-1 * [sigma*seebeck] = seebeck
      
      ! electron peltier thermalpower = alpah (sigma T)^-1 = alpah * sigma^-1 /T
      s_tmp2(1:2, 1:2) = matmul(c_aux(1:2, 1:2), tcs_tmp(1:2, 1:2))/tmpr

      ! phonon peltier thermalpower = alpah (sigma T)^-1 = alpah * sigma^-1 /T
      s_tmp3(1:2, 1:2) = matmul(c_aux(1:2, 1:2), tcs_tmp2(1:2, 1:2))/tmpr

      ! T * [sigma*seebeck] * seebeck or alpha * seebeck
      c_aux(1:2, 1:2) = matmul(tcs_tmp(1:2, 1:2), s_tmp(1:2, 1:2)) ! for Kappa = K - xxx
   else
      c_aux = mat_inv3( c_tmp ) ! sigma^-1
      ! seebeck = sigma^-1 * [sigma*seebeck]
      s_tmp = matmul(c_aux, cs_tmp) ! sigma^-1 * [sigma*seebeck] = seebeck

      ! electron peltier thermalpower = alpah (sigma T)^-1 = alpah * sigma^-1 /T
     !s_tmp2 = matmul(tcs_tmp, c_aux)/tmpr 
      s_tmp2 = matmul(c_aux, tcs_tmp)/tmpr 

      ! phonon peltier thermalpower = alpah (sigma T)^-1 = alpah * sigma^-1 /T
      s_tmp3 = matmul(c_aux, tcs_tmp2)/tmpr 

      ! T * [sigma*seebeck] * seebeck or alpha * seebeck
      c_aux = matmul(tcs_tmp, s_tmp) ! for Kappa = K - xxx
   endif
   !back to tensor
   call matrix2tensor(s_tmp, seebeck, ncomp2)
   !back to tensor
   if (present(seebeck_e)) call matrix2tensor(s_tmp2, seebeck_e, ncomp2)
   !back to tensor
   if (present(seebeck_ph)) call matrix2tensor(s_tmp3, seebeck_ph, ncomp2)
   !therm_cond = K - T * [sigma*seebeck] * seebeck
   call matrix2tensor(c_aux, t_cond, ncomp2)
   t_cond(1:ncomp2) = kk_coef(1:ncomp2) - t_cond(1:ncomp2)
   !
end subroutine 


! convert rank-2 tensor to matrix
subroutine tensor2matrix(tensor, matrix)
   implicit none
   real(dp), intent(in) :: tensor(:)
   real(dp), intent(out) :: matrix(3,3)
   !
   integer :: i, j, n

   matrix = 0.0_dp
   n = 0
   do j = 1, 3
   do i = 1, j
      !(1,1), (1,2), (2,2), (1,3), (2,3), (3,3)
      n = n + 1
      matrix(i,j) = tensor(n)
      !
      if((i .ne. j) .and. (size(tensor) .eq. 6)) matrix(j,i) = tensor(n)
   enddo; enddo

   if(size(tensor) .eq. 9) then
      matrix(2,1) = tensor(7)
      matrix(3,1) = tensor(8)
      matrix(3,2) = tensor(9)
   endif
end subroutine


! convert 2D matrix to rank-2 tensor
subroutine matrix2tensor(matrix, tensor, ncomp)
   implicit none
   integer,  intent(in) :: ncomp !Contains either 6 or 9 elements
   real(dp), intent(in) :: matrix(3,3)
   real(dp), intent(out) :: tensor(ncomp)
   !
   integer :: i, j, n
   real(dp) :: rtmp

   tensor = 0.0_dp
   n = 0
   do j = 1, 3
   do i = 1, j
      !(1,1), (1,2), (2,2), (1,3), (2,3), (3,3)
      n = n + 1
      if(ncomp .eq. 6) then
         tensor(n) = (matrix(i,j) + matrix(j,i)) * 0.5_dp
      else
         tensor(n) = matrix(i,j)
      !sanity check
      !rtmp = ( abs(matrix(i,j)) + abs(matrix(j,i)) ) * 0.5_dp
      !if( (i.ne.j) .and. rtmp>1.0E-16_dp .and. abs(matrix(i,j)-matrix(j,i)) > 1.0E-6_dp*rtmp ) &
      !   call errore('matrix2tensor','only symmetric matrix can be converted to tensor.',1)
      endif
   enddo; enddo

   if(ncomp .eq. 9) then
      tensor(7) = matrix(2,1) 
      tensor(8) = matrix(3,1) 
      tensor(9) = matrix(3,2)
   endif
end subroutine

end module boltz_trans_mod
