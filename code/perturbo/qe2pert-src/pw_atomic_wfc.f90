!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  adapted from PW/src/atomic_wfc.f90
!-----------------------------------------------------------------------
SUBROUTINE pw_atomic_wfc(xk, npw, igk_k, wfcatom )
  !-----------------------------------------------------------------------
  !! This routine computes the superposition of atomic wavefunctions
  !! for k-point "ik" - output in "wfcatom".
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi, fpi, pi
  USE cell_base,        ONLY : omega, tpiba
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,            ONLY : natomwfc
  USE gvect,            ONLY : mill, eigts1, eigts2, eigts3, g
  !USE klist,            ONLY : xk, igk_k, ngk
  USE wvfct,            ONLY : npwx
  USE uspp_data,        ONLY : tab_at, dq
  USE uspp_param,       ONLY : upf
  USE noncollin_module, ONLY : noncolin, domag, npol, angle1, angle2, &
                               starting_spin_angle
  USE upf_spinorb,      ONLY : rot_ylm, lmaxx
  USE mp_bands,         ONLY : inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  real(dp), intent(in) :: xk(3)  ! xk in cartesian coordinate
  INTEGER, INTENT(IN) :: npw
  INTEGER, INTENT(IN) :: igk_k(npw)
  !! k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom( npwx, npol, natomwfc )
  !! Superposition of atomic wavefunctions
  !
  ! ... local variables
  !
  INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, nwfcm
  REAL(DP),    ALLOCATABLE :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
  COMPLEX(DP), ALLOCATABLE :: sk (:), aux(:)
  COMPLEX(DP) :: kphase, lphase
  REAL(DP)    :: arg, px, ux, vx, wx
  INTEGER     :: ig_start, ig_end

  !CALL start_clock( 'atomic_wfc' )

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  DO nt = 1, ntyp
     lmax_wfc = MAX( lmax_wfc, MAXVAL( upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  END DO
  !
  nwfcm = MAXVAL( upf(1:ntyp)%nwfc )
  !npw = ngk(ik)
  !
  ALLOCATE( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  DO ig = 1, npw
     iig = igk_k (ig)
     gk (1,ig) = xk(1) + g(1,iig)
     gk (2,ig) = xk(2) + g(2,iig)
     gk (3,ig) = xk(3) + g(3,iig)
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  END DO
  !
  !  ylm = spherical harmonics
  !
  CALL ylmr2( (lmax_wfc+1)**2, npw, gk, qg, ylm )

  ! from now to the end of the routine the ig loops are distributed across bgrp
  CALL divide( inter_bgrp_comm,npw,ig_start,ig_end )
  !
  ! set now q=|k+G| in atomic units
  !
  DO ig = ig_start, ig_end
     qg(ig) = SQRT( qg(ig) )*tpiba
  END DO
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  DO nt = 1, ntyp
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc (nb) >= 0.d0 ) THEN
           DO ig = ig_start, ig_end
              px = qg (ig) / dq - INT(qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           END DO
        END IF
     END DO
  END DO

  DEALLOCATE( qg, gk )
  ALLOCATE( aux(npw) )
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  !
  DO na = 1, nat
     arg = (xk(1)*tau(1,na) + xk(2)*tau(2,na) + xk(3)*tau(3,na)) * tpi
     kphase = CMPLX( COS(arg), - SIN(arg) ,KIND=DP)
     !
     !     sk is the structure factor
     !
     DO ig = ig_start, ig_end
        iig = igk_k (ig)
        sk (ig) = kphase * eigts1 (mill (1,iig), na) * &
                           eigts2 (mill (2,iig), na) * &
                           eigts3 (mill (3,iig), na)
     END DO
     !
     nt = ityp (na)
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
           l = upf(nt)%lchi(nb)
           lphase = (0.d0,1.d0)**l
           !
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
           !
           IF ( noncolin ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !
                 IF (starting_spin_angle.OR..NOT.domag) THEN
                    CALL atomic_wfc_so( )
                 ELSE
                    CALL atomic_wfc_so_mag( )
                 END IF
                 !
              ELSE
                 !
                 CALL atomic_wfc_nc( )
                 !
              END IF
              !
           ELSE
              !
              CALL atomic_wfc___( )
              !
           END IF
           !
        END IF
        !
     END DO
     !
  END DO

  IF ( n_starting_wfc /= natomwfc) call errore ('atomic_wfc', &
       'internal error: some wfcs were lost ', 1 )

  DEALLOCATE( aux, sk, chiq, ylm )

  ! collect results across bgrp
  CALL mp_sum( wfcatom, inter_bgrp_comm )

  !CALL stop_clock( 'atomic_wfc' )
  RETURN

CONTAINS
!----------------------------------------------------------------
  SUBROUTINE atomic_wfc_so( )
   !------------------------------------------------------------
   !! Spin-orbit case.
   !
   REAL(DP) :: fact(2), j
   REAL(DP), EXTERNAL :: spinor
   INTEGER :: ind, ind1, n1, is, sph_ind
   !
   j = upf(nt)%jchi(nb)
   DO m = -l-1, l
      fact(1) = spinor(l,j,m,1)
      fact(2) = spinor(l,j,m,2)
      IF ( ABS(fact(1)) > 1.d-8 .OR. ABS(fact(2)) > 1.d-8 ) THEN
         n_starting_wfc = n_starting_wfc + 1
         IF (n_starting_wfc > natomwfc) CALL errore &
              ('atomic_wfc_so', 'internal error: too many wfcs', 1)
         DO is=1,2
            IF (abs(fact(is)) > 1.d-8) THEN
               ind=lmaxx+1+sph_ind(l,j,m,is)
               aux=(0.d0,0.d0)
               DO n1=1,2*l+1
                  ind1=l**2+n1
                  if (abs(rot_ylm(ind,n1)) > 1.d-8) &
                      aux(:)=aux(:)+rot_ylm(ind,n1)*ylm(:,ind1)
               ENDDO
               do ig = ig_start, ig_end
                  wfcatom(ig,is,n_starting_wfc) = lphase*fact(is)*&
                        sk(ig)*aux(ig)*chiq (ig, nb, nt)
               END DO
            ELSE
                wfcatom(:,is,n_starting_wfc) = (0.d0,0.d0)
            END IF
         END DO
      END IF
   END DO
   !
   END SUBROUTINE atomic_wfc_so
   ! 
   SUBROUTINE atomic_wfc_so_mag( )
   !
   !! Spin-orbit case, magnetization along "angle1" and "angle2"
   !! In the magnetic case we always assume that magnetism is much larger
   !! than spin-orbit and average the wavefunctions at l+1/2 and l-1/2
   !! filling then the up and down spinors with the average wavefunctions,
   !! according to the direction of the magnetization, following what is
   !! done in the noncollinear case.
   !
   REAL(DP) :: alpha, gamman, j
   COMPLEX(DP) :: fup, fdown  
   REAL(DP), ALLOCATABLE :: chiaux(:)
   INTEGER :: nc, ib
   !
   j = upf(nt)%jchi(nb)
   !
   !  This routine creates two functions only in the case j=l+1/2 or exit in the
   !  other case 
   !    
   IF (ABS(j-l+0.5_DP)<1.d-4) RETURN

   ALLOCATE(chiaux(npw))
   !
   !  Find the functions j=l-1/2
   !
   IF (l == 0)  THEN
      chiaux(:)=chiq(:,nb,nt)
   ELSE
      DO ib=1, upf(nt)%nwfc
         IF ((upf(nt)%lchi(ib) == l).AND. &
                      (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) THEN
            nc=ib
            EXIT
         ENDIF
      ENDDO
      !
      !  Average the two functions
      !
      chiaux(:)=(chiq(:,nb,nt)*(l+1.0_DP)+chiq(:,nc,nt)*l)/(2.0_DP*l+1.0_DP)
      !
   ENDIF 
   !
   !  and construct the starting wavefunctions as in the noncollinear case.
   !
   alpha = angle1(nt)
   gamman = - angle2(nt) + 0.5d0*pi
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc + 2*l+1 > natomwfc ) CALL errore &
            ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
      DO ig = ig_start, ig_end
         aux(ig) = sk(ig)*ylm(ig,lm)*chiaux(ig)
      END DO
      !
      ! now, rotate wfc as needed
      ! first : rotation with angle alpha around (OX)
      !
      DO ig = ig_start, ig_end
         fup = cos(0.5d0*alpha)*aux(ig)
         fdown = (0.d0,1.d0)*sin(0.5d0*alpha)*aux(ig)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         wfcatom(ig,1,n_starting_wfc) = (cos(0.5d0*gamman) &
                        +(0.d0,1.d0)*sin(0.5d0*gamman))*fup
         wfcatom(ig,2,n_starting_wfc) = (cos(0.5d0*gamman) &
                        -(0.d0,1.d0)*sin(0.5d0*gamman))*fdown
         !
         ! second: rotation with angle gamma around (OZ)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         fup = cos(0.5d0*(alpha+pi))*aux(ig)
         fdown = (0.d0,1.d0)*sin(0.5d0*(alpha+pi))*aux(ig)
         !
         ! second, rotation with angle gamma around (OZ)
         !
         wfcatom(ig,1,n_starting_wfc+2*l+1) = (cos(0.5d0*gamman) &
                  +(0.d0,1.d0)*sin(0.5d0 *gamman))*fup
         wfcatom(ig,2,n_starting_wfc+2*l+1) = (cos(0.5d0*gamman) &
                  -(0.d0,1.d0)*sin(0.5d0*gamman))*fdown
      END DO
   END DO
   n_starting_wfc = n_starting_wfc + 2*l+1
   DEALLOCATE( chiaux )
   !
   END SUBROUTINE atomic_wfc_so_mag
   !
   SUBROUTINE atomic_wfc_nc( )
   !
   !! noncolinear case, magnetization along "angle1" and "angle2"
   !
   REAL(DP) :: alpha, gamman
   COMPLEX(DP) :: fup, fdown  
   !
   alpha = angle1(nt)
   gamman = - angle2(nt) + 0.5d0*pi
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc + 2*l+1 > natomwfc) CALL errore &
            ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
      DO ig = ig_start, ig_end
         aux(ig) = sk(ig)*ylm(ig,lm)*chiq(ig,nb,nt)
      END DO
      !
      ! now, rotate wfc as needed
      ! first : rotation with angle alpha around (OX)
      !
      DO ig = ig_start, ig_end
         fup = cos(0.5d0*alpha)*aux(ig)
         fdown = (0.d0,1.d0)*sin(0.5d0*alpha)*aux(ig)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         wfcatom(ig,1,n_starting_wfc) = (cos(0.5d0*gamman) &
                        +(0.d0,1.d0)*sin(0.5d0*gamman))*fup
         wfcatom(ig,2,n_starting_wfc) = (cos(0.5d0*gamman) &
                        -(0.d0,1.d0)*sin(0.5d0*gamman))*fdown
         !
         ! second: rotation with angle gamma around (OZ)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         fup = cos(0.5d0*(alpha+pi))*aux(ig)
         fdown = (0.d0,1.d0)*sin(0.5d0*(alpha+pi))*aux(ig)
         !
         ! second, rotation with angle gamma around (OZ)
         !
         wfcatom(ig,1,n_starting_wfc+2*l+1) = (cos(0.5d0*gamman) &
                  +(0.d0,1.d0)*sin(0.5d0 *gamman))*fup
         wfcatom(ig,2,n_starting_wfc+2*l+1) = (cos(0.5d0*gamman) &
                  -(0.d0,1.d0)*sin(0.5d0*gamman))*fdown
      END DO
   END DO
   n_starting_wfc = n_starting_wfc + 2*l+1
   !
   END SUBROUTINE atomic_wfc_nc

   SUBROUTINE atomic_wfc___( )
   !
   ! ... LSDA or nonmagnetic case
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc > natomwfc) CALL errore &
         ('atomic_wfc___', 'internal error: too many wfcs', 1)
      !
      DO ig = ig_start, ig_end
         wfcatom (ig, 1, n_starting_wfc) = lphase * &
            sk (ig) * ylm (ig, lm) * chiq (ig, nb, nt)
      ENDDO
      !
   END DO
   !
   END SUBROUTINE atomic_wfc___
   !
END SUBROUTINE pw_atomic_wfc

!----------------------------------------------------------
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Set of subroutines needed for full LDA+U calculations
! after Liechtenstein and co-workers (PRB 52, R5467 (1995)). 
! Works with two-component spinor WFs and with fully-relativistic
! pseudopotentials. 
! In the last case the WFs are projected onto: 
!   real spherical harmonics * 
!   averaged j=l+1/2, l-1/2 radial WFs *
!   up/down spinor.
! 
! A. Smogunov, C. Barreteau
!
! adapted from PW/src/plus_u_full.f90/atomic_wfc_nc_updown
!
SUBROUTINE pw_atomic_wfc_nc_updown (xk, npw, igk_k, wfcatom)
  !-----------------------------------------------------------------------
  !
  ! For noncollinear case: builds up the superposition (for a k-point "ik") of 
  ! pure spin up or spin down atomic wavefunctions.
  ! 
  ! Based on atomic_wfc.f90

  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi, fpi, pi
  USE cell_base,  ONLY : tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : mill, eigts1, eigts2, eigts3, g
  !USE klist,      ONLY : xk, ngk, igk_k
  USE wvfct,      ONLY : npwx !, nbnd
  USE uspp_data,  ONLY : tab_at, dq
  USE uspp_param, ONLY : upf
  USE noncollin_module, ONLY : noncolin, domag, npol, angle1, angle2, &
                               starting_spin_angle
  USE upf_spinorb,   ONLY : rot_ylm, lmaxx
  !
  implicit none
  !
  real(dp), intent(in) :: xk(3) ! xk in cartesian coordinate
  integer, intent(in) :: npw
  integer, intent(in) :: igk_k(npw)
  complex(DP), intent(out) :: wfcatom (npwx, npol, natomwfc)
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, nwfcm
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
  complex(DP), allocatable :: sk (:), aux(:)
  complex(DP) :: kphase
  real(DP) :: arg, px, ux, vx, wx

  !call start_clock ('atomic_wfc')

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  do nt = 1, ntyp
     lmax_wfc = MAX ( lmax_wfc, MAXVAL (upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  enddo
  !
  nwfcm = MAXVAL ( upf(1:ntyp)%nwfc )
  !npw = ngk(ik)
  allocate ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  do ig = 1, npw
     iig = igk_k(ig)
     gk (1,ig) = xk(1) + g(1,iig)
     gk (2,ig) = xk(2) + g(2,iig)
     gk (3,ig) = xk(3) + g(3,iig)
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = 1, npw
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  do nt = 1, ntyp
     do nb = 1, upf(nt)%nwfc
        if ( upf(nt)%oc (nb) >= 0.d0) then
           do ig = 1, npw
              px = qg (ig) / dq - int (qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           enddo
        endif
     enddo
  enddo

  deallocate (qg, gk)
  allocate ( aux(npw) )
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  !
  do na = 1, nat
     arg = (xk(1)*tau(1,na) + xk(2)*tau(2,na) + xk(3)*tau(3,na)) * tpi
     kphase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     !
     !     sk is the structure factor
     !
     do ig = 1, npw
        iig = igk_k(ig)
        sk (ig) = kphase * eigts1 (mill (1,iig), na) * &
                           eigts2 (mill (2,iig), na) * &
                           eigts3 (mill (3,iig), na)
     enddo
     !
     nt = ityp (na)
     do nb = 1, upf(nt)%nwfc
        if (upf(nt)%oc(nb) >= 0.d0) then
           l = upf(nt)%lchi(nb)
           !
           !
              IF ( upf(nt)%has_so ) THEN
                 !
                 call wfc_atom ( .true. )
                 !
              ELSE
                 !
                 call wfc_atom ( .false. )
                 !
              ENDIF
              !
        END IF
        !
     END DO
     !
  END DO

  if (n_starting_wfc /= natomwfc) call errore ('atomic_wfc_nc_updown', &
       'internal error: some wfcs were lost ', 1)

  deallocate(aux, sk, chiq, ylm)

  call stop_clock ('atomic_wfc')
  return

CONTAINS

   SUBROUTINE wfc_atom ( soc )
   !
   !
   real(DP) :: j
   real(DP), ALLOCATABLE :: chiaux(:)
   integer :: nc, ib
   logical :: soc ! .true. if the fully-relativistic pseudo
   !

!  If SOC go on only if j=l+1/2
   if (soc) j = upf(nt)%jchi(nb)
   if (soc.and.ABS(j-l+0.5_DP)<1.d-4 ) return
!

   allocate (chiaux(npw))

   if (soc) then 

!
!  Find the index for j=l-1/2
!
     if (l == 0)  then
        chiaux(:)=chiq(:,nb,nt)
     else
        do ib=1, upf(nt)%nwfc
           if ((upf(nt)%lchi(ib) == l).and. &
                        (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) then
              nc=ib
              exit
           endif
        enddo
!
!  Average the two radial functions 
!
        chiaux(:)=(chiq(:,nb,nt)*(l+1.0_DP)+chiq(:,nc,nt)*l)/(2.0_DP*l+1.0_DP)
     endif

   else

     chiaux(:) = chiq(:,nb,nt)

   endif


   do m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1

      if (n_starting_wfc + 2*l+1 > natomwfc) call errore &
            ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
      do ig=1,npw
         aux(ig) = sk(ig)*ylm(ig,lm)*chiaux(ig)
      enddo
!
      do ig=1,npw
!
         wfcatom(ig,1,n_starting_wfc) = aux(ig) 
         wfcatom(ig,2,n_starting_wfc) = 0.d0 
!
         wfcatom(ig,1,n_starting_wfc+2*l+1) = 0.d0
         wfcatom(ig,2,n_starting_wfc+2*l+1) = aux(ig) 
!
      enddo
   enddo
   n_starting_wfc = n_starting_wfc + 2*l+1

   deallocate (chiaux)
   !
   END SUBROUTINE wfc_atom
   !
END SUBROUTINE pw_atomic_wfc_nc_updown
