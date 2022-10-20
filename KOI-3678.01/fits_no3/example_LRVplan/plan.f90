! =====================================================================================
! PLAN version 1.1
! v1.1 - Minor updates for speed enhancements
!
! AUTHOR: David Kipping
!         Harvard-Smithsonian Center for Astrophysics
!         Email: dkipping@cfa.harvard.edy
!
! =====================================================================================

MODULE planmod
use mandelmod
implicit none

CONTAINS

! ==============================================================
! === SUBROUTINE: PLAN ===
!
SUBROUTINE plan(t,Pb,T_b,p,ab,eb,wb,bb,u1,u2,fplan,nz,IpL)

implicit none

 INTEGER :: i, j, nz
 REAL(8) :: tstep
 REAL(8), INTENT(IN) :: Pb, T_b, p, ab, eb, wb, bb
 REAL(8), INTENT(IN) :: u1, u2, fplan
 REAL(8) :: rb_mid, cosib, temp
 REAL(8), DIMENSION(nz), INTENT(IN) :: t
 REAL(8), DIMENSION(nz) :: tb, fb, varrhob
 REAL(8), DIMENSION(nz) :: S_P
 REAL(8), DIMENSION(nz), INTENT(OUT) :: IpL       ! Limb darkened output
 REAL(8), DIMENSION(nz) :: IpN       ! Non-LDed output
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0

 ! Obtain periastrion planet separation and inclination
 rb_mid = ab*(1.0D0 - eb**2)/(1.0D0 + eb*DSIN(wb))
 cosib = bb/rb_mid     ! We use cosine here because ib=90 for bb=0

 ! Offset times
 DO i=1,nz
  tb(i) = t(i) - T_b
 END DO

 ! Kepler converts time array in f array
 call kepler(eb,wb,Pb,fplan,nz,tb,fb)

 ! Obtain instantaneous planet separation
 temp = 1.0D0 - eb**2
 DO i=1,nz
  varrhob(i) = temp/(1.0D0 + eb*DCOS(fb(i)))
 END DO

 ! Planetary motion
 DO i=1,nz
  ! Now S_P
  S_P(i) = ( ab*varrhob(i)*DCOS(fb(i)+wb) )**2 &
           + ( ab*cosib*varrhob(i)*DSIN(fb(i)+wb) )**2
  S_P(i) = DSQRT(S_P(i))
 END DO

 ! Get LD's flux...
 ! IpL = LDed flux; IpN = non-LDed flux
 call occultquad(S_P,u1,u2,p,IpL,IpN,nz)

 END SUBROUTINE plan

! ==============================================================
! === SUBROUTINE: KEPLER ===
!
 SUBROUTINE kepler(e,wrad,Pdays,f_ref,n,t,f_)
 ! Solves Kepler's equation and calculates f(t)

 implicit none
 INTEGER :: i, j, n!, tt
 REAL(8), INTENT(IN) :: e, wrad, Pdays, f_ref
 REAL(8) :: E_ref, M_ref, ek, Pk, toler
 INTEGER :: ok
 REAL(8), DIMENSION(n), INTENT(IN) :: t
 REAL(8), DIMENSION(n) :: M_, E_, E_p
 REAL(8), DIMENSION(n), INTENT(OUT) :: f_
 REAL(8), PARAMETER :: halfpi = 1.5707963267948966D0
 REAL(8), PARAMETER :: twopi = 6.2831853071795864769D0
     
 ! Time is in days, so all references to P must be Pdays
 ! The following values are predefined for speed
 ek = DSQRT((1.0D0+e)/(1.0D0-e))
 Pk = twopi/Pdays
 IF( e .LT. 0.9D-4 ) THEN
   E_ref = f_ref
   M_ref = E_ref
   ! 2) Solve Kepler's equation
   DO i=1,n
     M_(i) = Pk*t(i) + M_ref
     E_(i) = M_(i)
     f_(i) = E_(i)
   END DO
 ELSE
   ! If e>0, we must first solve f(t).
   E_ref = 2.0D0*DATAN((1.0D0/ek)*DTAN(f_ref*0.5D0))
   IF(E_ref .LT. -halfpi) THEN
     E_ref = E_ref + twopi
   END IF
   M_ref = E_ref - e*DSIN(E_ref)
   ! 2) Solve Kepler's equation
   toler = 1.1574074074074074D-4*Pk ! Set to 10 second tolerance, divide by 10 for 1 second
   DO i=1,n
     M_(i) = Pk*t(i) + M_ref
   END DO
   E_(:) = M_(:)
   DO i=1,n
     ok = 0
     DO WHILE ( ok == 0 )
       E_p(i) = E_(i) - ((E_(i) - e*DSIN(E_(i)) &
                - M_(i))/(1.0D0-e*DCOS(E_(i))))
       IF( DABS(E_p(i)-E_(i)) .LT. toler ) THEN
         ok = 1
       END IF
       E_(i) = E_p(i)
     END DO
     f_(i) = 2.0D0*DATAN(ek*DTAN(E_(i)*0.5D0))
   END DO
 END IF

 END SUBROUTINE kepler
! =======================================================

END MODULE planmod
