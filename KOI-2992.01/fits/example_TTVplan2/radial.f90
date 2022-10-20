! RADIAL version 3.0
! For use with JAMMI v3.0
! Handles 2 planets + gradient + curl
! Numerous re-writes to generalize code
! Adaptive Kepler's eqn
! Input parameter come from R now, rather than individual statements
!
! AUTHOR: David Kipping
!         Harvard-Smithsonian Center for Astrophysics/University College London
!         Email: dkipping@cfa.harvard.edu
!

MODULE radialmod
use params

 implicit none

 contains

 subroutine radial(Rin,chi2)

 ! IMPLICITS
 REAL(8), DIMENSION(sdim), INTENT(IN) :: Rin                 ! R parameter vector
 REAL(8) :: ttran, P1, K1, esinw1, ecosw1                   ! Planet 1 inputs
 REAL(8) :: rhomstar, bp, rhoP, aR, cosi_1, k_1, h_1
 REAL(8) :: tconj2, P2, K2, esinw2, ecosw2                  ! Planet 2 inputs
 REAL(8) :: gradsec, grad, curl, ttroj, jitter              ! System inputs
 REAL(8) :: fref1, fref2, etas
 REAL(8), INTENT(OUT) :: chi2
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8) :: e1, w1, tref1
 REAL(8) :: e2, w2, tref2
 REAL(8), DIMENSION(nrlen) :: t1, t2, vmod, trP
 INTEGER :: i
 ! Dense output parameters
 REAL(8), DIMENSION(:), ALLOCATABLE :: tr_dense, t1_dense, t2_dense, trP_dense
 REAL(8), DIMENSION(:), ALLOCATABLE :: vobs_dense, vsig_dense, vmod_dense
 INTEGER, DIMENSION(:), ALLOCATABLE :: rvmode_dense
 REAL(8) :: maxtime, mintime, steptime
 INTEGER :: n_dense

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! SECTION 1: Read-in and assign initial conditions/priors
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! Planet 1
 P1 = Rin(4)
 ttran = Rin(5) !pivotRV + Rin(5)*0.15915494309189535D0*P1
 K1 = 0.0D0 !Rin(9)
 e1 = 2.0D-8 !Rin(7)
 w1 = 0.7853981633974483D0 !Rin(8)

 ! Planet 1 inclination
 IF( priflag .EQ. 1 .AND. e1 .GT. 2.0D-8 ) THEN
   rhomstar = DSQRT(Rin(2)**3)  ! rho_{*}
   bp = Rin(3)        ! Barycentre's impact parameter
   k_1 = e1*DCOS(w1)
   h_1 = e1*DSIN(w1)
   rhoP = (1.0D0-e1**2)/(1.0D0+h_1)
   aR = rhomstar*Grvx*P1**2
   aR = aR**third
   cosi_1 = bp/(aR*rhoP)
 END IF

 ! Planet 2
 P2 = 1.0D0 !Rin(9)
 tconj2 = pivotRV !+ Rin(10)*0.15915494309189535D0*P2 ! = pivot point + (P2*phi2)/(2pi)
 K2 = 0.0D0 !Rin(11)
 e2 = 2.0D-8
 w2 = 0.7853981633974483D0

 ! System
 gradsec = 0.0D0 !Rin(7)*0.0027379092633269355D0
 grad = 0.0D0 !Rin(8)*0.0027379092633269355D0
 curl = 0.0D0 !Rin(3)
 ttroj = 0.0D0
 jitter = 0.0D0 !Rin(8)

 ! fref1, simply chose as the time of inferior conjunction of planet 1
 fref1 = 1.5707963267948966D0 - w1
 IF( fref1 .LT. 0.0D0 ) THEN
  fref1 = fref1 + 6.283185307179586D0
 END IF
 IF( priflag .EQ. 1 .AND. e1 .GT. 2.0D-8 ) THEN
   etas = 0.0D0
   etas = etas + ((k_1)/(1.0D0+h_1))*cosi_1**2				!eta1
   etas = etas + ((k_1)/(1.0D0+h_1))*(1.0D0/(1.0D0+h_1))*cosi_1**4	!eta2
   fref1 = fref1 - etas
 END IF

 ! fref2, simply chose as the time of inferior conjunction of planet 2
 fref2 = 1.5707963267948966D0 - w2
 IF( fref2 .LT. 0.0D0 ) THEN
  fref2 = fref2 + 6.283185307179586D0
 END IF

 ! Calculate reference points in time
 ! tref1 = time of conjunction of 1. The eta terms account for the slight offset
 ! between ttran and tconj.
 tref1 = ttran !+ ttroj
 tref2 = tconj2

 !! Set the the radial time stamp to be that from rammi (protects JAMMI from feedback)
 DO i=1,nrlen
  t1(i) = tr(i)
  t2(i) = tr(i)
  ! Pivot point (used for gradient/curl stuff)
  trP(i) = tr(i) - pivotRV
 END DO

 ! Subtract tref's
 DO i=1,nrlen
  t1(i) = t1(i) - tref1
  t2(i) = t2(i) - tref2
 END DO

 ! Call model generator
 call rasmine(jitter,grad,gradsec,curl,&
              P1,K1,e1,w1,fref1,&
              P2,K2,e2,w2,fref2,&
              rvmode,Vsyslen,&
              t1,t2,tr,trP,rv,sigrv,nrlen,vmod)

 ! Calculate chi^2
 chi2 = 0.0D0
 DO i=1,nrlen
  chi2 = chi2 + ( (rv(i)-vmod(i))**2 / (sigrv(i)**2+jitter**2) )
 END DO
 !!write(*,*) 'Chi^2 = ',chi2

 end subroutine radial

! ======================================================
 subroutine rasmine(jitter,grad,gradsec,curl,&
                    P1,K1,e1,w1,f_ref1,&
                    P2,K2,e2,w2,f_ref2,&
                    rvmode,Vsyslen,&
                    t1,t2,tr,trP,vobs,vsig,n,v)

 ! Rasmine is a radial velocity model generator for a two-planet system.
 ! The name Rasmine comes from its twinning to Jasmine,
 ! which is used to generate transit models (at least in the old days)
 !
 ! Rasmine is called by Radial,
 ! just like how Jasmine is called by Transit
 !
 ! Radial is ultimately called by Rammi,
 ! just like how Transit gets called by Jammi

 REAL(8), DIMENSION(Vsyslen) :: Vsys, Vsyswei
 REAL(8), INTENT(IN) :: jitter, grad, gradsec, curl 
 REAL(8), INTENT(IN) :: K1, P1, e1, w1, f_ref1
 REAL(8), INTENT(IN) :: K2, P2, e2, w2, f_ref2
 INTEGER, INTENT(IN) :: Vsyslen
 REAL(8), DIMENSION(n), INTENT(IN) :: t1, t2, tr, trP, vobs, vsig
 INTEGER, DIMENSION(n), INTENT(IN) :: rvmode
 REAL(8), DIMENSION(n) :: f1, f2, v1, v2, vx
 REAL(8), DIMENSION(n), INTENT(OUT) :: v
 INTEGER :: n, i, j, ok, m
 REAL(8) :: toler
 REAL(8) :: E_ref1, M_ref1, M_i1, E_i1, E_ip1   ! Planet 1 Kepler Eqn bits
 REAL(8) :: E_ref2, M_ref2, M_i2, E_i2, E_ip2   ! Planet 2 Kepler Eqn bits
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8), PARAMETER :: Grv = 6.67428D-11
 REAL(8), PARAMETER :: third = 0.333333333333333D0

 !!! PLANET 1 !!!
 IF( e1 .LT. 0.9D-4 ) THEN
  ! Calculate v1 array
  DO i=1,n
   v1(i) = -K1*DSIN( (2.0D0*pi*t1(i))/P1 )
  END DO
 ELSE
  ! Compute reference mean anomaly, based upon reference true anomaly
  E_ref1 = 2.0D0*DATAN(DSQRT((1.0D0-e1)/(1.0D0+e1))*DTAN(0.5D0*f_ref1))
  IF(E_ref1 < -0.5D0*pi) THEN
   E_ref1 = E_ref1 + 2.0D0*pi
  END IF
  M_ref1 = E_ref1 - e1*DSIN(E_ref1)
  ! Convert t1 -> f1
  toler = (pi*1.0D0*0.000023148148148148147)/(P1) ! 1.0D0 = 1 second tolerance
  DO i=1,n
   ! Begin iteration
   M_i1 = ((2.0D0*pi*t1(i))/P1) + M_ref1 ! Mean anomaly (M_i1 = M_i for planet 1)
   E_i1 = M_i1                           ! Eccentric anomaly first guess
   ok = 0 ! Logic flag which turns on when iteration is accurate enough
   DO WHILE ( ok == 0 )
    E_ip1 = E_i1 - ( (E_i1 - e1*DSIN(E_i1) - M_i1) / (1.0D0-e1*DCOS(E_i1)) )
    IF( DABS(E_ip1-E_i1) .LT. toler ) THEN
     ok = 1
    END IF
    E_i1 = E_ip1
   END DO
   f1(i) = 2.0D0*DATAN(DSQRT((1.0D0+e1)/(1.0D0-e1))*DTAN(0.5D0*E_i1))
  END DO
  ! Calculate v1 array
  DO i=1,n
   v1(i) = DCOS(f1(i)+w1) + e1*DCOS(w1)
   v1(i) = v1(i)*K1
  END DO
 END IF

 !!! PLANET 2 !!!
 IF( e2 .LT. 0.9D-4 ) THEN
  ! Calculate v2 array
  DO i=1,n
   v2(i) = -K2*DSIN( (2.0D0*pi*t2(i))/P2 )
  END DO
 ELSE
  ! Compute reference mean anomaly, based upon reference true anomaly
  E_ref2 = 2.0D0*DATAN(DSQRT((1.0D0-e2)/(1.0D0+e2))*DTAN(f_ref2/2.0D0))
  IF(E_ref2 < -0.5D0*pi) THEN
   E_ref2 = E_ref2 + 2.0D0*pi
  END IF
  M_ref2 = E_ref2 - e2*DSIN(E_ref2)
  ! Convert t2 -> f2
  toler = (pi*1.0D0*0.000023148148148148147)/(P2) ! 1.0D0 = 1 second tolerance
  DO i=1,n
   ! Begin iteration
   M_i2 = ((2.0D0*pi*t2(i))/P2) + M_ref2 ! Mean anomaly (M_i2 = M_i for planet 2)
   E_i2 = M_i2                           ! Eccentric anomaly first guess
   ok = 0 ! Logic flag which turns on when iteration is accurate enough
   DO WHILE ( ok == 0 )
    E_ip2 = E_i2 - ( (E_i2 - e2*DSIN(E_i2) - M_i2) / (1.0D0-e2*DCOS(E_i2)) )
    IF( DABS(E_ip2-E_i2) .LT. toler ) THEN
     ok = 1
    END IF
    E_i2 = E_ip2
   END DO
   f2(i) = 2.0D0*DATAN(DSQRT((1.0D0+e2)/(1.0D0-e2))*DTAN(0.5D0*E_i2))
  END DO
  ! Calculate v2 array
  DO i=1,n
   v2(i) = DCOS(f2(i)+w2) + e2*DCOS(w2)
   v2(i) = v2(i)*K2
  END DO 
 END IF

 !!! Planet X !!!
 ! (Planet X is really just the grad and curl components)
 DO i=1,n
    vx(i) = trP(i)*grad + trP(i)*trP(i)*curl
 END DO

 ! Combine TWO planets with offsets...
 DO i=1,n
   v(i) = vx(i) + v1(i) + v2(i)
 END DO

 ! Linearly optimized Vsys
 Vsys(:) = 0.0D0
 Vsyswei(:) = 0.0D0
 DO j=1,Vsyslen
   DO i=1,n
     IF( rvmode(i) .EQ. j ) THEN
       Vsyswei(j) = Vsyswei(j) + (1.0D0/(vsig(i)**2+jitter**2))
       Vsys(j) = Vsys(j) + (vobs(i)-v(i))/(vsig(i)**2+jitter**2)
     END IF
   END DO
   IF( Vsyswei(j) .NE. 0.0D0 ) THEN
     Vsys(j) = Vsys(j)/Vsyswei(j)
   ELSE
     Vsys(j) = 1984.0D0 ! error flag
   END IF
 END DO

 ! Offset RVs
 DO i=1,n
   v(i) = v(i) + Vsys(rvmode(i))
 END DO
 
 end subroutine rasmine
! =======================================================

END MODULE radialmod
