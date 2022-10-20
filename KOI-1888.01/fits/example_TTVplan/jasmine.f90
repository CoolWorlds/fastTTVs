MODULE jasminemod

! ===================
! Jasmine Version 2.2
! ===================
! David Kipping
! dkipping@ucl.ac.uk
!
! CITATION:
! If using this code, please cite:
! Kipping, D., 2009, MNRAS, 389, 1383
!
! VERSION:
! Version 2.0 is a complete re-write of Jasmine
! Version 2.1 uses RV fingerprint soln in all cases
! Version 2.2 uses a series expansion for fpri&fsec
!  which enables total reliability up to high precision.
!
! CHANGES:
! * Jasmine now works with Cos[f] angles to minimize number of computations
! * Jasmine also is now is more diffuse; split up into multiple subroutines
!   which perform various functions. This is done so make the code more readable
!   & fewer lines (~3 times fewer despite more comments) & easily checked for errors.
! * Jasmine now solves for the mid-transit and mid-eclipse anomalies (new feature)
! * Jasmine also solves for the RV minimum and finds the time offset (new feature)
!
! KNOWN ISSUES:
! * w -> Multiples of close to 90 degrees are very fishy; the quartic tends to yield 
!        solutions in extreme proximity (below the double precision in many cases).
!        Consequently, Jasmine has 'wsens' = w_sensitivity parameter which is fixed
!        in the main code. If w ~ 90 degrees, w is offset by wsens to try to prevent
!        these errors. wsens is chosen to be as small as possible, but obviously
!        increasing it avoids the risk of the w90 bug occuring. The biggest consequence
!        of these problems are that the advanced mode parameters are suspect.
!        A resolution is to map out the fingerprints over the whole parameter space
!        to obviate the whole process of the solution selection. This is in the works.
!
! BASIC PROCESS OF HOW JASMINE FIND THE TRANSIT DURATIONS:
! 1. Calculates coefficients of the quartic equation, which solves for Cos[f]
!    [ subroutine = durcoeff ]
! 2. Solves the quartic to give four candidate solutions for Cos[f]
!    [ subroutine = quartic ]
! 3. Converts the four Cos[f] solutions into eight f solutions. 
!    [ in-house ]
! 4. Each pair solutions contains one physical and one unphysical solution. 
!    These are distinguished by computing S for each solution.
!    [ Scheck ]
! 5. Distinguishes which solutions occur for primary and which for the secondary
!    [ in-house ]
! 6. Calculates guess solutions. This is usually switched off but is useful for
!    debugging.
! 7. Feeds final f values into subroutine to calculate the duration
!    [ subroutine = duration ]
! 8. Solution "fingerprint" - the solution can be charactized by a series of
!    logic controls. Usually switched off.
!
! SYMBOLS USED (DIMENSIONS GIVEN IN PARANTHESES):
!
! DurComF(4) = Complex candidate Cos[f] solutions
! DurCosF(4) = Real candidate Cos[f] solutions
! DurF(4,2) = Real candidate f solutions
! S2_DurF(4,2) = S^2 evaluated at the candidate f solutions
! DurChoice(2) = Logical array: 1=1st Dur solution accepted; 2=2nd
! S2_Acc(2) = S^2 evaluated at the accepted Cos[f] solutions
! PriF(k,2) = Primary f solutions for kth contact scenario
! SecF(k,2) = Secondary f solutions for kth contact scenario
! T(k,2) = Primary duration for kth contact scenario
! S(k,2) = Secondary duration for kth contact scenario
!
! Note - arrays are designed to have fastest moving components on the left

implicit none
CONTAINS
SUBROUTINE jasmine(bP,p,e,win,aR,Pdays,durexact,&
           tT,tF,DurP,T12,T34,rhostar,ideg,secpri,sT,RVoffset,fpri)

 ! Controls
 INTEGER :: i, j, y, k, vlogfile, verbose
 INTEGER :: calcpri, calcsec, calcmid, calcapprox, calcDFMT
 INTEGER, INTENT(IN) :: durexact
 INTEGER :: calcIEA, calcRV, calcSECOFF, durfinger, midfinger
 ! Physicals
 REAL(8), INTENT(IN) :: bP, p, e, win, aR, Pdays
 REAL(8) :: w, Pp
 ! Substitutions
 REAL(8) :: lam2, csc2i, bS, wsens
 REAL(8), DIMENSION(3) :: v          ! Contact scenarios
 REAL(8) :: rhoP, rhoS, irad
 ! Guess parameters
 REAL(8) :: fPguess, fSguess, delfP, delfS
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 2.0 DURATION
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 REAL(8), DIMENSION(5) :: Dcoeff     ! Coefficients for the duration quartic
 COMPLEX(8), DIMENSION(4) :: DurComF ! Complex candidate Cos[f] solutions
 REAL(8), DIMENSION(4) :: DurCosF    ! Real candidate Cos[f] solutions
 REAL(8), DIMENSION(4,2) :: DurF     ! Candidate f solutions
 REAL(8), DIMENSION(4,2) :: S2_DurF  ! S^2 evaluated at the candidate f solutions
 INTEGER, DIMENSION(4) :: DurChoice  ! Logical array: 1=1st Dur solution accepted; 2=2nd
 INTEGER :: prifound, secfound       ! Counters for number of pri/sec solutions
 INTEGER, DIMENSION(4) :: priflag    ! Logical array: 1=primary eclipse; 0 = secondary eclipse
 REAL(8), DIMENSION(2,3) :: PriF, SecF       ! Primary/Secondary f solutions
 REAL(8), DIMENSION(3) :: T, S         ! Final durations
 REAL(8), INTENT(OUT):: tT, DurP, tF, sT  ! Final durations
 REAL(8) :: DurS, sF
 REAL(8), DIMENSION(3) :: T_approx, S_approx ! Final durations
 REAL(8), DIMENSION(2,2) :: pricombo, seccombo     ! Different combinations to make the ingress/egress
 REAL(8) :: angwidthP, angwidthS                   ! Full-half-angular widths of tT & sT
 REAL(8), DIMENSION(2) :: slopedurP, slopeangP     ! Possible combinations for ingress/egress
 REAL(8), DIMENSION(2) :: slopedurS, slopeangS     ! Possible combinations for ingress/egress
 REAL(8), INTENT(OUT) :: T12, T34        ! Final ingress/egress durations
 REAL(8) :: S12, S34        ! Final ingress/egress durations
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 3.0 MID -POINT
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 REAL(8), INTENT(OUT) :: fpri ! Mid-transit angles
 REAL(8) :: fsec              ! Mid-transit angles
 REAL(8) :: kk, hh, hhP1, hhM1, prefac, errT, errF, cos2i
 REAL(8) :: eta                ! eta is the offset from fPguess
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 4.0 ADVANCED
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 REAL(8) :: halfT, DFMT, IEA 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 5.0 RV TIME OFFSET
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 REAL(8) :: fRVmin
 REAL(8), INTENT(OUT) :: RVoffset             ! Final values
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 6.0 REFORMATTING
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 REAL(8), INTENT(OUT):: rhostar, ideg, secpri   ! Some extra outputs
 REAL(8) :: secpriapprox
 ! Constants
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8), PARAMETER :: Grv = 6.67428D-11

 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! SECTION 1.0                                             !
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! Definitions                                             !
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 ! Open the log file
 vlogfile = 0    ! Write output to verbose log file?
 IF( vlogfile == 1 ) THEN
 open(11,FILE='jasmine.vlog',FORM='FORMATTED',STATUS='UNKNOWN')
 END IF

 ! Controls
 verbose = 0     ! Output durations to screen?
 calcpri = 1     ! Calculate primary durations?
 calcsec = 1     ! Calculate secondary durations?
 calcmid = 1     ! Calculate mid-anomalies?
 calcapprox = 1  ! Calculate approx formulae?
 durfinger = 0   ! Print the duration fingerprint?
 midfinger = 0   ! Print the mid fingerprint?
 calcDFMT = 0    ! Calculate Delta-Flux-Minima-Timing?
 calcIEA = 0     ! Calculate Ingress-Egress-Asymmetry?
 calcRV = 1      ! Calculate RV time offset?
 calcSECOFF = 1  ! Calculate pri-sec offset time?
 ! Contacts to use
 v(1) = (1.0D0 + p)**2        ! Contact scenario 1
 v(2) = 1.0D0                 ! Contact scenario 2
 v(3) = (1.0D0 - p)**2        ! Contact scenario 3
 ! P correction
 Pp = Pdays*86400.0D0

 IF( durfinger == 1 ) THEN
  open(21,FILE='durfinger.dat',FORM='FORMATTED',ACCESS='APPEND')
 END IF
 IF( midfinger == 1 ) THEN
  open(22,FILE='midfinger.dat',FORM='FORMATTED',ACCESS='APPEND')
 END IF

 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 IF( vlogfile == 1 ) THEN
  write(11,*) '================ JASMINE ================'
 END IF
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 ! Physical Values
 IF ( aR*(1.0D0-e) < 1.0D0 ) THEN
  !write(*,*) 'Orbit is inside star!'
 END IF

 ! w corrections
 w = win
 IF( w .GE. 2.0D0*pi ) THEN
  w = w - 2.0D0*pi
 END IF

 ! Substitutions
 lam2 = (aR*(1.0D0-e*e))**2
 rhoP = (1.0D0-e*e)/(1.0D0+e*DSIN(w))
 rhoS = (1.0D0-e*e)/(1.0D0-e*DSIN(w))
 bS = (rhoS/rhoP)*bP
 !IF( bS > 1.0D0 ) THEN
 ! IF( verbose == 2 ) THEN
 !  write(*,*) 'Secondary eclipse does not occur, b_S = ',bS
 ! END IF
 !END IF
 csc2i = 1.0D0 - (bP/(aR*rhoP))**2
 csc2i = 1.0D0/csc2i

 ! Assign inclination
 irad = DASIN(DSQRT(1.0D0/csc2i))

 ! fP guess
 fPguess = 0.5D0*pi - w
 IF( fPguess < 0.0D0 ) THEN
  fPguess = fPguess + 2.0D0*pi
 END IF
 ! fS guess
 fSguess = fPguess - pi
 IF( fSguess < 0.0D0 ) THEN
  fSguess = fSguess + 2.0D0*pi
 END IF

!==================================================================!
!                                                                  !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! SECTION 2                                               !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! Duration calculations                                   !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!                                                                  !
!==================================================================!

 IF( durexact == 1 ) THEN

 DO k=1,3

  ! If verbose, delcare what calculation you are doing first
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) '======================================='
   IF( k == 1 ) THEN
    write(11,*) '=========Calculating T_{1,4}==========='
   ELSE IF( k == 2 ) THEN
    write(11,*) '=======Calculating T_{1.5,3.5}========='
   ELSE IF( k == 3 ) THEN
    write(11,*) '=========Calculating T_{2,3}==========='
   END IF
   write(11,*) '======================================='
  END IF

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.1: Call coefficients for quartic
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  call durcoeffs(lam2,e,w,csc2i,v(k),Dcoeff)

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.2: Call quartic eqn
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  call quartic(Dcoeff(1),Dcoeff(2),Dcoeff(3),Dcoeff(4),Dcoeff(5),DurComF(1),DurComF(2),DurComF(3),DurComF(4))
  ! Assign real parts
  DO i = 1,4
   DurCosF(i) = DurComF(i)
  END DO
  
  ! Write solutions
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) 'Candidate Cos[f] solutions are...'
   write(11,*) '------------------------------------------------------------'
   DO i=1,4
    write(11,*) 'Cos[f] = ',DurCosF(i)
   END DO
   write(11,*) '------------------------------------------------------------'
  END IF

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.3: Eight candidate f values
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  DO i=1,4
   DurF(i,1) = DACOS(DurCosF(i))
   DurF(i,2) = -DACOS(DurCosF(i)) + 2.0D0*pi
  END DO

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.4: SubStage A of solution selection
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! SubStage A: Find physical/unphysical solutions using Scheck
  call Scheck(DurF,w,e,csc2i,lam2,v(k),4,S2_DurF,DurChoice)

  ! Write candidate solutions
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) 'Candidate f solutions...'
   write(11,*) '------------------------------------------------------------'
   DO i=1,4
    write(11,*) 'PAIR ',i
    write(11,*) 'f1 = ',DurF(i,1)*(180.0D0/pi),' deg --> Delta(S^2) = ',ABS(S2_DurF(i,1)-v(k))
    write(11,*) 'f2 = ',DurF(i,2)*(180.0D0/pi),' deg --> Delta(S^2) = ',ABS(S2_DurF(i,2)-v(k))
    write(11,*) 'Choice made = ',DurChoice(i)
   END DO
   write(11,*) '------------------------------------------------------------'
  END IF

  ! Write accepted solutions
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) 'Accepted f solutions occur at...'
   write(11,*) '------------------------------------------------------------'
   DO i=1,4
    write(11,*) 'f = ',DurF(i,DurChoice(i))*(180.0/pi),' deg --> Delta(S^2) = ',&
    ABS(S2_DurF(i,DurChoice(i)) - v(k))
    END DO
   write(11,*) '------------------------------------------------------------'
   write(11,*) ' '
  END IF

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.5: SubStage B of solution selection
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! SubStage B: Distinguish solutions
  prifound = 0
  secfound = 0
  DO i=1,4
   ! If Z > 0 then we have a primary solution
   IF( DSIN(w+DurF(i,DurChoice(i)))/(1.0D0+e*DCOS(DurF(i,DurChoice(i)))) > 0.0D0 ) THEN
    prifound = prifound + 1
    PriF(prifound,k) = DurF(i,DurChoice(i))
    priflag(i) = 1  ! Logic flag to mark as a primary eclipse
   ELSE
    secfound = secfound + 1
    SecF(secfound,k) = DurF(i,DurChoice(i))
    priflag(i) = 0  ! Logic flag to mark as a secondary eclipse
   END IF
  END DO

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.6: Guess solutions
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  IF( calcapprox == 1 ) THEN
   ! Delta(f) guess [defined as full width of transit]
   delfP = (v(k) - bP*bP)/(aR*rhoP*aR*rhoP - bP*bP)
   delfP = 2.0D0*DASIN(DSQRT(delfP))
   IF( bP*bP > v(k) ) THEN
    delfP = 0.0D0
   END IF
   ! Secondary eclipse
   delfS = (v(k) - bS*bS)/(aR*rhoS*aR*rhoS - bS*bS)
   delfS = 2.0D0*DASIN(DSQRT(delfS))
   IF( bS*bS > v(k) ) THEN
    delfS = 0.0D0
   END IF
  END IF

  ! Final angles
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) 'Primary solutions...'
   write(11,*) '------------------------------------------------------------'
   write(11,*) PriF(1,k)*(180.0D0/pi),' and ',PriF(2,k)*(180.0D0/pi),' degrees'
   IF( calcapprox == 1 ) THEN
    write(11,*) '~',(fPguess-0.5D0*delfP)*(180.0D0/pi),&
    ' and ~',(fPguess+0.5D0*delfP)*(180.0D0/pi),' degrees'
   END IF
   write(11,*) ' '
   write(11,*) 'Secondary solutions...'
   write(11,*) '------------------------------------------------------------'
   write(11,*) SecF(1,k)*(180.0D0/pi),' and ',SecF(2,k)*(180.0D0/pi),' degrees'
   IF( calcapprox == 1 ) THEN
    write(11,*) '~',(fSguess-0.5D0*delfS)*(180.0D0/pi),&
    ' and ~',(fSguess+0.5D0*delfS)*(180.0D0/pi),' degrees'
   END IF
   write(11,*) ' '
  END IF
  
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.7: Calculate duration
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Primary transit
  IF( calcpri == 1 ) THEN
   call duration(e,Pp,PriF(1,k),PriF(2,k),T(k))
   ! Non-transiting override
   IF( bp > DSQRT(v(k)) ) THEN
    IF( vlogfile == 1 ) THEN
     write(11,*) 'OVERRIDE: bp > cutoff => T = 0'
    END IF
    T(k) = 0.0D0
   END IF
  END IF

  ! Secondary eclipse
  IF( calcsec == 1 ) THEN
   call duration(e,Pp,SecF(1,k),SecF(2,k),S(k))
   ! Non-transiting override
   IF( bs > DSQRT(v(k)) ) THEN
    IF( vlogfile == 1 ) THEN
     write(11,*) 'OVERRIDE: bs > cutoff => S = 0'
    END IF
    S(k) = 0.0D0
   END IF
  END IF

  ! Output the primary duration
  IF( vlogfile == 1 ) THEN
   IF( calcpri == 1 ) THEN
    write(11,*) ' '
    IF( k == 1 ) THEN
     write(11,*) 'Duration T_{1,4}'
    ELSE IF( k == 2 ) THEN
     write(11,*) 'Duration T_{1.5,3.5}'
    ELSE IF( k == 3 ) THEN
     write(11,*) 'Duration T_{2,3}'
    END IF
    write(11,*) '------------------------------------------------------------'
    write(11,*) 'duration = ',T(k),' seconds'
   END IF
  END IF
  IF( verbose == 1 ) THEN
   IF( calcpri == 1 ) THEN
    write(*,*) ' '
    IF( k == 1 ) THEN
     write(*,*) 'Duration T_{1,4}'
    ELSE IF( k == 2 ) THEN
     write(*,*) 'Duration T_{1.5,3.5}'
    ELSE IF( k == 3 ) THEN
     write(*,*) 'Duration T_{2,3}'
    END IF
    write(*,*) 'duration = ',T(k),' seconds'
   END IF
  END IF
  ! Approx primary duration formula
  IF( calcapprox == 1 ) THEN
   IF( calcpri == 1 ) THEN
    call durapprox(e,Pp,rhoP,bP,aR,v(k),T_approx(k))
    IF( vlogfile == 1 ) THEN
     write(11,*) 'duration ~ ',T_approx(k),' seconds'
     write(11,*) '------------------------------------------------------------'
    END IF
   END IF
  END IF

  ! Output the secondary duration
  IF( vlogfile == 1 ) THEN
   IF( calcsec == 1 ) THEN
    write(11,*) ' '
    IF( k == 1 ) THEN
     write(11,*) 'Duration S_{1,4}'
    ELSE IF( k == 2 ) THEN
     write(11,*) 'Duration S_{1.5,3.5}'
    ELSE IF( k == 3 ) THEN
     write(11,*) 'Duration S_{2,3}'
    END IF
    write(11,*) '------------------------------------------------------------'
    write(11,*) 'duration = ',S(k),' seconds'
   END IF
  END IF
  IF( verbose == 1 ) THEN
   IF( calcsec == 1 ) THEN
    write(*,*) ' '
    IF( k == 1 ) THEN
     write(*,*) 'Duration S_{1,4}'
    ELSE IF( k == 2 ) THEN
     write(*,*) 'Duration S_{1.5,3.5}'
    ELSE IF( k == 3 ) THEN
     write(*,*) 'Duration S_{2,3}'
    END IF
    write(*,*) 'duration = ',S(k),' seconds'
   END IF
  END IF
  ! Approx secondary duration formula
  IF( calcapprox == 1 ) THEN
   IF( calcsec == 1 ) THEN
    call durapprox(e,Pp,rhoS,bS,aR,v(k),S_approx(k))
    IF( vlogfile == 1 ) THEN
     write(11,*) 'duration ~ ',S_approx(k),' seconds'
     write(11,*) '------------------------------------------------------------'
    END IF
   END IF
  END IF

 END DO  ! This is the end of the contact scenario loop

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.8: Obtain T_{1,2} and T_{3,4}
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 ! SubSubSection 2.8.1: Primary transits

 IF( vlogfile == 1 .AND. calcpri == 1 ) THEN
  write(11,*) ' '
  write(11,*) '======================================='
  write(11,*) '=========Calculating T_{1,2}==========='
  write(11,*) '======================================='
  write(11,*) ' '
 END IF
 ! Can only be done if we have a full transit
 IF( bp < DSQRT(v(3)) .AND. calcpri .EQ. 1 ) THEN

  ! Angular width of tT
  angwidthP = ABS(PriF(1,1) - PriF(2,1))
  IF( angwidthP > pi ) THEN
   angwidthP = 0.5D0*(2.0D0*pi - angwidthP)
  ELSE
   angwidthP = 0.5D0*angwidthP
  END IF
  IF( vlogfile == 1 ) THEN
   write(11,*) 'Full-half-angular-width = ',angwidthP*(180.0D0/pi),' deg'
  END IF

  ! Four combinations could be the ingress/egress
  DO j=1,2
   DO y=1,2
    ! Define angular width for each combo
    pricombo(j,y) = ABS(PriF(j,3) - PriF(y,1))
    ! Test for a flip
    IF( pricombo(j,y) > pi ) THEN
     pricombo(j,y) = 2.0D0*pi - pricombo(j,y)
    END IF
    IF( vlogfile == 1 ) THEN
     write(11,*) ' '
     write(11,*) 'Slope trial angles = ',PriF(j,3)*(180.0D0/pi),' & ',PriF(y,1)*(180.0D0/pi),' degs'
     write(11,*) 'Angular width = ',pricombo(j,y)*(180.0D0/pi),' deg'
    END IF
    ! Now test if less than the tT width or not
    IF( pricombo(j,y) < angwidthP ) THEN
     ! ACCEPTED...
     IF( vlogfile == 1 ) THEN
      write(11,*) 'Accepted'
     END IF
     ! Get duration
     call duration(e,Pp,PriF(j,3),PriF(y,1),slopedurP(j))
     slopeangP(j) = ABS(PriF(j,3) - PriF(y,1))
     ! Test for a flip
     IF( slopeangP(j) > pi ) THEN
      slopeangP(j) = 2.0D0*pi - slopeangP(j)
     END IF
     IF( ABS(PriF(j,3)-PriF(y,1)) < pi ) THEN
      ! Normal case, assign slopeang as the addition on top of the lowest angle
      IF( PriF(j,3) < PriF(y,1) ) THEN
       slopeangP(j) = PriF(j,3) + slopeangP(j)
      ELSE
       slopeangP(j) = PriF(y,1) + slopeangP(j)
      END IF
     ELSE
      ! Flipped case, one of the angles <360 and one is >0
      IF( PriF(j,3) < PriF(y,1) ) THEN
       slopeangP(j) = PriF(y,1) + slopeangP(j)
      ELSE
       slopeangP(j) = PriF(j,3) + slopeangP(j)
      END IF
      IF( slopeangP(j) > 2.0D0*pi ) THEN
       slopeangP(j) = slopeangP(j) - 2.0D0*pi
      END IF
     END IF
    ELSE
     ! REJECTED...
     IF( vlogfile == 1 ) THEN
      write(11,*) 'Not accepted'
     END IF
    END IF
   END DO ! y loop ends
   ! Output so far...
   IF( vlogfile == 1 ) THEN
    write(11,*) ' '
    write(11,*) 'slopedur(',j,') = ',slopedurP(j),' with midangle = ',slopeangP(j)*(180.D0/pi), 'deg'
   END IF
  END DO ! jloop
  ! Assign t12 and t34
  ! Whatever angle is lowest is assumed to be t12
  IF( ABS(slopeangP(2) - slopeangP(1)) < pi ) THEN
   ! Normal case
   IF( slopeangP(1) < slopeangP(2) ) THEN
    T12 = slopedurP(1); T34 = slopedurP(2)
   ELSE
    T12 = slopedurP(2); T34 = slopedurP(1)
   END IF
  ELSE
   ! Flipped case
   IF( slopeangP(1) < slopeangP(2) ) THEN
    T12 = slopedurP(2); T34 = slopedurP(1)
   ELSE
    T12 = slopedurP(1); T34 = slopedurP(2)
   END IF
  END IF
 ELSE
  T12 = 0.0D0; T34 = 0.0D0
 END IF ! Full transit IF statement

 ! Output primary results
 IF( calcpri == 1 .AND. vlogfile == 1 ) THEN
  write(11,*) ' '
  write(11,*) 'T12 = ',T12,' seconds'
  write(11,*) 'T34 = ',T34,' seconds'
 END IF
 IF( calcpri == 1 .AND. verbose == 1 ) THEN
  write(*,*) ' '
  write(*,*) 'T12 = ',T12,' seconds'
  write(*,*) 'T34 = ',T34,' seconds'
 END IF

! ! SubSubSection 2.8.2: Secondary eclipses

 IF( vlogfile == 1 .AND. calcsec == 1 ) THEN
  write(11,*) ' '
  write(11,*) '======================================='
  write(11,*) '=========Calculating S_{1,2}==========='
  write(11,*) '======================================='
  write(11,*) ' '
 END IF
 ! Can only be done if we have a full transit
 IF( bs < DSQRT(v(3)) .AND. calcsec .EQ. 1 ) THEN

  ! Angular width of tT
  angwidthS = ABS(SecF(1,1) - SecF(2,1))
  IF( angwidthS > pi ) THEN
   angwidthS = 0.5D0*(2.0D0*pi - angwidthS)
  ELSE
   angwidthS = 0.5D0*angwidthS
  END IF
  IF( vlogfile == 1 ) THEN
   write(11,*) 'Full-half-angular-width = ',angwidthS*(180.0D0/pi),' deg'
  END IF

  ! Four combinations could be the ingress/egress
  DO j=1,2
   DO y=1,2
    ! Define angular width for each combo
    seccombo(j,y) = ABS(SecF(j,3) - SecF(y,1))
    ! Test for a flip
    IF( seccombo(j,y) > pi ) THEN
     seccombo(j,y) = 2.0D0*pi - seccombo(j,y)
    END IF
    IF( vlogfile == 1 ) THEN
     write(11,*) ' '
     write(11,*) 'Slope trial angles = ',SecF(j,3)*(180.0D0/pi),' & ',SecF(y,1)*(180.0D0/pi),' degs'
     write(11,*) 'Angular width = ',seccombo(j,y)*(180.0D0/pi),' deg'
    END IF
    ! Now test if less than the tT width or not
    IF( seccombo(j,y) < angwidthS ) THEN
     ! ACCEPTED...
     IF( vlogfile == 1 ) THEN
      write(11,*) 'Accepted'
     END IF
     ! Get duration
     call duration(e,Pp,SecF(j,3),SecF(y,1),slopedurS(j))
     slopeangS(j) = ABS(SecF(j,3) - SecF(y,1))
     ! Test for a flip
     IF( slopeangS(j) > pi ) THEN
      slopeangS(j) = 2.0D0*pi - slopeangS(j)
     END IF
     IF( ABS(SecF(j,3)-SecF(y,1)) < pi ) THEN
      ! Normal case, assign slopeang as the addition on top of the lowest angle
      IF( SecF(j,3) < SecF(y,1) ) THEN
       slopeangS(j) = SecF(j,3) + slopeangS(j)
      ELSE
       slopeangS(j) = SecF(y,1) + slopeangS(j)
      END IF
     ELSE
      ! Flipped case, one of the angles <360 and one is >0
      IF( SecF(j,3) < SecF(y,1) ) THEN
       slopeangS(j) = SecF(y,1) + slopeangS(j)
      ELSE
       slopeangS(j) = SecF(j,3) + slopeangS(j)
      END IF
      IF( slopeangS(j) > 2.0D0*pi ) THEN
       slopeangS(j) = slopeangS(j) - 2.0D0*pi
      END IF
     END IF
    ELSE
     ! REJECTED...
     IF( vlogfile == 1 ) THEN
      write(11,*) 'Not accepted'
     END IF
    END IF
   END DO ! y loop ends
   ! Output so far...
   IF( vlogfile == 1 ) THEN
    write(11,*) ' '
    write(11,*) 'slopedur(',j,') = ',slopedurS(j),' with midangle = ',slopeangS(j)*(180.D0/pi), 'deg'
   END IF
  END DO ! jloop
  ! Assign S12 and S34
  ! Whatever angle is lowest is assumed to be t12
  IF( ABS(slopeangS(2) - slopeangS(1)) < pi ) THEN
   ! Normal case
   IF( slopeangS(1) < slopeangS(2) ) THEN
    S12 = slopedurS(1); S34 = slopedurS(2)
   ELSE
    S12 = slopedurS(2); S34 = slopedurS(1)
   END IF
  ELSE
   ! Flipped case
   IF( slopeangS(1) < slopeangS(2) ) THEN
    S12 = slopedurS(2); S34 = slopedurS(1)
   ELSE
    S12 = slopedurS(1); T34 = slopedurS(2)
   END IF
  END IF
 ELSE
  S12 = 0.0D0; S34 = 0.0D0
 END IF ! Full transit IF statement

 ! Output secondary results
 IF( calcsec == 1 .AND. vlogfile == 1 ) THEN
  write(11,*) ' '
  write(11,*) 'S12 = ',S12,' seconds'
  write(11,*) 'S34 = ',S34,' seconds'
 END IF
 IF( calcsec == 1 .AND. verbose == 1 ) THEN
  write(*,*) ' '
  write(*,*) 'S12 = ',S12,' seconds'
  write(*,*) 'S34 = ',S34,' seconds'
 END IF

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 2.9: Solution fingerprint
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! [used to probe general properties of solution conditions]

 ! No need to give secflag because it is just the opposite of priflag
 IF( durfinger == 1 ) THEN
  IF ( calcpri == 1 .OR. calcsec == 1 ) THEN
   !write(21,*) e,w,bp,Pp,aR,csc2i,DurChoice,priflag
   write(21,*) e,(w*180.0D0/pi),DurChoice,priflag
  END IF
 END IF

 ELSE ! This is the durexact IF
  
  ! Approx primary duration formula
  IF( calcpri == 1 ) THEN
   DO k=1,3
    call durapprox(e,Pp,rhoP,bP,aR,v(k),T(k))
   END DO
   T12 = 0.5D0*(T(1) - T(3))
   T34 = T12
  END IF
  ! Approx secondary duration formula
  IF( calcsec == 1 ) THEN
   DO k=1,3
    call durapprox(e,Pp,rhoS,bS,aR,v(k),S(k))
   END DO
   S12 = 0.5D0*(S(1) - S(3))
   S34 = S12
  END IF

 END IF

!==================================================================!
!                                                                  !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! SECTION 3                                               !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! Midtime calculations                                    !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!                                                                  !
!==================================================================!

 IF( calcmid == 1 ) THEN
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) '======================================='
   write(11,*) '==========Calculating Mids============='
   write(11,*) '======================================='
  END IF

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 3.1: Define substitutions & precision
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! Jasmine v2.0 used to work by solving the quartic but we found
! that results were numerically unstable surrounding w=0 and w=180.
! To ammend this, Jasmine v2.2 and greater use a series expansion
! to calculate fpri & fsec.

 ! Required substitutions
 cos2i = DCOS(irad)**2
 kk = e*DCOS(w)
 hh = e*DSIN(w)
 hhP1 = hh + 1.0D0
 hhM1 = hh - 1.0D0

 ! Precision aim
 errT = 1.0D-2 ! Error is set to 10ms precision by default
 errF = (2.0D0*pi/Pp)*(DSQRT(1.0D0-e*e)/(rhoP*rhoP))*errT ! Error in terms of f

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 3.2: Primary transit iteration
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

 ! Initial guess, 0th order
 fpri = fPguess
 prefac = kk/hhP1

 ! Go to 1st order
 eta = prefac*cos2i
 fpri = fpri - eta
 IF( vlogfile == 1 ) THEN
 write(11,*) ' '
 write(11,*) 'For primary transit...'
 write(11,*) '1st [O] used, 0th [O] error = ',DABS((eta/errF)*errT),' seconds'
 END IF
 IF( DABS(eta) .GT. errF ) THEN
  ! Go to 2nd order
  eta = (prefac/hhP1)*cos2i**2
  fpri = fpri - eta
  IF( vlogfile == 1 ) THEN
  write(11,*) '2nd [O] used, 1st [O] error = ',DABS((eta/errF)*errT),' seconds'
  END IF
  IF( DABS(eta) .GT. errF ) THEN
   ! Go to 3rd order
   eta = -prefac*( (-6.0D0*hhP1 + kk*kk*(-1.0D0+2.0D0*hh))/&
           (6.0D0*hhP1**3) )*cos2i**3
   fpri = fpri - eta
   IF( vlogfile == 1 ) THEN
   write(11,*) '3rd [O] used, 2nd [O] error = ',DABS((eta/errF)*errT),' seconds'
   END IF
   IF( DABS(eta) .GT. errF ) THEN
    ! Go to 4th order
    eta = -prefac*( (-2.0D0*hhP1 + kk*kk*(-1.0D0+3.0D0*hh))/(2.0D0*hhP1**4) )*cos2i**4
    fpri = fpri - eta
    IF( vlogfile == 1 ) THEN
    write(11,*) '4th [O] used, 3rd [O] error = ',DABS((eta/errF)*errT),' seconds'
    END IF
    IF( DABS(eta) .GT. errF ) THEN
     ! Go to 5th order
     eta = prefac*( (40.0D0*hhP1**2-40.0D0*kk*kk*(-1.0D0+3.0D0*hh+4.0D0*hh*hh) + &
            kk**4*(3.0D0-19.0D0*hh+8.0D0*hh*hh))/(40.0D0*hhP1**6) )*cos2i**5
     fpri = fpri - eta
     IF( vlogfile == 1 ) THEN
     write(11,*) '5th [O] used, 4th [O] error = ',DABS((eta/errF)*errT),' seconds'
     END IF
     IF( DABS(eta) .GT. errF ) THEN
      ! Go to 6th order
      eta = prefac*( (24.0D0*hhP1**2+9.0D0*kk**4*(1.0D0-8.0D0*hh+5.0D0*hh*hh) - &
             40.0D0*kk*kk*(-1.0D0+4.0D0*hh+5.0D0*hh*hh))/(24.0D0*hhP1**7) )*cos2i**6
      fpri = fpri - eta
      IF( vlogfile == 1 ) THEN
      write(11,*) '6th [O] used, 5th [O] error = ',DABS((eta/errF)*errT),' seconds'
      END IF
      IF( DABS(eta) .GT. errF .AND. vlogfile == 1) THEN
       write(11,*) 'Unable to meet precision requirement, stopping at 6th order'
      END IF
     END IF ! End 6th order
    END IF ! End 5th order
   END IF ! End 4th order
  END IF ! End 3rd order
 END IF ! End 2nd order

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 3.3: Secondary transit iteration
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

 ! Initial guess, 0th order
 fsec = fSguess
 prefac = kk/hhM1

 ! Go to 1st order
 eta = -prefac*cos2i
 fsec = fsec - eta
 IF( vlogfile == 1 ) THEN
 write(11,*) ' '
 write(11,*) 'For secondary transit...'
 write(11,*) '1st [O] used, 0th [O] error = ',DABS((eta/errF)*errT),' seconds'
 END IF
 IF( DABS(eta) .GT. errF ) THEN
  ! Go to 2nd order
  eta = -(prefac/hhM1)*cos2i**2
  fsec = fsec - eta
  IF( vlogfile == 1 ) THEN
  write(11,*) '2nd [O] used, 1st [O] error = ',DABS((eta/errF)*errT),' seconds'
  END IF
  IF( DABS(eta) .GT. errF ) THEN
   ! Go to 3rd order
   eta = -prefac*( (-6.0D0*hhM1 + kk*kk*(1.0D0+2.0D0*hh))/&
           (6.0D0*hhM1**3) )*cos2i**3
   fsec = fsec - eta
   IF( vlogfile == 1 ) THEN
   write(11,*) '3rd [O] used, 2nd [O] error = ',DABS((eta/errF)*errT),' seconds'
   END IF
   IF( DABS(eta) .GT. errF ) THEN
    ! Go to 4th order
    eta = prefac*( (-2.0D0*hhM1 + kk*kk*(1.0D0+3.0D0*hh))/(2.0D0*hhM1**4) )*cos2i**4
    fsec = fsec - eta
    IF( vlogfile == 1 ) THEN
    write(11,*) '4th [O] used, 3rd [O] error = ',DABS((eta/errF)*errT),' seconds'
    END IF
    IF( DABS(eta) .GT. errF ) THEN
     ! Go to 5th order
     eta = prefac*( (40.0D0*hhM1**2-40.0D0*kk*kk*(-1.0D0+3.0D0*hh+4.0D0*hh*hh) + &
            kk**4*(3.0D0+19.0D0*hh+8.0D0*hh*hh))/(40.0D0*hhM1**6) )*cos2i**5
     fsec = fsec - eta
     IF( vlogfile == 1 ) THEN
     write(11,*) '5th [O] used, 4th [O] error = ',DABS((eta/errF)*errT),' seconds'
     END IF
     IF( DABS(eta) .GT. errF ) THEN
      ! Go to 6th order
      eta = -prefac*( (24.0D0*hhM1**2+9.0D0*kk**4*(1.0D0+8.0D0*hh+5.0D0*hh*hh) - &
             40.0D0*kk*kk*(-1.0D0-4.0D0*hh+5.0D0*hh*hh))/(24.0D0*hhM1**7) )*cos2i**6
      fsec = fsec - eta
      IF( vlogfile == 1 ) THEN
      write(11,*) '6th [O] used, 5th [O] error = ',DABS((eta/errF)*errT),' seconds'
      END IF
      IF( DABS(eta) .GT. errF .AND. vlogfile == 1 ) THEN
       write(11,*) 'Unable to meet precision requirement, stopping at 6th order'
      END IF
     END IF ! End 6th order
    END IF ! End 5th order
   END IF ! End 4th order
  END IF ! End 3rd order
 END IF ! End 2nd order

 END IF ! calcmid flag


!==================================================================!
!                                                                  !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! SECTION 4                                               !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! Advanced calculations                                   !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!                                                                  !
!==================================================================!

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 4.1: Flux minima timing (Delta-FMT)
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

! |DFMT| is proportional to |Cos[w]|

 IF( calcDFMT == 1 ) THEN
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) '======================================='
   write(11,*) '==========Calculating DFMT============='
   write(11,*) '======================================='
  END IF
  
  ! Choose lowest angle
  IF( ABS(PriF(1,2) - PriF(2,2)) < pi ) THEN
   ! Normal case
   IF( PriF(1,2) < PriF(2,2) ) THEN
    call duration(e,Pp,fpri,PriF(1,2),halfT)
   ELSE
    call duration(e,Pp,fpri,PriF(2,2),halfT)
   END IF
  ELSE
   ! Flip case
   IF( PriF(1,2) < PriF(2,2) ) THEN
    call duration(e,Pp,fpri,PriF(2,2),halfT)
   ELSE
    call duration(e,Pp,fpri,PriF(1,2),halfT)
   END IF
  END IF
  DFMT = 0.5D0*T(2)-halfT
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) 'DFMT = ',DFMT,' seconds'
  END IF
  IF( verbose == 1 ) THEN
   write(*,*) ' '
   write(*,*) 'DFMT = ',DFMT,' seconds'
  END IF
 END IF

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! SUBSECTION 4.2: Ingress/Egress asymmetry (IEA)
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

! |IEA| is proportional to |ecosw| too

 IF( calcIEA == 1 ) THEN
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) '======================================='
   write(11,*) '==========Calculating IEA=============='
   write(11,*) '======================================='
  END IF
  IEA = (T12-T34)/(T12+T34)
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) 'IEA = ',IEA*100.0D0,' %'
  END IF
  IF( verbose == 1 ) THEN
   write(*,*) ' '
   write(*,*) 'IEA = ',IEA*100.0D0,' %'
  END IF
 END IF

!==================================================================!
!                                                                  !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! SECTION 5                                               !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! RV-Transit offset                                     !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!                                                                  !
!==================================================================!

 IF( calcRV == 1 ) THEN
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) '======================================='
   write(11,*) '==========Calculating RV offset========'
   write(11,*) '======================================='
   write(11,*) ' '
  END IF

  ! RV minimum occurs when dZdt is minimized
  ! Officially, four solns exist, two of which are equivalent
  ! The two remaining solns correspond to the two RV minima
  ! From fingerprint testing of v2.0, we know which soln we require.
  fRVmin = -w + ( DACOS(-e*DCOS(w)) )
  ! Add 360 degrees if needed
  DO WHILE( fRVmin .LT. 0.0D0 )
   fRVmin = fRVmin + 2.0D0*pi
  END DO
  IF( vlogfile == 1 ) THEN
   write(11,*) 'fRV = ',fRVmin*(180.0D0/pi),' deg'
   write(11,*) 'c.f. fpri = ',fpri*(180.0D0/pi),' deg'
  END IF

  ! Now calculate time difference, t(RV min) = t(transit) + delta
  ! So if fRVmin < fpri, make it negative
  IF( calcmid == 1 ) THEN
   call duration(e,Pp,fpri,fRVmin,RVoffset)
   IF( DABS(fRVmin-fpri) < pi ) THEN
    ! Normal case
    IF( fRVmin < fpri ) THEN
     RVoffset = - RVoffset
    END IF
   ELSE
    ! Flip case
    IF( fRVmin > fpri ) THEN
     RVoffset = - RVoffset
    END IF
   END IF
  ELSE
   RVoffset = 0.0D0 ! fpri was not calculated
  END IF

  ! Output final 1
  IF( vlogfile == 1 ) THEN
   write(11,*) 't(RV) - t(transit) = ',RVoffset,' s'
   write(11,*) 'Approximate eqn gives ~ ',-(e*DCOS(w)*Pp)/(2.0D0*pi),' s'
  END IF
  IF( verbose == 1 ) THEN
   write(*,*) ' '
   write(*,*) 't(RV) - t(transit) = ',RVoffset,' s'
   write(*,*) 'Approximate eqn gives ~ ',-(e*DCOS(w)*Pp)/(2.0D0*pi),' s'
  END IF

 END IF

!==================================================================!
!                                                                  !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! SECTION 6                                               !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! Sec-Transit Offset                                      !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!                                                                  !
!==================================================================!

 ! Calculate secpri
 IF( calcSECOFF .EQ. 1 ) THEN
  IF( vlogfile == 1 ) THEN
   write(11,*) ' '
   write(11,*) '======================================='
   write(11,*) '=========Calculating SEC offset========'
   write(11,*) '======================================='
   write(11,*) ' '
   write(11,*) 'fpri = ',fpri*(180.0D0/pi),' degrees'
   write(11,*) 'fsec = ',fsec*(180.0D0/pi),' degrees'
   write(11,*) ' '
  END IF

  ! Approximate version
  secpriapprox = DATAN( (e*DCOS(w)) / (DSQRT(1.0D0-e*e)) )
  secpriapprox = secpriapprox + (e*DCOS(w)*DSQRT(1.0D0-e*e))/(1.0D0-e*e*DSIN(w)*DSIN(w))
  secpriapprox = ((secpriapprox/pi) + 0.5D0)*Pp
  ! Exact version
  call duration(e,Pp,fpri,fsec,secpri)
  IF( vlogfile == 1 ) THEN
   write(11,*) 'Sphas = ',secpri/Pp
   write(11,*) 'Sphas ~ ',secpriapprox/Pp
  END IF
  ! Correction
  IF( (secpri/Pp) < 0.5D0 .AND. (secpriapprox/Pp) > 0.5D0 ) THEN
   secpri = Pp - secpri
   IF( vlogfile == 1 ) THEN
    write(11,*) 'Correcting secpri...'
   END IF
  ELSE IF( (secpri/Pp) > 0.5D0 .AND. (secpriapprox/Pp) < 0.5D0 ) THEN
   secpri = Pp - secpri
   IF( vlogfile == 1 ) THEN
    write(11,*) 'Correcting secpri...'
   END IF
  END IF
  ! Output
  IF( vlogfile == 1 ) THEN
   write(11,*) 'Secpri = ',secpri
   write(11,*) 'Secpri ~ ',secpriapprox
  END IF
 ELSE
  secpri = 0.0D0
 END IF

!==================================================================!
!                                                                  !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! SECTION 7                                               !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!   ! Reformatting                                            !    !
!   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!    !
!                                                                  !
!==================================================================!

! This section makes some simple conversions so that Jasmine provides
! an output which Jammi will understand

 ! Assign durations
 tT = T(1); DurP = T(2); tF = T(3)
 sT = S(1); DurS = S(2); sF = S(3)

 ! Assign inclination
 ideg = (180.0D0/pi)*DASIN(DSQRT(1.0D0/csc2i))

 ! Assign rhostar
 rhostar = (Pp/(2.0D0*pi))*(Pp/(2.0D0*pi))*Grv*(1.0D0/(aR**3))
 rhostar = 0.001D0/((4.0D0/3.0D0)*pi*rhostar)

 ! Change RVoffset into days
 RVoffset = 0.000011574074074074073D0*RVoffset

 ! Close the vlog file
 close(11)

! ======================================================
END SUBROUTINE jasmine
! ======================================================

! ======================================================
SUBROUTINE durcoeffs(lam2,e,w,csc2i,S2,Q)

!
! durcoeffs simply saves us having to calculate the quartic 
! coefficients for tT, tF and DurP every time in the main
! code. Instead, we may simply call this code changing
! the parameter S2 = S^2 (e.g. S^2 = (1+p)^2 for tT)
!

implicit none

 ! Inputs
 REAL(8), INTENT(IN) :: lam2,e,w,csc2i,S2
 ! Coefficients
 REAL(8) :: Q0, Q1, Q2, Q3, Q4
 REAL(8), DIMENSION(5), INTENT(OUT) :: Q

 ! Q0
 Q0 = csc2i*(S2-lam2) + lam2*DCOS(w)*DCOS(w)
 Q0 = Q0*Q0
 ! Q1
 Q1 = 2.0D0*csc2i*(S2-lam2) + lam2 + lam2*DCOS(2.0D0*w)
 Q1 = 2.0D0*csc2i*e*S2*Q1
 ! Q2
 Q2 = lam2*(csc2i*(-2.0D0+e*e)*S2 - lam2 + 2.0D0*csc2i*lam2)*DCOS(2.0D0*w)
 Q2 = csc2i*e*e*S2*lam2 - lam2*lam2 + Q2
 Q2 = 2.0D0*csc2i*csc2i*e*e*S2*(3.0D0*S2 - lam2) + Q2
 ! Q3
 Q3 = csc2i*e*e*S2 - lam2*DCOS(2.0D0*w)
 Q3 = 4.0D0*csc2i*e*S2*Q3
 ! Q4
 Q4 = lam2*lam2 - 2.0D0*csc2i*e*e*S2*lam2*DCOS(2.0D0*w)
 Q4 = csc2i*csc2i*e*e*e*e*S2*S2 + Q4
 ! Array it
 Q(1) = Q0
 Q(2) = Q1
 Q(3) = Q2
 Q(4) = Q3
 Q(5) = Q4

END SUBROUTINE durcoeffs
!============================================================

! ======================================================
SUBROUTINE midcoeffs(e,w,csc2i,R)

!
! durcoeffs simply saves us having to calculate the quartic 
! coefficients for tT, tF and DurP every time in the main
! code. Instead, we may simply call this code changing
! the parameter S2 = S^2 (e.g. S^2 = (1+p)^2 for tT)
!

implicit none

 ! Inputs
 REAL(8), INTENT(IN) :: e,w,csc2i
 REAL(8) :: sin2i, cosi2
 ! Coefficients
 REAL(8) :: R0, R1, R2, R3, R4
 REAL(8), DIMENSION(5), INTENT(OUT) :: R

 ! Get sin^2i
 sin2i = 1.0D0/csc2i
 cosi2 = 2.0D0*(1.0D0 - sin2i) - 1.0D0

 ! R0
 R0 = -4.0D0*e*e*(DCOS(w)*DCOS(w)*sin2i - 1.0D0)**2
 R0 = R0 + sin2i*sin2i*DSIN(2.0D0*w)*DSIN(2.0D0*w)
 ! R1
 R1 = -DCOS(w)*DCOS(w) + (DCOS(w)**4)*sin2i + DSIN(w)*DSIN(w)
 R1 = -8.0D0*e*sin2i*R1
 ! R2
 R2 = e*e*(3.0D0+cosi2)*DCOS(2.0D0*w) + cosi2*(-2.0D0+e*e)
 R2 = 2.0D0 + 3.0D0*e*e + R2
 R2 = 4.0D0*(e*e - 0.25D0*R2*sin2i)
 ! R3
 R3 = -2.0D0*e*(3.0D0 + cosi2)*DCOS(2.0D0*w)*sin2i
 R3 = R3 + 4.0D0*e*sin2i*sin2i
 ! R4
 R4 = 4.0D0*sin2i*sin2i
 ! Array it
 R(1) = R0
 R(2) = R1
 R(3) = R2
 R(4) = R3
 R(5) = R4

END SUBROUTINE midcoeffs
!============================================================

! ======================================================
SUBROUTINE quartic(a,b,c,d,e,x1,x2,x3,x4)

! David Kipping, CfA
!
! Quartic solves the quartic eqn analytically.
! The form is assumed to be
!
! 0 == a + b x + c x^2 + d x^3 + e x^4
!
! Quartic returns x1, x2, x3, x4
!

implicit none

 ! Coefficients
 REAL(8), INTENT(IN) :: a, b, c, d, e
 ! First level
 COMPLEX(8) :: p, q, r, s
 REAL(8) :: be, ce, de
 ! Final solutions
 COMPLEX(8), INTENT(OUT) :: x1, x2, x3, x4
 ! Constants
 REAL(8), PARAMETER :: third = 0.3333333333333D0
 REAL(8), PARAMETER :: tom = 1.2599210498948732D0

 ! First level
 p = -2.0D0*c*c*c + 9.0D0*b*c*d - 27.0D0*a*d*d 
 p = p - 27.0D0*b*b*e + 72.0D0*a*c*e
 q = c*c - 3.0D0*b*d + 12.0D0*a*e
 r = p + SQRT(p*p - 4.0D0*q*q*q)
 r = r**third
 s = ((tom*q)/(3.0D0*e*r)) + (r/(3.0D0*tom*e))

 be = b/e
 ce = c/e
 de = d/e

 ! x1 solution
 x1 = -2.0D0*third*ce + 0.25D0*de*de - s
 x1 = 4.0D0*SQRT(x1)
 x1 = ( -8.0D0*be + 4.0D0*ce*de - de*de*de )/(x1)
 x1 = 0.5D0*SQRT(s - x1 - 4.0D0*third*ce + 0.5D0*de*de)
 x1 = -0.25D0*de - x1
 x1 = x1 - 0.5D0*SQRT(-2.0D0*third*ce + 0.25D0*de*de - s)

 ! x2 solution
 x2 = 2.0D0*SQRT(-8.0D0*third*ce + de*de - 4.0D0*s)
 x2 = (8.0D0*be - 4.0D0*ce*de + de*de*de)/x2
 x2 = 0.5D0*SQRT( x2 + s + 0.5D0*de*de - 4.0D0*third*ce )
 x2 = x2 -0.25D0*de
 x2 = x2 - 0.25D0*SQRT(-8.0D0*third*ce + de*de - 4.0D0*s)

 ! x3 solution
 x3 = 2.0D0*SQRT(-8.0D0*third*ce + de*de - 4.0D0*s)
 x3 = (8.0D0*be - 4.0D0*ce*de + de*de*de)/x3
 x3 = 0.5D0*SQRT( -x3 + s + 0.5D0*de*de - 4.0D0*third*ce )
 x3 = -x3 -0.25D0*de
 x3 = x3 + 0.25D0*SQRT(-8.0D0*third*ce + de*de - 4.0D0*s)

 ! x4 solution
 x4 = 2.0D0*SQRT(-8.0D0*third*ce + de*de - 4.0D0*s)
 x4 = (8.0D0*be - 4.0D0*ce*de + de*de*de)/x4
 x4 = 0.5D0*SQRT( -x4 + s + 0.5D0*de*de - 4.0D0*third*ce )
 x4 = x4 -0.25D0*de
 x4 = x4 + 0.25D0*SQRT(-8.0D0*third*ce + de*de - 4.0D0*s)

 ! Check that x is a solution
 !write(*,*) 'Sol x1 gives ',a + b*x1 + c*x1**2 + d*x1**3 + e*x1**4
 !write(*,*) 'Sol x2 gives ',a + b*x2 + c*x2**2 + d*x2**3 + e*x2**4
 !write(*,*) 'Sol x3 gives ',a + b*x3 + c*x3**2 + d*x3**3 + e*x3**4
 !write(*,*) 'Sol x4 gives ',a + b*x4 + c*x4**2 + d*x4**3 + e*x4**4
 !write(*,*) ' '

END SUBROUTINE quartic
! ======================================================

! ======================================================
SUBROUTINE Scheck(f,w,e,csc2i,lam2,vv,n,S2,sol)
! Scheck calculates S^2 as a function of f
! Scheck expects to receive a (n,2)-dimension array of f

implicit none

 ! Controls
 INTEGER :: i, j, n
 ! Inputs
 REAL(8), INTENT(IN) :: w,e,csc2i,lam2,vv
 REAL(8), DIMENSION(n,2), INTENT(IN) :: f
 ! Output
 REAL(8), DIMENSION(n,2), INTENT(OUT) :: S2
 INTEGER, DIMENSION(n), INTENT(OUT) :: sol

  ! Calculate S2 for all
  DO i=1,n
   DO j=1,2
    S2(i,j) = 1.0D0 - ((DSIN(w+f(i,j))*DSIN(w+f(i,j)))/csc2i)
    S2(i,j) = (lam2*S2(i,j))/((1.0D0 + e*DCOS(f(i,j)))**2)
   END DO
  END DO

  ! Now choose the favourite solution (e.g. value closest to 1 for T)
  DO i=1,n
   IF( ABS(S2(i,1) - vv) < ABS(S2(i,2) - vv) ) THEN
    sol(i) = 1
   ELSE
    sol(i) = 2
   END IF
  END DO

END SUBROUTINE Scheck
!============================================================

! ======================================================
SUBROUTINE dScheck(f,w,e,csc2i,lam2,n,dSdf,sol)
! dScheck calculates dS/df as a function of f
! dScheck expects to receive a (n,2)-dimension array of f

implicit none

 ! Controls
 INTEGER :: i, j, n
 ! Inputs
 REAL(8), INTENT(IN) :: w,e,csc2i,lam2
 REAL(8) :: lam
 REAL(8), DIMENSION(n,2), INTENT(IN) :: f
 ! Output
 REAL(8), DIMENSION(n,2), INTENT(OUT) :: dSdf
 INTEGER, DIMENSION(n), INTENT(OUT) :: sol

  lam = DSQRT(lam2)
  ! Calculate dSdF for all
  DO i=1,n
   DO j=1,2
    dSdf(i,j) = -(lam*((e - 2.0D0*csc2i*e)*DSIN(f(i,j)) + &
       DSIN(2.0D0*(f(i,j)+w)) + e*DSIN(f(i,j) + 2.0D0*w)))/&
       (2.0D0*csc2i*(1.0D0 + e*DCOS(f(i,j)))**2*DSQRT(1.0D0 - DSIN(f(i,j)+w)**2/csc2i))
   END DO
  END DO

  ! Now choose the parts of each pair which give lowest dSdf
  DO i=1,n
   IF( ABS(dSdf(i,1)) < ABS(dSdf(i,2)) ) THEN
    sol(i) = 1
   ELSE
    sol(i) = 2
   END IF
  END DO

END SUBROUTINE dScheck
!============================================================

! ======================================================
SUBROUTINE d2Scheck(f,w,e,csc2i,lam2,n,Q)!,sol)
! d2Scheck calculates d^2S/df^2 as a function of f
! d2Scheck expects to receive a (n,2)-dimension array of f

implicit none

 ! Controls
 INTEGER :: i, j, n
 ! Inputs
 REAL(8), INTENT(IN) :: w,e,csc2i,lam2
 REAL(8) :: lam
 REAL(8), DIMENSION(n,2), INTENT(IN) :: f
 ! Output
 REAL(8), DIMENSION(n,2), INTENT(OUT) :: Q

  lam = DSQRT(lam2)
  ! Calculate dSdSdF for all
  DO i=1,n
!    j = take(i)
   DO j=1,2
    Q(i,j) = (lam*((-1.0D0 + csc2i)*csc2i*e**2*DCOS(f(i,j))**2 - csc2i*DCOS(f(i,j)+w)**2 + &
    2.0D0*csc2i**2*e**2.0D0*DSIN(f(i,j))**2 + csc2i*DSIN(f(i,j) + w)**2 - &
    4.0D0*csc2i*e**2*DSIN(f(i,j))**2*DSIN(f(i,j) + w)**2 + &
    2.0D0*e*(1.0D0 + e*DCOS(f(i,j)))*DCOS(f(i,j) + w)*DSIN(f(i,j))*DSIN(f(i,j) + w)**3 - &
    DSIN(f(i,j) + w)**4 + & 
    2.0D0*e**2*DSIN(f(i,j))**2*DSIN(f(i,j) + w)**4 - csc2i*e*DSIN(f(i,j))*&
    DSIN(2.0D0*(f(i,j) + w)) - &
    e*DCOS(f(i,j))*(2.0D0*csc2i*DCOS(f(i,j) + w)**2 + DSIN(f(i,j) + w)**4 + &
    csc2i*(-csc2i + e*DSIN(f(i,j))*DSIN(2.0D0*(f(i,j) + w))))))/&
    (csc2i*(1.0D0 + e*DCOS(f(i,j)))**3*(csc2i - DSIN(f(i,j) + w)**2)*&
    DSQRT(1.0D0 - DSIN(f(i,j) + w)**2/csc2i))
   END DO
  END DO

END SUBROUTINE d2Scheck
!============================================================

! =======================================================
SUBROUTINE durapprox(e,Pp,rho,b,aR,vv,T_1)
! Uses T1 eqn to get duration
! This will work for both primary and secondary eclipses

implicit none

 REAL(8), INTENT(IN) :: e, Pp, rho, b, aR, vv
 REAL(8), INTENT(OUT) :: T_1
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0

 IF( b < DSQRT(vv) ) THEN
  T_1 = (vv - b*b)/(rho*rho*aR*aR - b*b)
  T_1 = rho*rho*Pp*DASIN(DSQRT(T_1))
  T_1 = T_1/(pi*DSQRT(1.0D0 - e*e))
 ELSE
  T_1 = 0.0D0
 END IF

END SUBROUTINE durapprox
!============================================================

! =======================================================
SUBROUTINE duration(e,Pp,f1,f2,T)
! Uses T eqn to get duration

implicit none

 INTEGER :: i
 REAL(8), INTENT(IN) :: e, Pp, f1, f2
 REAL(8), DIMENSION(2) :: f
 REAL(8) :: Doffset, J
 REAL(8), DIMENSION(2) :: D
 REAL(8), INTENT(OUT) :: T
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0

 f(1) = f1
 f(2) = f2
 J = (2.0D0*pi*DSQRT(1.0D0-e*e))/Pp
 Doffset = 2.0D0*DSQRT(1.0D0-e*e)*pi

 ! Get D
 DO i=1,2
  D(i) = (DSQRT(1.0D0-e)/DSQRT(1.0D0+e))*DTAN(0.5D0*f(i))
  D(i) = 2.0D0*DSQRT(1.0D0-e*e)*DATAN(D(i))
  D(i) = D(i) - ((e*(1.0D0-e*e)*DSIN(f(i)))/(1.0D0+e*DCOS(f(i))))
  IF(f(i)>pi) THEN
   D(i) = D(i) + Doffset
  END IF
 END DO

 ! Calculate T
 T = (ABS(D(2) - D(1)))/J

 ! Correct if overflows
 IF( T > 0.5D0*Pp ) THEN
  T = Pp - T
 END IF

END SUBROUTINE duration
!============================================================

!============================================================
END MODULE jasminemod
!============================================================
