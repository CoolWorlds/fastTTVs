MODULE transitmod

use params
use planmod
use jasminemod
implicit none
      
contains
      
      
!=======================================================================
SUBROUTINE transit(Rin,nz,ressy,show,&
                   tobs,fobs,sigfobs,epoch_seq,fobswei,sigfobswei,tdeviant,&
                   Nresam,integ_big,seccy)
         
 implicit none
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 REAL(8), DIMENSION(nest_nPar) :: Rin           ! R in vector
 INTEGER, INTENT(IN) :: nz, Nresam, show
 LOGICAL, INTENT(IN) :: seccy          
 REAL(8), INTENT(IN) :: integ_big
 REAL(8) :: chi2                                ! merit function of fit
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 INTEGER :: i, j
 ! Data terms
 REAL(8), DIMENSION(nz), INTENT(IN) :: tobs, fobs, sigfobs
 REAL(8), DIMENSION(nz), INTENT(IN) :: fobswei, sigfobswei, tdeviant
 INTEGER, DIMENSION(nz), INTENT(IN) :: epoch_seq
 REAL(8), DIMENSION(nz) :: tdiff, flux, ressy
 REAL(8), DIMENSION(OOTlen) :: OOTvec, LINvec, QADvec
 REAL(8), DIMENSION(taulen) :: tauvec
 !!REAL(8), DIMENSION(OOSlen) :: OOSvec
 ! Linear minimization terms
 REAL(8), DIMENSION(OOTlen) :: obsmod1, mod2
 REAL(8), DIMENSION(OOTlen) :: obsmod1tim1, mod2tim1, mod2tim2
 REAL(8), DIMENSION(OOTlen) :: obsmod2, obsmod2tim1, obsmod2tim2, mod4
 REAL(8), DIMENSION(OOTlen) :: mod4tim1, mod4tim2, mod4tim3, mod4tim4
 REAL(8), DIMENSION(OOTlen) :: X1, X2, X3, X4, X5, X6, X7
 INTEGER, DIMENSION(OOTlen) :: epochlength
 ! Parameters
 REAL(8) ::  p, rhomstar, bp, Pdays, gamglobal, samrecip, u1, tmid, w1, w2!, u1pu2
 ! Derived parameters
 REAL(8) :: u2, wrad, wavP, wavS, e, aR, rhoP, rhoS
 REAL(8) :: secpri, tsec, secphas, fpri
 REAL(8) :: tT, tF, DurP, t12, t34, rhostar, ideg, RVoffset, s14
 LOGICAL :: process
 ! Blending stuff
 REAL(8), DIMENSION(nz) :: gammy
 REAL(8), DIMENSION(nz*Nresam) :: gammy_big
 REAL(8) :: gamrecip, sam
 ! Explosion variables
 INTEGER :: k_big, k_max, nz_big
 REAL(8) :: Ndash, integ_small
 ! Explosion flag array
 INTEGER, DIMENSION(nz) :: explo
 ! Unflattened arrays
 REAL(8), DIMENSION(Nresam,nz) :: t_arr, flux_arr
 ! Flattened arrays
 REAL(8), DIMENSION(nz*Nresam) :: t_big
 REAL(8), DIMENSION(nz*Nresam) :: mulimb0!, b0
 INTEGER, DIMENSION(nz*Nresam) :: epoch_seq_big
 ! show=1 variables
 REAL(8), DIMENSION(nz) :: res
 REAL(8) :: time, epoch, epochmid, tmidfold, tsecfold

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 1.0 DECLARATIONS ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 p = Rin(1)   ! Planet's size
 rhomstar = Rin(2)  ! rho_{*}
 bp = Rin(3)  ! Barycentre's impact parameter
 Pdays = Rin(4)     ! Barycentre's period [days]
 gamglobal = 1.0D0  ! Blending factor
 samrecip = 1.0D0 !DABS(Rin(7)) ! Fp/F*
 w1 = Rin(6)   ! Limb darkening w1
 w2 = Rin(7) ! Limb darkening w2
 tmid = Rin(5)     ! Barycentre's transit time
 e = 2.0D-8  ! Barycentre's e
 wrad = 0.7853981633974483D0 ! Barycentre's w DATAN2(hb,kb)
 !write(*,*) 'e = ',e

 ! tauarray
 IF( globalflag ) THEN
   ! donothing
 ELSE
   DO i=1,taulen
     tauvec(i) = Rin(nparamorig+i)
   END DO
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 2.0 CONVERT FITTED PARAMETERS ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 nz_big = nz*Nresam

 ! Calculate wavP
 wavP = fivehalfpi - wrad
 IF( wavP .GT. twopi ) THEN
   wavP = wavP - twopi
 END IF
 wavS = wavP + pi
 IF( wavS .GT. twopi ) THEN
   wavS = wavS - twopi
 END IF
 ! Calculate varrhoP and varrhoS
 rhoP = 1.0D0 - e**2
 rhoS = rhoP/(1.0D0+e*DCOS(wavS))
 rhoP = rhoP/(1.0D0+e*DCOS(wavP))

 ! Convert rhomstar to aR (keep units in days)
 !aR = (DSQRT(rhomstar**3)*Grv*(Pdays*86400.0D0)**2)/(3.0D0*pi)
 aR = rhomstar*Grvx*Pdays**2
 aR = aR**third

 ! Get u2
! u2 = u1pu2 - u1
! u1 = w1*costheta + w2*sintheta
! u2 = -w1*sintheta + w2*costheta
 u2 = DSQRT(w1) !u1+u2
 u1 = 2.0D0*u2*w2
 u2 = u2 - u1
 ! Override secondary eclipses to have no LD
 IF( seccy ) THEN
   u1 = 0.0D0
   u2 = 0.0D0
 END IF

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 2.1 REJECT TRIALS ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 process = .TRUE.

 ! Contact star-planet binaries                                                
 IF( aR*(1.0D0-e) .LT. 2.0D0 ) THEN
   process = .FALSE.
 END IF

 ! Delete excessive eccentricities                                               
 IF( e .GT. 0.95D0 ) THEN
   process = .FALSE.
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 3.0 TIME ARRAY ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 IF( process ) THEN

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 3.1 Offset time array for mid-transit time
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 DO i = 1,nz
   IF( globalflag .OR. seccy ) THEN
     ! Use a global tmid time
     tdiff(i) = tobs(i) - tmid
   ELSE
     ! Use individual transit times
     tdiff(i) = tobs(i) - tauvec(epoch_seq(i))
     !tdiff(i) = tobs(i) - tmid ! Temp holding until intelligence.o included
   END IF
 END DO

 ! Assign quarter-to-quarter blending factors
 gammy(:) = 1.0D0
 DO i=1,nz
   DO j=1,NQ
     IF( tobs(i) .GE. QXstart(j) .AND. tobs(i) .LT. QXend(j) ) THEN
       gammy(i) = gam(j)
     END IF
   END DO
 END DO

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 3.2 Jasmine Call
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~    

 ! We call Jasmine here so that we may implement
 ! selecive resampling

 ! Jasmine calculates various durations in an exact manner
 IF( seccy ) THEN
   ! Use Jasmine's secpri value to get tsec & secphas
   call jasmine(bp,p,e,wrad,aR,Pdays,0,&
                tT,tF,DurP,t12,t34,rhostar,&
                ideg,secpri,s14,RVoffset,fpri)
 ELSE
   secpri = 0.5D0*86400.0D0*Pdays
   fpri = 1.5707963267948966D0 - wrad
 END IF
 tsec = (secpri/86400.0D0) + tmid
 secphas = secpri/(86400.0D0*Pdays)

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 3.3 Explode the time array
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 ! We now have to 'explode' the time array
 ! This is part of the process for accounting for the integration time
 ! of each time stamp.
 ! CAVEATS:
 ! * If m=0, no integration is required.
 ! * If we are in OOT, no integation is required (we define intransit times stamps
 !   as those < 1.1*0.5*(t_T + integration_time) (the 1.1 gives a 10% safety net)
 ! * Points which are exploded are assigned a logic flag 1 in the exploded(i) array

 IF( Nresam .GT. 1 ) THEN ! We only need to do this if Nresam > 1
  ! Stage one, create 2D time array
  k_big = 0
  Ndash = 0.5D0*(Nresam+1.0D0)     ! These two are commonly used
  integ_small = integ_big/Nresam   ! later, so easier to define here
  DO i=1,nz
   ! You add a 2nd condition here eg selective resampling, SC/LC mixed data, etc
   IF( sigfobs(i) .LE. flick ) THEN ! All data before this point is LC
    explo(i) = 1 ! Explosion is occuring, flag it
    DO j=1,Nresam
     ! Performing explosion
     k_big = k_big + 1
     t_arr(j,i) = tdiff(i) + (j-Ndash)*integ_small
     ! Stage two, flatten the 2D array into a 1D array
     t_big(k_big) = t_arr(j,i)
     gammy_big(k_big) = gammy(i)
     epoch_seq_big(k_big) = epoch_seq(i)
    END DO
   ELSE
    k_big = k_big + 1
    explo(i) = 0 ! No explosion occured for this timestamp
    t_big(k_big) = tdiff(i)
    gammy_big(k_big) = gammy(i)
    epoch_seq_big(k_big) = epoch_seq(i)
   END IF
  END DO
  k_max = k_big
 ELSE
  ! Infinitessimal integration time => go as normal
  DO i=1,nz
   t_big(i) = tdiff(i)
   gammy_big(i) = gammy(i)
   explo(i) = 0
   epoch_seq_big(i) = epoch_seq(i)
  END DO
  k_max = nz
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 4.0 GENERATE LIGHTCURVE ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 4.1 Main call to PLAN
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 call plan(t_big,Pdays,0.0D0,p,aR,e,wrad,bp,&
      u1,u2,fpri,k_max,mulimb0)

 !  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 4.2 Transformations
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 ! Define reciprocals of various parameters to speed up computation
 sam = 1.0D0/samrecip
 gamrecip = 1.0D0/gamglobal ! gam = blending factor
 ! PERFORM MODEL -> OBSERVATIONS TRANSFORMATIONS
 IF( seccy ) THEN
   DO i=1,k_max
   ! i) Secondary Eclipse Transformation
     mulimb0(i) = mulimb0(i) + sam - 1.0D0
     mulimb0(i) = mulimb0(i)*samrecip
   END DO
 END IF
 DO i=1,k_max
   ! ii) Blending Transformation
   mulimb0(i) = mulimb0(i) + gammy_big(i) - 1.0D0
   mulimb0(i) = mulimb0(i)/gammy_big(i)
   mulimb0(i) = mulimb0(i) + gamglobal - 1.0D0
   mulimb0(i) = mulimb0(i)*gamrecip
 END DO

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 4.3 Implode the flux array
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 ! Implode the flux array, using standard binning
 ! First stage is to un-flatten the flux array from a
 ! a 1D vector to a 2D array

 IF( Nresam .GT. 1 ) THEN  ! Nresam>1 => implosions required
  k_big = 0
  DO i=1,nz
   ! For ith point, is an implosion needed?
   IF( explo(i) == 1 ) THEN
    ! Stage one, un-flatten the flux array
    DO j=1,Nresam
     k_big = k_big + 1
     flux_arr(j,i) = mulimb0(k_big)
    END DO
    ! Stage two, perform binning
    flux(i) = 0.0D0
    DO j=1,Nresam
     flux(i) = flux(i) + flux_arr(j,i)
    END DO
    flux(i) = flux(i)/Nresam
   ELSE
    k_big = k_big + 1
    flux(i) = mulimb0(k_big)
   END IF
  END DO
 ELSE
 ! Infinitessimal integration time => go as normal
  DO i=1,nz
   flux(i) = mulimb0(i)
  END DO
 END IF

 ! Linear optimization of OOTsfobswei, sigfobswei
 ! Constant terms
 IF( linmode .GE. 0 ) THEN
   obsmod1(:) = 0.0D0
   mod2(:) = 0.0D0
 END IF
 ! Linear terms
 IF( linmode .GE. 1 ) THEN
   obsmod1tim1(:) = 0.0D0
   mod2tim1(:) = 0.0D0
   mod2tim2(:) = 0.0D0
   epochlength(:) = 0
 END IF
 ! Quadratic terms
 IF( linmode .GE. 2 ) THEN
   obsmod2(:) = 0.0D0
   obsmod2tim1(:) = 0.0D0
   obsmod2tim2(:) = 0.0D0
   mod4tim1(:) = 0.0D0
   mod4tim2(:) = 0.0D0
   mod4tim3(:) = 0.0D0
   mod4tim4(:) = 0.0D0
   mod4(:) = 0.0D0
 END IF
 DO j=1,OOTlen
   DO i=1,nz
     IF( epoch_seq(i) .EQ. j ) THEN
       ! Constant terms
       IF( linmode .GE. 0 ) THEN
         obsmod1(j) = obsmod1(j) + fobswei(i)*flux(i)
         mod2(j) = mod2(j) + flux(i)**2*sigfobswei(i)
       END IF
       ! Linear terms
       IF( linmode .GE. 1 ) THEN
         epochlength(j) = epochlength(j) + 1
         obsmod1tim1(j) = obsmod1tim1(j) + fobswei(i)*flux(i)*tdeviant(i)
         mod2tim1(j) = mod2tim1(j) + flux(i)**2*tdeviant(i)*sigfobswei(i)
         mod2tim2(j) = mod2tim2(j) + flux(i)**2*tdeviant(i)**2*sigfobswei(i)
       END IF
       ! Quadratic terms
       IF( linmode .GE. 2 ) THEN
         obsmod2(j) = obsmod2(j) + fobswei(i)*flux(i)**2
         obsmod2tim1(j) = obsmod2tim2(j) + fobswei(i)*flux(i)**2*tdeviant(i)
         obsmod2tim2(j) = obsmod2tim2(j) + fobswei(i)*flux(i)**2*tdeviant(i)**2
         mod4(j) = mod4(j) + flux(i)**4*sigfobswei(i)
         mod4tim1(j) = mod4tim1(j) + flux(i)**4*tdeviant(i)*sigfobswei(i)
         mod4tim2(j) = mod4tim2(j) + flux(i)**4*tdeviant(i)**2*sigfobswei(i)
         mod4tim3(j) = mod4tim3(j) + flux(i)**4*tdeviant(i)**3*sigfobswei(i)
         mod4tim4(j) = mod4tim4(j) + flux(i)**4*tdeviant(i)**4*sigfobswei(i)
       END IF
     END IF
   END DO
   ! Combination terms
   IF( linmode .GE. 2 ) THEN
     X1(j) = obsmod2(j)*mod4tim3(j) - obsmod2tim1(j)*mod4tim2(j)
     X2(j) = mod4tim2(j)**2 - mod4tim1(j)*mod4tim3(j)
     X3(j) = mod4tim1(j)*mod4tim3(j) - mod4tim2(j)**2
     X4(j) = mod4tim1(j)*mod4tim2(j) - mod4(j)*mod4tim3(j)
     X5(j) = obsmod2(j)*mod4tim4(j) - obsmod2tim2(j)*mod4tim2(j)
     X6(j) = mod4tim2(j)*mod4tim3(j) - mod4tim1(j)*mod4tim4(j)
     X7(j) = mod4tim2(j)**2 - mod4(j)*mod4tim4(j)
   END IF
   IF( linmode .EQ. 0 ) THEN
     ! Constant minimization
     OOTvec(j) = obsmod1(j)/mod2(j)
   ELSE IF( linmode .EQ. 1 ) THEN
     IF( epochlength(j) .GE. 2 ) THEN
       ! Linear minimization
       OOTvec(j) = ( obsmod1tim1(j)*mod2tim1(j) - obsmod1(j)*mod2tim2(j) ) &
                   / ( mod2tim1(j)**2 - mod2(j)*mod2tim2(j) )
       LINvec(j) = ( obsmod1(j)/mod2tim1(j) ) &
                   - ( OOTvec(j)*mod2(j)/mod2tim1(j) )
     ELSE IF( epochlength(j) .EQ. 1 ) THEN
       OOTvec(j) = obsmod1(j)/mod2(j)
       LINvec(j) = 0.0D0
     END IF
   ELSE IF( linmode .EQ. 2 ) THEN
     IF( epochlength(j) .GE. 3 ) THEN
       OOTvec(j) = (X3(j)*X5(j) + X1(j)*X6(j))/(X2(j)*X7(j)-X4(j)*X6(j))
       LINvec(j) = OOTvec(j)*(X4(j)/X3(j)) - (X1(j)/X2(j))
       QADvec(j) = ( (X2(j)*X5(j)-X1(j)*X6(j))/(X2(j)*X7(j)-X4(j)*X6(j)) )*&
                   ( (mod4(j)/mod4tim2(j)) - &
                     ( (X4(j)*mod4tim1(j))/(X2(j)*mod4tim2(j)) ) ) + &
                   (obsmod2(j)/mod4tim2(j)) + &
                   ( (X1(j)*mod4tim1(j))/(X2(j)*mod4tim2(j)) )
     ELSE IF( epochlength(j) .EQ. 2 ) THEN
       OOTvec(j) = ( obsmod1tim1(j)*mod2tim1(j) - obsmod1(j)*mod2tim2(j) ) &
                   / ( mod2tim1(j)**2 - mod2(j)*mod2tim2(j) )
       LINvec(j) = ( obsmod1(j)/mod2tim1(j) ) &
                   - ( OOTvec(j)*mod2(j)/mod2tim1(j) )
       QADvec(j) = 0.0D0
     ELSE IF( epochlength(j) .EQ. 1 ) THEN
       OOTvec(j) = obsmod1(j)/mod2(j)
       LINvec(j) = 0.0D0
       QADvec(j) = 0.0D0
     END IF
   END IF
 END DO

 ! Now normalize
 DO i=1,nz
   IF( linmode .EQ. 0 ) THEN
     ! Constant minimization
     flux(i) = flux(i)*OOTvec(epoch_seq(i))
   ELSE IF( linmode .EQ. 1 ) THEN
     ! Linear minimization
     flux(i) = flux(i)*( OOTvec(epoch_seq(i)) &
               + LINvec(epoch_seq(i))*tdeviant(i) )
   ELSE IF( linmode .EQ. 2 ) THEN
     ! Linear minimization
     flux(i) = flux(i)*( OOTvec(epoch_seq(i)) &
               + LINvec(epoch_seq(i))*tdeviant(i) &
               + QADvec(epoch_seq(i))*tdeviant(i)**2 )
   END IF
 END DO

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 5.0 PRINT RESULTS ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 res(:) = 0.0D0
 ! PRIMARY TRANSIT
 IF (show == 1) then
 ! First we do full lightcurve
   open(unit=91,file='PRI_full.jam')
   DO i=1,nz
     IF( globalflag .OR. seccy ) THEN
       time = tdiff(i) + tmid
     ELSE
       time = tdiff(i) + tauvec(epoch_seq(i))
     END IF
     res(i) = ( fobs(i) - flux(i) )*1.0D6
     write(91,*) time,fobs(i),sigfobs(i),flux(i)
   END DO
   close(91)
 ! Next we do the the folded LC
   epochmid = tmid/Pdays
   tmidfold = tmid - epochmid*Pdays
   open(unit=92,file='PRI_fold.jam')
   DO i=1,nz
     IF( globalflag .OR. seccy ) THEN
       time = tdiff(i) + tmid
     ELSE
       time = tdiff(i) + tauvec(epoch_seq(i))
     END IF
     epoch = time/(Pdays)
     time = time - epoch*Pdays
     time = time - tmidfold
     IF( time .GE. 0.5D0*Pdays) THEN
       time = time - Pdays
     ELSE IF( time .LE. (-0.5D0*Pdays)) THEN
       time = time + Pdays
     END IF
     write(92,*) time,fobs(i),sigfobs(i),flux(i)
   END DO
   close(92) 
 END IF
 ! SECONDARY ECLIPSE
 IF (show == 2) THEN
   ! First we do full lightcurve
   open(unit=93,file='SEC_full.jam')
   DO i=1,nz
     time = tdiff(i) + tmid
     res(i) = ( fobs(i) - flux(i) )*1.0D6
     write(93,*) time,fobs(i),sigfobs(i),flux(i)
   END DO
   close(93)
   ! Next we do the the folded LC
   epochmid = tsec/Pdays
   tsecfold = tsec - epochmid*Pdays
   open(unit=94,file='SEC_fold.jam')
   DO i=1,nz
     time = tdiff(i)+tmid
     epoch = time/(Pdays)
     time = time - epoch*Pdays
     time = time - tsecfold
     IF( time .GE. 0.5D0*Pdays) THEN
       time = time - Pdays
     ELSE IF( time .LE. (-0.5D0*Pdays)) THEN
       time = time + Pdays
     END IF
     write(94,*) time,fobs(i),sigfobs(i),flux(i)
   ENDDO
   close(94) 
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 6.0 COMPUTE CHI^2 ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 ! 6.1 Basic chi^2
 ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

 ! Calculated flux = flux(i)
 ! Observed flux = fobs(i)
 ! Observed error = sigfobs(i)
 chi2=0.0D0
 DO i=1,nz
   chi2 = chi2 + ((flux(i) - fobs(i))/sigfobs(i))**2
 END DO
 !write(*,*) 'chi2 = ',chi2,MAXVAL(epoch_seq(:))

 DO i=1,nz
   ressy(i) = flux(i) - fobs(i)
 END DO

 ELSE 
   ! process = .FALSE.
   ressy(:) = 1.0D0
 END IF

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 ! === 7.0 CLOSE PROGRAM ===
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

END SUBROUTINE transit
!=======================================================================

END MODULE transitmod
