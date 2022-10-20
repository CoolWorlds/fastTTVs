MODULE like

use params
use transitmod
use radialmod
use inversebetamod

implicit none
      
contains
      
!=======================================================================

SUBROUTINE slikelihood(R,slhood)
         
	implicit none
      
	REAL(8), DIMENSION(nest_nPar) :: R
    REAL(8) :: slhood, loglike
	INTEGER :: i, ifault
	REAL(8) :: aRtemp, rhotemp, Pdaystemp, logpenalty

    ! Scaling of the parameters from hypercube
	DO i = 1, sdim
          IF( Rflag(i) .EQ. 0 .OR. Rflag(i) .EQ. 2 ) THEN ! Modes 0 and 2
            ! Uniform: Rmax = max, Rmin = min
	    R(i) = Rmin(i) + (Rmax(i)-Rmin(i))*R(i)
          ELSE IF( Rflag(i) .EQ. 3 ) THEN ! Mode 3
            ! Gaussian: Rmax = mean, Rmin = stdev
            R(i) = Rmax(i) + roottwo*Rmin(i)*inverf(-1.0D0+2.0D0*R(i))
          ELSE IF( Rflag(i) .EQ. 4 ) THEN ! Mode 4
            ! Jeffrey's: Rmax = max, Rmin = min
            R(i) = ( Rmax(i)**R(i) )*( Rmin(i)**(1.0D0-R(i)) )
          ELSE IF( Rflag(i) .EQ. 5 ) THEN ! Mode 5
            ! Modified Jeffrey's: Rmax = max, Rmin = inflection point
            R(i) = -( Rmin(i)**(1.0D0-R(i)) )*( Rmin(i)**R(i) - ( Rmin(i)+Rmax(i) )**R(i) )
          ELSE IF( Rflag(i) .EQ. 6 ) THEN ! Mode 6
            ! Beta prior: Rmax = b, Rmin = a
            beta_log = alngam(Rmin(i),ifault) + alngam(Rmax(i),ifault) &
                       - alngam(Rmin(i)+Rmax(i),ifault)
            write(*,*) 'beta_log = ',beta_log
            R(i) = xinbta ( Rmin(i), Rmax(i), beta_log, R(i), ifault )
          END IF
	END DO

    ! Call transit to get chi^2
    call models(R,loglike,0)

	rhotemp = R(2)
	Pdaystemp = R(4)
    aRtemp = rhotemp*Grvx*Pdaystemp**2
    aRtemp = aRtemp**third
	logpenalty = (1.0D0/aRtemp)/( DLOG(20000.0D0) - DLOG(2.0D0) )
	logpenalty = DLOG(logpenalty)
	IF( aRtemp .GE. 20000.0D0 ) THEN
		logpenalty = -HUGE(1.0D0)
	END IF

	slhood = loglike + logpenalty

END SUBROUTINE slikelihood
      
!=======================================================================

!! ======================================================================
SUBROUTINE models(Rvec,loglike,showpri)

!implicit none

 REAL(8), DIMENSION(nest_nPar), INTENT(IN) :: Rvec   ! Fitted-parameter vector
 INTEGER :: showpri, showrv
 REAL(8), DIMENSION(nplen) :: resP
 REAL(8) :: loglike, loglikeP, loglikeP_lc, loglikeP_sc, loglikeR
 INTEGER :: i, j, nplen_lc, nplen_sc
 REAL(8) :: jitter
 REAL(8) :: chi2RV

 ! === Call transit to primary transit ===

 call transit(Rvec,nplen,resP,showpri,&
              tp,fp,sigfp,epochp_seq,fpwei,sigfpwei,tpdeviant,&
              NresamP,integ_bigP,.FALSE.)

 IF( priflag .EQ. 1 ) THEN
   IF( NresamP .EQ. 1 ) THEN
     ! Compute log likelihood under simple conditions of global noise terms
     loglikeP = 0.0D0
     DO i=1,nplen
       loglikeP = loglikeP + resP(i)**2*sigfpwei(i) + logsigfp(i)
     END DO
     loglikeP = -0.5D0*nplen*LogTwoPi - 0.5D0*loglikeP
   ELSE
     ! == LC ==
     loglikeP_lc = 0.0D0
     nplen_lc = 0
     DO i=1,nplen
       IF( sigfp(i) .LE. flick ) THEN !Take only the LC data
         loglikeP_lc = loglikeP_lc + resP(i)**2*sigfpwei(i) + logsigfp(i)
         nplen_lc = nplen_lc + 1
       END IF
     END DO
     loglikeP_lc = -0.5D0*nplen_lc*LogTwoPi - 0.5D0*loglikeP_lc
     ! == SC ==
     loglikeP_sc = 0.0D0
     nplen_sc = 0
     DO i=1,nplen
       IF( sigfp(i) .GT. flick ) THEN !Take only the SC data
         loglikeP_sc = loglikeP_sc + resP(i)**2*sigfpwei(i) + logsigfp(i)
         nplen_sc = nplen_sc + 1
       END IF
     END DO
     loglikeP_sc = -0.5D0*nplen_sc*LogTwoPi - 0.5D0*loglikeP_sc
     ! == BOTH ==
     loglikeP = loglikeP_lc + loglikeP_sc
   END IF
 ELSE
   ! Log[like] = 1 => like = 1 i.e. all trials look good!
   loglikeP = 1.0D0 
 END IF

 ! === Call radial ===

 IF ( rvflag .EQ. 1 ) THEN

 jitter = 0.0D0 !Rvec(8)
 call radial(Rvec,chi2RV)

 ! Compute log likelihood under simple conditions of global noise terms
 loglikeR = chi2RV
 DO i=1,nrlen
   loglikeR = loglikeR + DLOG( sigrv(i)**2 + jitter**2 )
 END DO
 loglikeR = -0.5D0*nrlen*LogTwoPi - 0.5D0*loglikeR
 ELSE
   loglikeR = 0.0D0
 END IF

 ! === Sum up all loglikes (product of likes) ===
 loglike = loglikeP + loglikeR

END SUBROUTINE models
! ======================================================================

! ======================================================================
FUNCTION inverf(x)

 implicit none

 REAL(8) :: x
 REAL(8), PARAMETER :: awil = 0.14001228868666646D0
 REAL(8), PARAMETER :: bwil = 4.546884979448289D0
 REAL(8) :: factor, xsq, inverf

 IF( x .LT. 0.0D0 ) THEN
  factor = -1.0D0
 ELSE
  factor = 1.0D0
 END IF

 xsq = 1.0D0 - x**2
 x = bwil + 0.5D0*DLOG(xsq)
 x = DSQRT( x**2 - (DLOG(xsq)/awil) ) - x
 inverf = factor*DSQRT(x)

END FUNCTION
! ======================================================================

END MODULE like

