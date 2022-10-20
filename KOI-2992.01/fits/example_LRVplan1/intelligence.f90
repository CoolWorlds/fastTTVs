MODULE intelmod

! INTELLIGENCE MODULE V1.0
! This module contain all subroutines for intelligence.
! intelligence refers to automatic assigment of Vsys, OOT, OOS parameters
!
! AUTHOR: David Kipping, Harvard-Smithsonian Center for Astrophysics
! CONTACT: dkipping@cfa.harvard.edu
! UPDATES: www.davidkipping.co.uk
!

 implicit none

 CONTAINS

! =======================================================
SUBROUTINE iVsys1(rvmode,nrlen,Vsyslen,allowedrvmodes)
! covariance computes the covariance matrix of Z

implicit none

 INTEGER, DIMENSION(nrlen), INTENT(IN) :: rvmode
 INTEGER, INTENT(IN) :: nrlen
 INTEGER :: i, j
 INTEGER, INTENT(OUT) :: Vsyslen
 INTEGER, DIMENSION(2,10), INTENT(OUT) :: allowedrvmodes

  ! Count total number of RV modes present (maximum allowed is 10)
  DO j=0,9
   allowedrvmodes(1,j) = j ! Mode number
   allowedrvmodes(2,j) = 0 ! Counter for number of times this mode has been spotted
  END DO
  DO i=1,nrlen
    DO j=0,9
      IF( rvmode(i) .EQ. allowedrvmodes(1,j) ) THEN
        allowedrvmodes(2,j) = allowedrvmodes(2,j) + 1
      END IF
    END DO
  END DO
  ! We need to know how many Vsys's we require
  Vsyslen = 0
  DO j=0,9
    IF( allowedrvmodes(2,j) .GT. 0 ) THEN
      Vsyslen = Vsyslen + 1
    END IF
  END DO
  IF( Vsyslen .LT. 10 ) THEN
    write(*,*) 'There are ',Vsyslen,' RV modes'
  ELSE
    write(*,*) 'CRITICAL ERROR! JAMMI CAN ONLY HANDLE <10 RV MODES (Vsyslen = ',&
    Vsyslen,')'
  END IF
  
END SUBROUTINE iVsys1
! =======================================================

! =======================================================
SUBROUTINE iVsys2(rvmode,rv,sigrv,allowedrvmodes,nrlen,Vsyslen,&
           Rpar_Vsys,Rflag_Vsys,Rmin_Vsys,Rmax_Vsys,Rsol_Vsys,Rdel_Vsys)
! covariance computes the covariance matrix of Z

implicit none

 INTEGER, DIMENSION(nrlen), INTENT(INOUT) :: rvmode
 REAL(8), DIMENSION(nrlen), INTENT(IN) :: rv, sigrv
 INTEGER, INTENT(IN) :: nrlen, Vsyslen
 INTEGER, DIMENSION(2,10), INTENT(IN) :: allowedrvmodes
 INTEGER :: i, j, m
 CHARACTER(LEN=9), DIMENSION(Vsyslen), INTENT(OUT) :: Rpar_Vsys
 INTEGER, DIMENSION(Vsyslen), INTENT(OUT) :: Rflag_Vsys
 REAL(8), DIMENSION(Vsyslen), INTENT(OUT) :: Rmin_Vsys, Rmax_Vsys, Rsol_Vsys, Rdel_Vsys
 CHARACTER(LEN=9) :: Vsyschar
 CHARACTER(LEN=1) :: numeral1

  ! Now allocate a Rpar_Vsys array
  ! For each mode with >0 instances, we define a unique Vsys
  Vsyschar = 'Vsys____ '
  ! Go through and assign Rpar_Vsys values
  m = 0
  DO j=0,9
    IF( allowedrvmodes(2,j) .GT. 0 ) THEN
       m = m + 1
       write( numeral1, '(i1)' ) j ! Let numeral = "j"
       Rpar_Vsys(m) = TRIM(Vsyschar)//numeral1
       write(*,*) 'Creating parameter ',Rpar_Vsys(m)
    END IF
  END DO
  ! 1. Here we assign Rmin,Rmax & Rflag values for the Vsys terms
  ! 2. Then, we assign Rsol & Rdel values for the Vsys terms
  ! 3. Finally, we rename the mode numbers to a sequential order (for radial.f90)
  ! 
  ! 1. First, Rflag, Rmin and Rmax (easy ones)
  DO j=1,Vsyslen
    Rflag_Vsys(j) = 2
    Rmin_Vsys(j) = -1.0D9 ! Greater than speed of light, so should be OK!
    Rmax_Vsys(j) = 1.0D9  ! Greater than speed of light, so should be OK!
  END DO
  ! 2. Next, Rsol&Rdel: we assume to be the mean of all RV points in the specified mode
  ! Go through and assign Rpar_Vsys values
  m = 0
  DO j=0,9 ! Look through 10 allowed modes...
    IF( allowedrvmodes(2,j) .GT. 0 ) THEN ! => We have found a valid mode...
       m = m + 1
       ! Find the weighted mean of this mode...(= Rsol)
       Rsol_Vsys(m) = 0.0D0
       Rdel_Vsys(m) = 0.0D0
       DO i=1,nrlen
         IF( rvmode(i) .EQ. allowedrvmodes(1,j) ) THEN
           Rsol_Vsys(m) = Rsol_Vsys(m) + rv(i)/(sigrv(i)*sigrv(i))
           Rdel_Vsys(m) = Rdel_Vsys(m) + 1.0D0/(sigrv(i)*sigrv(i))
         END IF
       END DO
       Rsol_Vsys(m) = Rsol_Vsys(m)/Rdel_Vsys(m)
       Rdel_Vsys(m) = DSQRT(1.0D0/Rdel_Vsys(m))
    END IF
  END DO
  ! 3. Finally, redub mode numbers sequentially
  m = 0
  DO j=0,9   ! Go through each allowed mode number
    IF( allowedrvmodes(2,j) .GT. 0 ) THEN ! => j^th mode is a valid mode
      m = m + 1
      DO i=1,nrlen ! Go through each data point and find those equal to this mode
        IF( rvmode(i) .EQ. allowedrvmodes(1,j) ) THEN
          rvmode(i) = m ! Redub the rvmode as m
        END IF
      END DO
    END IF
  END DO

END SUBROUTINE iVsys2
! =======================================================

! =======================================================
SUBROUTINE iOOT1(Pguess,tauguess,nplen,epochp,tp,OOTlen)

implicit none

 INTEGER, INTENT(IN) :: nplen
 INTEGER, INTENT(OUT) :: OOTlen
 REAL(8), INTENT(IN) :: Pguess, tauguess
 INTEGER :: i, j, m, k
 INTEGER, DIMENSION(nplen), INTENT(OUT) :: epochp
 REAL(8), DIMENSION(nplen), INTENT(IN) :: tp
 INTEGER, DIMENSION(:), ALLOCATABLE :: history, historytemp
 INTEGER :: nonmatches, historysize

  ! ------------------ INTELLIGENT OOT ASSIGNMENT Prt.1 -----------------------
  ! Find how many epochs exist in your data set.
  ! (This part requires decent tmid and Pdays estimates from Rsol)
  ! Rsol has now been assigned yet, but we can do a 'cheat' by using the dat
  ! array.
  ! 1] Create epoch array
  DO i=1,nplen
    epochp(i) = NINT( (tp(i) - tauguess)/Pguess )
  END DO
  ! 2] Find how many unique epochs there are
  ALLOCATE (history(1)) ! history is a working list of all identified epochs
  history(1) = epochp(1)
  write(*,*) 'Epoch ',epochp(1),' is unique'
  DO i=2,nplen
    historysize = SIZE(history)
    nonmatches = 0
    DO j=1,historysize
      IF( history(j) .NE. epochp(i) ) THEN
         nonmatches = nonmatches + 1
      END IF
    END DO
    ! Now see if the number of nonmatches indicates a new epoch...
    IF( nonmatches .GE. (historysize) ) THEN
      ! => a new epoch has been found
      write(*,*) 'Epoch ',epochp(i),' is unique!'
      ALLOCATE (historytemp(SIZE(history)))
      historytemp(:) = history(:)
      DEALLOCATE (history)
      ALLOCATE (history(SIZE(historytemp)+1))
      DO k=1,SIZE(historytemp)
         history(k) = historytemp(k)
      END DO
      history(SIZE(historytemp)+1) = epochp(i)
      DEALLOCATE (historytemp)
    END IF
  END DO
  OOTlen = SIZE(history)
  write(*,*) 'JAMMI has found ',OOTlen,' unique primary epochs...'

END SUBROUTINE iOOT1
! =======================================================

! =======================================================
SUBROUTINE iOOT2(nplen,OOTlen,sigfp,epochp,epochp_seq,&
           Rpar_OOT,Rflag_OOT,Rmin_OOT,Rmax_OOT,Rsol_OOT,Rdel_OOT)

implicit none

 INTEGER, INTENT(IN) :: nplen, OOTlen
 INTEGER, DIMENSION(nplen), INTENT(IN) :: epochp
 REAL(8), DIMENSION(nplen), INTENT(IN) :: sigfp
 INTEGER, DIMENSION(nplen), INTENT(OUT) :: epochp_seq
 INTEGER, DIMENSION(:), ALLOCATABLE :: history, historytemp
 INTEGER :: nonmatches, historysize
 INTEGER :: i, j, m, k
 CHARACTER(LEN=1) :: numeral1
 CHARACTER(LEN=2) :: numeral2
 CHARACTER(LEN=3) :: numeral3
 CHARACTER(LEN=4) :: numeral4
 CHARACTER(LEN=9), DIMENSION(OOTlen), INTENT(OUT) :: Rpar_OOT
 INTEGER, DIMENSION(OOTlen), INTENT(OUT) :: Rflag_OOT
 REAL(8), DIMENSION(OOTlen), INTENT(OUT) :: Rmin_OOT, Rmax_OOT, Rsol_OOT, Rdel_OOT

  ! 2] Find how many unique epochs there are
  ALLOCATE (history(1)) ! history is a working list of all identified epochs
  history(1) = epochp(1)
  !!write(*,*) 'Epoch ',epochp(1),' is unique'
  DO i=2,nplen
    historysize = SIZE(history)
    nonmatches = 0
    DO j=1,historysize
      IF( history(j) .NE. epochp(i) ) THEN
         nonmatches = nonmatches + 1
      END IF
    END DO
    ! Now see if the number of nonmatches indicates a new epoch...
    IF( nonmatches .GE. (historysize) ) THEN
      ! => a new epoch has been found
      !!write(*,*) 'Epoch ',epochp(i),' is unique!'
      ALLOCATE (historytemp(SIZE(history)))
      historytemp(:) = history(:)
      DEALLOCATE (history)
      ALLOCATE (history(SIZE(historytemp)+1))
      DO k=1,SIZE(historytemp)
         history(k) = historytemp(k)
      END DO
      history(SIZE(historytemp)+1) = epochp(i)
      DEALLOCATE (historytemp)
    END IF
  END DO
  ! 3] Create a unique OOT parameter for each OOTlen
  DO m=1,OOTlen
    ! Positive numbers
    IF( history(m) .GE. 0 .AND. history(m) .LE. 9 ) THEN
      write( numeral1, '(i1)' ) history(m) ! Let numeral1 = "j"
      Rpar_OOT(m) = 'OOT_p'//'0'//'0'//'0'//numeral1
      write(*,*) 'Creating parameter ',Rpar_OOT(m)
    ELSE IF( history(m) .GE. 10 .AND. history(m) .LE. 99 ) THEN
      write( numeral2, '(i2)' ) history(m) ! Let numeral2 = "j"
      Rpar_OOT(m) = 'OOT_p'//'0'//'0'//numeral2
      write(*,*) 'Creating parameter ',Rpar_OOT(m)
    ELSE IF( history(m) .GE. 100 .AND. history(m) .LE. 999 ) THEN
      write( numeral3, '(i3)' ) history(m) ! Let numeral3 = "j"
      Rpar_OOT(m) = 'OOT_p'//'0'//numeral3
      write(*,*) 'Creating parameter ',Rpar_OOT(m)
    ELSE IF( history(m) .GE. 1000 .AND. history(m) .LE. 9999 ) THEN
      write( numeral4, '(i4)' ) history(m) ! Let numeral4 = "j"
      Rpar_OOT(m) = 'OOT_p'//numeral4
      write(*,*) 'Creating parameter ',Rpar_OOT(m)
    ! Negative numbers
    ELSE IF( history(m) .LE. -1 .AND. history(m) .GE. -9 ) THEN
      write( numeral1, '(i1)' ) ABS(history(m)) ! Let numeral1 = "j"
      Rpar_OOT(m) = 'OOT_n'//'0'//'0'//'0'//numeral1
      write(*,*) 'Creating parameter ',Rpar_OOT(m)
    ELSE IF( history(m) .LE. -10 .AND. history(m) .GE. -99 ) THEN
      write( numeral2, '(i2)' ) ABS(history(m)) ! Let numeral2 = "j"
      Rpar_OOT(m) = 'OOT_n'//'0'//'0'//numeral2
      write(*,*) 'Creating parameter ',Rpar_OOT(m)
    ELSE IF( history(m) .LE. -100 .AND. history(m) .GE. -999 ) THEN
      write( numeral3, '(i3)' ) ABS(history(m)) ! Let numeral3 = "j"
      Rpar_OOT(m) = 'OOT_n'//'0'//numeral3
      write(*,*) 'Creating parameter ',Rpar_OOT(m)
    ELSE IF( history(m) .LE. -1000 .AND. history(m) .GE. -9999 ) THEN
      write( numeral4, '(i4)' ) ABS(history(m)) ! Let numeral4 = "j"
      Rpar_OOT(m) = 'OOT_n'//numeral4
      write(*,*) 'Creating parameter ',Rpar_OOT(m)
    ELSE IF( history(m) .LE. -10000 .OR. history(m) .GE. 10000 ) THEN
      write(*,*) 'CRITICAL ERROR: JAMMI CANNOT HANDLE PRI EPOCHS >10,000'
    END IF
  END DO
  ! 1. Here we assign Rmin,Rmax & Rflag values for the OOT terms
  ! 2. Then, we assign Rsol & Rdel values for the OOT terms
  ! 3. Finally, we rename the mode numbers to a sequential order (for models)
  ! 
  ! 1. First, Rflag, Rmin, Rmax & Rsol (easy ones)
  DO j=1,OOTlen
    Rflag_OOT(j) = 2
    Rmin_OOT(j) = 0.95D0
    Rmax_OOT(j) = 1.05D0
    Rsol_OOT(j) = 1.0D0
  END DO
  ! 2. Next, Rdel: we assume to be the error on the weighted mean
  DO j=1,OOTlen
    Rdel_OOT(j) = 0.0D0
    DO i=1,nplen
      IF( epochp(i) .EQ. history(j) ) THEN
        Rdel_OOT(j) = Rdel_OOT(j) + (1.0D0/(sigfp(i)**2))
      END IF
    END DO
    Rdel_OOT(j) = DSQRT( 1.0D0/Rdel_OOT(j) )
  END DO
  ! 3. Finally, overwrite history with sequential flags
  DO i=1,nplen
    DO j=1,OOTlen
      IF( epochp(i) .EQ. history(j) ) THEN
        epochp_seq(i) = j
      END IF
    END DO
  END DO

END SUBROUTINE iOOT2
! =======================================================

! =======================================================
SUBROUTINE itau2(nplen,taulen,Rpar_OOT,Pdays,tau0,Rjumptau0,&
           Rpar_tau,Rflag_tau,Rmin_tau,Rmax_tau,Rsol_tau,Rdel_tau)

implicit none

 INTEGER, INTENT(IN) :: nplen, taulen
 CHARACTER(LEN=9), DIMENSION(taulen), INTENT(IN) :: Rpar_OOT
 REAL(8), INTENT(IN) :: Pdays, tau0, Rjumptau0
 INTEGER :: i, epoch
 CHARACTER(LEN=9), DIMENSION(taulen), INTENT(OUT) :: Rpar_tau
 INTEGER, DIMENSION(taulen), INTENT(OUT) :: Rflag_tau
 REAL(8), DIMENSION(taulen), INTENT(OUT) :: Rmin_tau, Rmax_tau, Rsol_tau, Rdel_tau

  ! The purpose of itau2 is to:
  ! 1] Create OOTlen parameters for individual tau times
  ! 2] Assign a starting point, jump size, min, max and flag for each tau

  ! Step 1] Create tau parameters
  write(*,*) 'Global flag is off, JAMMI will create indiv transit times...'
  DO i=1,taulen
   ! We need to trim last 5 letters away since Rpar_OOT = 'OOT_n1234'
   Rpar_tau(i) = 'tau_'//Rpar_OOT(i)(5:9)
   write(*,*) 'Creating parameter ',Rpar_tau(i)
  END DO

  ! Step 2] Assign simple flags
  DO i=1,taulen
    Rflag_tau(i) = 2
    Rdel_tau(i) = Rjumptau0
    IF( Rpar_OOT(i)(5:5) .EQ. 'p' ) THEN
      read( Rpar_OOT(i)(6:9), '(i4)' ) epoch
      epoch = epoch + 0
    ELSE IF( Rpar_OOT(i)(5:5) .EQ. 'n' ) THEN
      read( Rpar_OOT(i)(6:9), '(i4)' ) epoch
      epoch = -epoch + 0
    END IF
    Rsol_tau(i) = tau0 + epoch*Pdays
    Rmin_tau(i) = Rsol_tau(i) - 1.0D0
    Rmax_tau(i) = Rsol_tau(i) + 1.0D0
  END DO
 
  !DO i=1,taulen
  !  write(*,*) 'Rpar = ',Rpar_tau(i)
  !  write(*,*) Rsol_tau(i),Rdel_tau(i),Rmin_tau(i),Rmax_tau(i),Rflag_tau(i)
  !END DO

END SUBROUTINE itau2
! =======================================================

! =======================================================
SUBROUTINE iOOS1(Pguess,tauguess,hguess,kguess,nslen,epochs,ts,OOSlen)

implicit none

 INTEGER, INTENT(IN) :: nslen
 INTEGER, INTENT(OUT) :: OOSlen
 REAL(8), INTENT(IN) :: Pguess, tauguess, hguess, kguess
 REAL(8) :: tsecguess
 INTEGER :: i, j, m, k
 INTEGER, DIMENSION(nslen), INTENT(OUT) :: epochs
 REAL(8), DIMENSION(nslen), INTENT(IN) :: ts
 INTEGER, DIMENSION(:), ALLOCATABLE :: history, historytemp
 INTEGER :: nonmatches, historysize
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0

  ! ------------------ INTELLIGENT OOS ASSIGNMENT Prt.1 -----------------------
  ! Find how many epochs exist in your data set.
  ! (This part requires decent tmid and Pdays estimates from Rsol)
  ! Rsol has now been assigned yet, but we can do a 'cheat' by using the dat
  ! array.
  ! 1] Create epoch array
  tsecguess = kguess*DSQRT(1.0D0-hguess**2+kguess*2)
  tsecguess = tsecguess/(1.0D0 - hguess**2)
  tsecguess = tsecguess + DATAN(kguess/DSQRT(1.0D0-hguess**2-kguess**2))
  tsecguess = (tsecguess*(Pguess/pi)) + tauguess + 0.5D0*Pguess
  DO i=1,nslen
    epochs(i) = NINT( (ts(i) - tsecguess)/Pguess )
  END DO
  ! 2] Find how many unique epochs there are
  ALLOCATE (history(1)) ! history is a working list of all identified epochs
  history(1) = epochs(1)
  write(*,*) 'Epoch ',epochs(1),' is unique'
  DO i=2,nslen
    historysize = SIZE(history)
    nonmatches = 0
    DO j=1,historysize
      IF( history(j) .NE. epochs(i) ) THEN
         nonmatches = nonmatches + 1
      END IF
    END DO
    ! Now see if the number of nonmatches indicates a new epoch...
    IF( nonmatches .GE. (historysize) ) THEN
      ! => a new epoch has been found
      write(*,*) 'Epoch ',epochs(i),' is unique!'
      ALLOCATE (historytemp(SIZE(history)))
      historytemp(:) = history(:)
      DEALLOCATE (history)
      ALLOCATE (history(SIZE(historytemp)+1))
      DO k=1,SIZE(historytemp)
         history(k) = historytemp(k)
      END DO
      history(SIZE(historytemp)+1) = epochs(i)
      DEALLOCATE (historytemp)
    END IF
  END DO
  OOSlen = SIZE(history)
  write(*,*) 'JAMMI has found ',OOSlen,' unique secondary epochs...'

END SUBROUTINE iOOS1
! =======================================================

! =======================================================
SUBROUTINE iOOS2(nslen,OOSlen,sigfs,epochs,epochs_seq,&
           Rpar_OOS,Rflag_OOS,Rmin_OOS,Rmax_OOS,Rsol_OOS,Rdel_OOS)

implicit none

 INTEGER, INTENT(IN) :: nslen, OOSlen
 INTEGER, DIMENSION(nslen), INTENT(IN) :: epochs
 REAL(8), DIMENSION(nslen), INTENT(IN) :: sigfs
 INTEGER, DIMENSION(nslen), INTENT(OUT) :: epochs_seq
 INTEGER, DIMENSION(:), ALLOCATABLE :: history, historytemp
 INTEGER :: nonmatches, historysize
 INTEGER :: i, j, m, k
 CHARACTER(LEN=1) :: numeral1
 CHARACTER(LEN=2) :: numeral2
 CHARACTER(LEN=3) :: numeral3
 CHARACTER(LEN=4) :: numeral4
 CHARACTER(LEN=9), DIMENSION(OOSlen), INTENT(OUT) :: Rpar_OOS
 INTEGER, DIMENSION(OOSlen), INTENT(OUT) :: Rflag_OOS
 REAL(8), DIMENSION(OOSlen), INTENT(OUT) :: Rmin_OOS, Rmax_OOS, Rsol_OOS, Rdel_OOS

  ! 2] Find how many unique epochs there are
  ALLOCATE (history(1)) ! history is a working list of all identified epochs
  history(1) = epochs(1)
  !!write(*,*) 'Epoch ',epochs(1),' is unique'
  DO i=2,nslen
    historysize = SIZE(history)
    nonmatches = 0
    DO j=1,historysize
      IF( history(j) .NE. epochs(i) ) THEN
         nonmatches = nonmatches + 1
      END IF
    END DO
    ! Now see if the number of nonmatches indicates a new epoch...
    IF( nonmatches .GE. (historysize) ) THEN
      ! => a new epoch has been found
      !!write(*,*) 'Epoch ',epochs(i),' is unique!'
      ALLOCATE (historytemp(SIZE(history)))
      historytemp(:) = history(:)
      DEALLOCATE (history)
      ALLOCATE (history(SIZE(historytemp)+1))
      DO k=1,SIZE(historytemp)
         history(k) = historytemp(k)
      END DO
      history(SIZE(historytemp)+1) = epochs(i)
      DEALLOCATE (historytemp)
    END IF
  END DO
  ! 3] Create a unique OOS parameter for each OOSlen
  DO m=1,OOSlen
    ! Positive numbers
    IF( history(m) .GE. 0 .AND. history(m) .LE. 9 ) THEN
      write( numeral1, '(i1)' ) history(m) ! Let numeral1 = "j"
      Rpar_OOS(m) = 'OOS_p'//'0'//'0'//'0'//numeral1
      write(*,*) 'Creating parameter ',Rpar_OOS(m)
    ELSE IF( history(m) .GE. 10 .AND. history(m) .LE. 99 ) THEN
      write( numeral2, '(i2)' ) history(m) ! Let numeral2 = "j"
      Rpar_OOS(m) = 'OOS_p'//'0'//'0'//numeral2
      write(*,*) 'Creating parameter ',Rpar_OOS(m)
    ELSE IF( history(m) .GE. 100 .AND. history(m) .LE. 999 ) THEN
      write( numeral3, '(i3)' ) history(m) ! Let numeral3 = "j"
      Rpar_OOS(m) = 'OOS_p'//'0'//numeral3
      write(*,*) 'Creating parameter ',Rpar_OOS(m)
    ELSE IF( history(m) .GE. 1000 .AND. history(m) .LE. 9999 ) THEN
      write( numeral4, '(i4)' ) history(m) ! Let numeral4 = "j"
      Rpar_OOS(m) = 'OOS_p'//numeral4
      write(*,*) 'Creating parameter ',Rpar_OOS(m)
    ! Negative numbers
    ELSE IF( history(m) .LE. -1 .AND. history(m) .GE. -9 ) THEN
      write( numeral1, '(i1)' ) ABS(history(m)) ! Let numeral1 = "j"
      Rpar_OOS(m) = 'OOS_n'//'0'//'0'//'0'//numeral1
      write(*,*) 'Creating parameter ',Rpar_OOS(m)
    ELSE IF( history(m) .LE. -10 .AND. history(m) .GE. -99 ) THEN
      write( numeral2, '(i2)' ) ABS(history(m)) ! Let numeral2 = "j"
      Rpar_OOS(m) = 'OOS_n'//'0'//'0'//numeral2
      write(*,*) 'Creating parameter ',Rpar_OOS(m)
    ELSE IF( history(m) .LE. -100 .AND. history(m) .GE. -999 ) THEN
      write( numeral3, '(i3)' ) ABS(history(m)) ! Let numeral3 = "j"
      Rpar_OOS(m) = 'OOS_n'//'0'//numeral3
      write(*,*) 'Creating parameter ',Rpar_OOS(m)
    ELSE IF( history(m) .LE. -1000 .AND. history(m) .GE. -9999 ) THEN
      write( numeral4, '(i4)' ) ABS(history(m)) ! Let numeral4 = "j"
      Rpar_OOS(m) = 'OOS_n'//numeral4
      write(*,*) 'Creating parameter ',Rpar_OOS(m)
    ELSE IF( history(m) .LE. -10000 .OR. history(m) .GE. 10000 ) THEN
      write(*,*) 'CRITICAL ERROR: JAMMI CANNOT HANDLE SEC EPOCHS >10,000'
    END IF
  END DO
  ! 1. Here we assign Rmin,Rmax & Rflag values for the OOS terms
  ! 2. Then, we assign Rsol & Rdel values for the OOS terms
  ! 3. Finally, we rename the mode numbers to a sequential order (for models)
  ! 
  ! 1. First, Rflag, Rmin, Rmax & Rsol (easy ones)
  DO j=1,OOSlen
    Rflag_OOS(j) = 2
    Rmin_OOS(j) = 0.5D0
    Rmax_OOS(j) = 1.5D0
    Rsol_OOS(j) = 1.0D0
  END DO
  ! 2. Next, Rdel: we assume to be the error on the weighted mean
  DO j=1,OOSlen
    Rdel_OOS(j) = 0.0D0
    DO i=1,nslen
      IF( epochs(i) .EQ. history(j) ) THEN
        Rdel_OOS(j) = Rdel_OOS(j) + (1.0D0/(sigfs(i)**2))
      END IF
    END DO
    Rdel_OOS(j) = DSQRT( 1.0D0/Rdel_OOS(j) )
  END DO
  ! 3. Finally, overwrite history with sequential flags
  DO i=1,nslen
    DO j=1,OOSlen
      IF( epochs(i) .EQ. history(j) ) THEN
        epochs_seq(i) = j
      END IF
    END DO
  END DO

END SUBROUTINE iOOS2
! =======================================================

END MODULE intelmod

