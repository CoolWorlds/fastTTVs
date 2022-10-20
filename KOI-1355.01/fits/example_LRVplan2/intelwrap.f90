PROGRAM intelwrap
use intelmod

implicit none

 INTEGER :: i, j
 INTEGER :: priflag, nplen, nparamorig, globalflag
 INTEGER :: OOTlen, taulen
 REAL(8), DIMENSION(:,:), ALLOCATABLE :: dat
 REAL(8), DIMENSION(:), ALLOCATABLE :: tp, fp, sigfp
 INTEGER, DIMENSION(:), ALLOCATABLE :: epochp, epochp_seq
 INTEGER, DIMENSION(:), ALLOCATABLE :: Rflag_OOT, Rflag_tau
 REAL(8), DIMENSION(:), ALLOCATABLE :: Rmin_OOT, Rmax_OOT, Rsol_OOT, Rdel_OOT
 REAL(8), DIMENSION(:), ALLOCATABLE :: Rmin_tau, Rmax_tau, Rsol_tau, Rdel_tau
 CHARACTER(LEN=9), DIMENSION(:), ALLOCATABLE :: Rpar_OOT, Rpar_tau

 nparamorig = 7

 ! Allocate dat array
 ALLOCATE (dat(5,nparamorig+2))

 ! Open up jammi_in.txt and read file
 open(11,FILE='jammi_in.txt',FORM='FORMATTED',STATUS='UNKNOWN')
 DO j=1,nparamorig+2
   read(11,*) dat(1,j),dat(2,j),dat(3,j),dat(4,j),dat(5,j)
 END DO
 close(11)

 ! Assign control information
 nplen = dat(1,nparamorig+1)   ! Data length of seriesP.dat
 priflag = dat(1,nparamorig+2) ! Pri flag
 globalflag = dat(5,nparamorig+2) ! Global flag
 write(*,*) 'Global flag = ',globalflag

 ALLOCATE (tp(nplen))
 ALLOCATE (fp(nplen))
 ALLOCATE (sigfp(nplen))
 write(*,*) 'nplen = ',nplen
 ! Read in seriesP.dat
 open(12,FILE='seriesP.dat',FORM='FORMATTED',STATUS='UNKNOWN')
 DO i=1,nplen
   read(12,*) tp(i),fp(i),sigfp(i)
 END DO
 close(12)

! iOOT
 IF( priflag .EQ. 1 ) THEN
   ! Call iOOT1
   !Pguess = dat(1,4)
   !tauguess = dat(1,12)
   ALLOCATE (epochp(nplen))
   call iOOT1(dat(1,4),dat(1,5),nplen,epochp,tp,OOTlen)
   ! Now assign and determine OOT relevant terms
   ALLOCATE (Rflag_OOT(OOTlen))
   ALLOCATE (Rmin_OOT(OOTlen))
   ALLOCATE (Rmax_OOT(OOTlen))
   ALLOCATE (Rsol_OOT(OOTlen))
   ALLOCATE (Rdel_OOT(OOTlen))
   ALLOCATE (Rpar_OOT(OOTlen))
   ALLOCATE (epochp_seq(nplen))
   call iOOT2(nplen,OOTlen,sigfp,epochp,epochp_seq,&
              Rpar_OOT,Rflag_OOT,Rmin_OOT,Rmax_OOT,Rsol_OOT,Rdel_OOT)
   IF( globalflag .EQ. 0 ) THEN
     ! Globalflag = 0 => we need individual transit times
     taulen = OOTlen
     ALLOCATE (Rflag_tau(taulen))
     ALLOCATE (Rmin_tau(taulen))
     ALLOCATE (Rmax_tau(taulen))
     ALLOCATE (Rsol_tau(taulen))
     ALLOCATE (Rdel_tau(taulen))
     ALLOCATE (Rpar_tau(taulen))
     call itau2(nplen,taulen,Rpar_OOT,dat(1,4),dat(1,5),dat(2,5),&
                Rpar_tau,Rflag_tau,Rmin_tau,Rmax_tau,Rsol_tau,Rdel_tau)
   END IF
   ! Output seriesP.jam
   open(13,FILE='seriesP.jam',FORM='FORMATTED',STATUS='UNKNOWN')
   DO i=1,nplen
     write(13,*) tp(i),fp(i),sigfp(i),epochp_seq(i)
   END DO
   close(13)
   ! Output jammi_in.jam
   open(14,FILE='jammi_in.jam',FORM='FORMATTED',STATUS='UNKNOWN')
   DO j=1,nparamorig
     write(14,*) dat(1,j),dat(2,j),dat(3,j),dat(4,j),NINT(dat(5,j))
   END DO
   !DO j=1,OOTlen
   !  write(14,*) Rsol_OOT(j),Rdel_OOT(j),Rmin_OOT(j),Rmax_OOT(j),Rflag_OOT(j)
   !END DO
   IF( globalflag .EQ. 0 ) THEN
     DO j=1,taulen
       write(14,*) Rsol_tau(j),Rdel_tau(j),Rmin_tau(j),Rmax_tau(j),Rflag_tau(j)
     END DO
   END IF
   close(14)
   IF( globalflag .EQ. 1 ) THEN
     write(*,*) 'sdim = ',nparamorig!+OOTlen
     !write(*,*) 'OOTlen = ',OOTlen
     !write(*,*) 'taulen = ',0
   ELSE
     write(*,*) 'sdim = ',nparamorig+taulen!+OOTlen
     !write(*,*) 'OOTlen = ',OOTlen
     write(*,*) 'taulen = ',taulen
   END IF
 END IF

END PROGRAM
