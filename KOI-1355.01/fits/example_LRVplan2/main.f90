PROGRAM main

	use params
	use nestwrapper
      
	implicit none
	
	INTEGER :: i, j

      	!no parameters to wrap around
      	nest_pWrap = 0

        ! OBTAIN PRIORS
        !--------------------------------------------------------------------------
        ! Open up jammi_in.txt and read file
        open(1,FILE='example_LRVplan/jammi_in.jam',FORM='FORMATTED',STATUS='UNKNOWN')
        DO j=1,sdim
          read(1,*) Rsol(j),Rdel(j),Rmin(j),Rmax(j),Rflag(j)
          write(*,*) Rsol(j),Rdel(j),Rmin(j),Rmax(j),Rflag(j)
        END DO
        close(1)

        ! OBTAIN BLENDS
        !--------------------------------------------------------------------------
        ! Open up jammi_in.txt and read file
        open(2784,FILE='example_LRVplan/blends.dat',FORM='FORMATTED',STATUS='UNKNOWN')
        DO j=1,NQ
          read(2784,*) gam(j)
					gam(j) = 1.0D0/gam(j)
        END DO
        close(2784)

        ! Mode 0 override
        DO j=1,sdim
          IF( Rflag(j) .EQ. 0 ) THEN
            Rmin(j) = Rsol(j) !- 1.0D-9
            Rmax(j) = Rsol(j) !+ 1.0D-9
          END IF
        END DO
        
        ! Write out
        DO j=1,sdim
        write(*,*) Rmin(j),Rmax(j)
        END DO
        !--------------------------------------------------------------------------

        ! Read-in the seriesP.dat
        IF( priflag .EQ. 1 ) THEN
          write(*,*) 'Reading in primary data...'
          open(unit=30301,file='example_LRVplan/seriesP.jam')
          DO i = 1,nplen
            read(30301,*) tp(i),fp(i),sigfp(i),epochp_seq(i)
            fpwei(i) = fp(i)/(sigfp(i)**2)
            sigfpwei(i) = 1.0D0/(sigfp(i)**2)
            logsigfp(i) = DLOG(sigfp(i)**2)
            ! tpdeviant(i) = time from Rsol ephemeris
            tpdeviant(i) = tp(i) - Rsol(5) &
                           - Rsol(4)*DNINT( (tp(i)-Rsol(5)) / Rsol(4) )
          END DO
          close(30301)
        END IF

        ! Read-in the seriesR.dat
        IF( rvflag .EQ. 1 ) THEN
          write(*,*) 'Reading in rv data...'
          open(unit=30302,file='example_LRVplan/seriesR.dat')
          DO i = 1,nrlen
            read(30302,*) tr(i),rv(i),sigrv(i),rvmode(i)
          END DO
          close(30302)
        END IF

      	call nest_Sample
END
