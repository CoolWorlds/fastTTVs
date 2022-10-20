! Include file for example MultiNest program 'Gaussian' (see arXiv:1001.0719)

module params
implicit none


! ==============================================================================
! Model Parameters

	! sdim = Number of parameters in the model fit
      	INTEGER, PARAMETER :: sdim = 7
      	INTEGER, PARAMETER :: OOTlen = 9
      	INTEGER, PARAMETER :: taulen = 0
      	INTEGER, PARAMETER :: Vsyslen = 0
      	INTEGER, PARAMETER :: nparamorig = 7

        ! seriesP.dat controls
      	INTEGER, PARAMETER :: nplen = 972
      	INTEGER, PARAMETER :: priflag = 1
        INTEGER, PARAMETER :: NresamP = 30
        REAL(8), PARAMETER :: integ_bigP = 0.02043361111111111D0
        REAL(8), PARAMETER :: flick = 1.0D0
        INTEGER, PARAMETER :: NQ = 17			! Number of quarters
      	INTEGER, PARAMETER :: linmode = 0
        !REAL(8), PARAMETER :: theta = 0.5649093506153341D0 ! 722
        !REAL(8), PARAMETER :: costheta = 0.844637132067425D0 ! koi722
        !REAL(8), PARAMETER :: sintheta = 0.5353392523745248D0 ! koi722
        REAL(8), DIMENSION(NQ) :: gam
        REAL(8) :: beta_log

        ! seriesR.dat controls
      	INTEGER, PARAMETER :: nrlen = 12
      	INTEGER, PARAMETER :: rvflag = 0
      	INTEGER, PARAMETER :: showrv = 0
        REAL(8), PARAMETER :: pivotRV = 55784.7922045D0

        ! Control flags
        LOGICAL, PARAMETER :: globalflag = .TRUE.

 ! Quarter-to-quarter start/end points
 REAL(8), DIMENSION(NQ), PARAMETER :: QXstart = (/ 54953.52862037D0, &		! 1
                                                55002.50923614D0, &		! 2
						55093.2141461D0, &		! 3
						55185.36716152D0, &		! 4
						55276.48039581D0, &		! 5
						55372.43908262D0, &		! 6
						55463.1644707D0, &		! 7
						55568.35386174D0, &		! 8
						55641.50607599D0, &		! 9
						55739.83509146D0, &		! 10
						55834.1977742D0, &		! 11
						55932.39901141D0, &		! 12
						56015.72702732D0, &		! 13
						56107.12890281D0, &		! 14
						56206.47750326D0, &		! 15
						56305.08740097D0, &		! 16
						56390.97969746D0 /)		! 17

 REAL(8), DIMENSION(NQ), PARAMETER :: QXend = (/ 54997.99334075D0, &		! 1
						55091.47732973D0, &		! 2
						55182.50650574D0, &		! 3
						55275.21348059D0, &		! 4
						55371.17221239D0, &		! 5
						55462.30628212D0, &		! 6
						55552.55903682D0, &		! 7
						55635.35547074D0, &		! 8
						55738.93604628D0, &		! 9
						55833.27831544D0, &		! 10
						55931.33647373D0, &		! 11
						56015.03229749D0, &		! 12
						56106.0663811D0, &		! 13
						56204.33204943D0, &		! 14
						56304.14747395D0, &		! 15
						56390.96969746D0, &		! 16
						66390.96969746D0 /)		! 17

! ==============================================================================
! Model Variables

        REAL(8), DIMENSION(nplen) :: tp, fp, sigfp
        REAL(8), DIMENSION(nplen) :: fpwei, sigfpwei, logsigfp, tpdeviant
        INTEGER, DIMENSION(nplen) :: epochp_seq
        REAL(8), DIMENSION(nplen) :: fpmod

        REAL(8), DIMENSION(nrlen) :: tr, rv, sigrv
        INTEGER, DIMENSION(nrlen) :: rvmode

      	REAL(8), DIMENSION(sdim) :: Rmin, Rmax, Rsol, Rdel
        INTEGER, DIMENSION(sdim) :: Rflag

        REAL(8), PARAMETER :: radian = 0.017453292519943295D0
        REAL(8), PARAMETER :: third = 0.333333333333333333D0
        REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
        REAL(8), PARAMETER :: twopi = 6.283185307179586D0
        REAL(8), PARAMETER :: fivehalfpi = 7.853981633974483D0
        REAL(8), PARAMETER :: LogTwoPi = 1.8378770664093453D0
        REAL(8), PARAMETER :: roottwo = 1.4142135623730951D0
        REAL(8), PARAMETER :: Grv = 6.67428D-11
        REAL(8), PARAMETER :: Grvx = 5.2864092327892624D-2 ! Grvx = (86400^2*Grv)/(3pi)
      
! ==============================================================================

! Parameters for Nested Sampler
	
      	!whether to do multimodal sampling
	logical nest_mmodal 
 	parameter(nest_mmodal=.true.)
	
      	!sample with constant efficiency
	logical nest_ceff
 	parameter(nest_ceff=.false.)
	
      	!max no. of live points
      	integer nest_nlive
	parameter(nest_nlive=4000)
      
      	!tot no. of parameters, should be sdim in most cases but if you need to
      	!store some additional parameters with the actual parameters then
      	!you need to pass them through the likelihood routine
	integer nest_nPar 
	parameter(nest_nPar=sdim)
      
      	!seed for nested sampler, -ve means take it from sys clock
	integer nest_rseed 
	parameter(nest_rseed=-1)
      
      	!evidence tolerance factor
      	double precision nest_tol 
      	parameter(nest_tol=1.0)
      
      	!enlargement factor reduction parameter
      	double precision nest_efr
      	parameter(nest_efr=0.1d0)
      
      	!root for saving posterior files
      	character*100 nest_root
	parameter(nest_root='chains/LRVplan1-')
	
	!after how many iterations feedback is required & the output files should be updated
	!note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	integer nest_updInt
	parameter(nest_updInt=1000)
	
	!null evidence (set it to very high negative no. if null evidence is unknown)
	double precision nest_Ztol
	parameter(nest_Ztol=-1.d90)
      
      	!max modes expected, for memory allocation
      	integer nest_maxModes 
      	parameter(nest_maxModes=100)
      
      	!no. of parameters to cluster (for mode detection)
      	integer nest_nClsPar
      	parameter(nest_nClsPar=sdim)
      
      	!whether to resume from a previous run
      	logical nest_resume
      	parameter(nest_resume=.true.)
      
      	!whether to write output files
      	logical nest_outfile
      	parameter(nest_outfile=.true.)
      
      	!initialize MPI routines?, relevant only if compiling with MPI
	!set it to F if you want your main program to handle MPI initialization
      	logical nest_initMPI
      	parameter(nest_initMPI=.true.)
      
      	!points with loglike < nest_logZero will be ignored by MultiNest
      	double precision nest_logZero
      	parameter(nest_logZero=-huge(1d0))

      	!max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
	!has done max no. of iterations or convergence criterion (defined through nest_tol) has been satisfied
      	integer nest_maxIter
      	parameter(nest_maxIter=0)
	
	!parameters to wrap around (0 is F & non-zero T)
	integer nest_pWrap(sdim)
	
      	!feedback on the sampling progress?
      	logical nest_fb 
      	parameter(nest_fb=.true.)
!=======================================================================


end module params

! =============== COMMENTS ===============
!
! ---------------- PURPOSE OF THIS FILE ----------------
!
! Seems this is analogous to jammi_in.txt in many ways. Gives the control
! inputs for the actual code.
