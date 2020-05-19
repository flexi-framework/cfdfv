MODULE MOD_TimeDisc
!===================================================================================================================================
! TimeDisc
! Selection of implicit or explicit time integration
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitTimedisc
   MODULE PROCEDURE InitTimedisc
END INTERFACE
INTERFACE Timedisc
   MODULE PROCEDURE Timedisc
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC   :: InitTimedisc,Timedisc
!===================================================================================================================================

CONTAINS

SUBROUTINE InitTimedisc()
!===================================================================================================================================
! Init TimDisc
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars
USE MOD_Output_Vars  ,ONLY:IOIterInterval,IOTimeInterval
USE MOD_ReadinTools  ,ONLY:GETREAL,GETINT,GETLOGICAL
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iRK
!===================================================================================================================================

WRITE(*,*)
WRITE(*,*) '-[TimeDisc]------------------------------------------------'
! Timestepping scheme (stability)
CFL = GETREAL('CFL','0.9')
DFL = GETREAL('DFL','0.9')
TimeStep1D = GETLOGICAL('TimeStep1D','.false.')
! temporal discretization
TimeOrder = GETINT('TimeOrder','1')
!WRITE(*,'(a,I2)') '   Order of temporal discretization: ', TimeOrder
IF (TimeOrder .gt. 3) THEN
  WRITE(*,*) 'Error: Temporal discretization order must be 1 or 2 ... or maybe 3'
  STOP
END IF
! time integration method
implicit=GETLOGICAL('implicit','.false.')
! Readin of the RK time stepping
nRKStages = GETINT('nRKStages','1')
RKCoeff = 0.
IF (nRKStages.GT.5) THEN
  WRITE(*,*) 'Error: No more then 5 RK stages are possible'
  STOP
END IF
IF((nRKStages.EQ.1).AND.(TimeOrder.EQ.1))THEN
  WRITE(*,'(a)') ' |   time stepping scheme: Euler time stepping'
  RKCoeff(1)   = 1.
  RKCoeff(2:5) = 0.
ELSE
  WRITE(*,'(a)') ' |   time stepping scheme: RK time stepping'
  
  ! Assign Coefficients for the Runge-Kutta scheme according to order of convergence and number of RK stages

  ! Insert your Code here

  DO iRK=1,nRKStages
    WRITE(UNIT_StdOut,'(a3,a29,I1,a3,e33.5,a3,a10)')' | ','RKCoeff ',iRK,' | ', RKCoeff(iRK),' | ',' | '
  END DO
  IF(RKCoeff(1).LE.0.)THEN
    WRITE(*,*) ' Wrong RK Coefficient!'
    STOP
  END IF
END IF
! stationary Computation
IF(Stationary)THEN
  AbortResidualCl=.FALSE.
  AbortResidualCd=.FALSE.
  AbortResidual=GETREAL('AbortResidual','1e-6')
  AbortVariable=GETINT('AbortVariable','1')
  SELECT CASE(AbortVariable)
  CASE(1)
    AbortVarName='rho'
  CASE(2)
    AbortVarName='M1'
  CASE(3)
    AbortVarName='M2'
  CASE(4)
    AbortVarName='E'
  END SELECT
  Cl_AbortResidual=GETREAL('Cl_AbortResidual','0')
  Cd_AbortResidual=GETREAL('Cd_AbortResidual','0')
  WRITE(*,*) '|  Stationary Problem'
  IF((Cl_AbortResidual.EQ.0.).AND.(Cd_AbortResidual.EQ.0.))THEN
    WRITE(*,'(A,E8.2)') '   Abort Residual: ', AbortResidual
  ELSE IF(Cl_AbortResidual.GT.0.)THEN
    WRITE(*,'(A,E8.2)') ' Cl Abort Residual: ', Cl_AbortResidual
    AbortResidualCl=.TRUE.
  ELSE IF(Cd_AbortResidual.GT.0.)THEN
    WRITE(*,'(A,E8.2)') ' Cd Abort Residual: ', Cd_AbortResidual
    AbortResidualCd=.TRUE.
  ELSE
    WRITE(*,'(A)') ' Wrong definition of abort residual'
    STOP
  END IF
ELSE
  WRITE(*,*) '|  Transient Problem'
END IF
! Maximum iteration number
MaxIter = GETINT('MaxIter','100000')
!WRITE(*,'(a,I10)') '|   Maximum Iteration Number: ', MaxIter
! Stoptime
Stoptime = GETREAL('tEnd')
!WRITE(*,'(a,E10.4)') '|   Final simulation time: ', Stoptime
!-----------------------------------------------------------------------------------------------------------------------------------
! Initialization of other variables
IF (Restart) THEN
  t = RestartTime
ELSE 
  t = 0.
END IF
IF (stationary) THEN
  WRITE(*,*) '| IniIterationNumber: ', IniIterationNumber
ELSE
  WRITE(*,*) '| Start Time: ', t
END IF
printiter = (IniIterationNumber / IOIterInterval + 1) * IOIterInterval
printtime = (INT(t / IOTimeInterval) + 1) * IOTimeInterval
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE InitTimedisc


SUBROUTINE TimeDisc()
!===================================================================================================================================
! Selection of temporal integration method
! Manegement of data output and analyze tools
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Output
USE MOD_Analyze
USE MOD_Output_Vars           ,ONLY:IOTimeInterval,IOIterInterval
USE MOD_Analyze_Vars          ,ONLY:WING,CalcWing,ExactSolution
USE MOD_TimeDisc_Vars         ,ONLY:t,MaxIter,implicit,Stationary,StopTime,IniIterationNumber,AbortResidual,StartTime, &
                                    AbortVariable,AbortVarName,AbortResidualCl,AbortResidualCd,Cl_AbortResidual,Cd_AbortResidual
USE MOD_TimeDisc_Vars         ,ONLY:PrintIter,PrintTime,nRKStages,TimeOrder
USE MOD_LinearSolver_Vars     ,ONLY:nNewtonIterGlobal,nGMRESIterGlobal
USE MOD_ExplicitTimeStep      ,ONLY:ExplicitTimeStepEuler,ExplicitTimeStepRK
USE MOD_ImplicitTimeStep      ,ONLY:ImplicitTimeStep
USE MOD_CalcTimeStep          ,ONLY:CalcTimeStep
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!REAL              :: printtime          ! Next printtime          
REAL              :: dt                 ! Timestep                
INTEGER           :: iter               ! Iteration number        
!INTEGER           :: printiter          ! next iteration for I/O  
INTEGER           :: start              ! Start-Value for Itteration
REAL              :: res_iter(NVAR + 2)
REAL              :: tStart, tEnd, tIOStart, tIOEnd
LOGICAL           :: viscous_timestep, convergenz
!$ DOUBLE PRECISION :: omp_get_wtime
!===================================================================================================================================


IF(stationary) convergenz=.FALSE.

! Write initial condition to disk
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) ' Writing Initial Condition to Disk:'
CALL DataOutput( t, IniIterationNumber)
WRITE(*,*) '      ... done.'
!-----------------------------------------------------------------------------------------------------------------------------------

! Main program loop
WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) ' Starting Computation:'
start = IniIterationNumber + 1
CALL CPU_Time(tStart)
!$ tStart=omp_get_wtime()
tIOStart = tStart
CALL CalcTimeStep(printtime,dt,viscous_timestep)
WRITE(*,*) ' Initial Time Step: ',dt
IF (viscous_timestep) WRITE(*,*) '  Viscous time step dominates!'
!-----------------------------------------------------------------------------------------------------------------------------------
! loop over all iterations
!-----------------------------------------------------------------------------------------------------------------------------------
DO iter = start, MaxIter
  ! Determine Timestep
  CALL CalcTimeStep(printtime,dt,viscous_timestep)
  ! Main computation loop
  ! Normal first order computation and RK (only one timestep needed)
  ! explizit Euler
  IF (.NOT.(implicit)) THEN
    IF ((TimeOrder == 1).AND.(nRKStages == 1)) THEN
      CALL ExplicitTimeStepEuler(t,dt,iter,res_iter)  ! Euler time integration
    ELSE
      CALL ExplicitTimeStepRK(t,dt,iter,res_iter)  ! RK time integration
    END IF
  ELSE
    CALL ImplicitTimeStep(t,dt,iter,res_iter)  ! solving a linear system with GMRES
  END IF
  t = t + dt 
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Analyze results
  !---------------------------------------------------------------------------------------------------------------------------------
  CALL Analyze(t, iter, res_iter)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! EndTime Abort Criterion
  !---------------------------------------------------------------------------------------------------------------------------------
  IF (Stoptime - t <= 1.E-15 ) THEN
    WRITE (*,*) '------------------------------------------------------------'
    WRITE (*,*) ' Time limit reached -  computation complete!'
    WRITE (*,*) ' Final time      : ', t
    WRITE (*,*) ' Iteration number: ', iter
    WRITE (*,*) ' Writing final state to disk...'
    CALL DataOutput(t, iter)
    CALL FinalizeDataOutput()
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Error Calculation (for exact functions)
    !-------------------------------------------------------------------------------------------------------------------------------
    IF (ExactSolution) THEN
      CALL CalcErrors(t)
    END IF
    EXIT
  ENDIF 
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Residual Abort Criterion
  !---------------------------------------------------------------------------------------------------------------------------------
  IF ((stationary))THEN
    IF(ABS(res_iter(AbortVariable)) .LE. AbortResidual ) THEN
      IF(AbortResidualCl)THEN
        IF(ABS(res_iter(5)).LE. Cl_AbortResidual)THEN
          WRITE (*,'(A)')      '------------------------------------------------------------'
          WRITE (*,'(A,A3,A)') ' Converged in ','C_L',' - computation complete!'
          WRITE (*,'(A,I12)')      ' Iteration number: ', iter
          WRITE (*,'(A)')      ' Writing final state to disk...'
          IF (CalcWing) THEN
            WRITE(*,'(A)')     ' Final lift and drag coefficients:'
            WRITE(*,'(A,F9.6)')     '   cl : ', WING%CL
            WRITE(*,'(A,F9.6)')     '   cd : ', WING%CD
          END IF
          convergenz=.TRUE.
        END IF
      ELSE IF(AbortResidualCd)THEN
        IF(ABS(res_iter(6)).LE. Cd_AbortResidual)THEN
          WRITE (*,'(A)')      '------------------------------------------------------------'
          WRITE (*,'(A,A3,A)') ' Converged in ','C_D',' - computation complete!'
          WRITE (*,'(A,I12)')      ' Iteration number: ', iter
          WRITE (*,'(A)')      ' Writing final state to disk...'
          IF (CalcWing) THEN
            WRITE(*,'(A)')     ' Final lift and drag coefficients:'
            WRITE(*,'(A,F9.6)')     '   cl : ', WING%CL
            WRITE(*,'(A,F9.6)')     '   cd : ', WING%CD
          END IF
          convergenz=.TRUE.
        END IF
      ELSE
        WRITE (*,'(A)')      '------------------------------------------------------------'
        WRITE (*,'(A,A3,A)') ' Converged in ',AbortVarName,' - computation complete!'
        WRITE (*,'(A,I12)')      ' Iteration number: ', iter
        WRITE (*,'(A)')      ' Writing final state to disk...'
        IF (CalcWing) THEN
          WRITE(*,'(A)')     ' Final lift and drag coefficients:'
          WRITE(*,'(A,F9.6)')     '   cl : ', WING%CL
          WRITE(*,'(A,F9.6)')     '   cd : ', WING%CD
        END IF
        convergenz=.TRUE.
      END IF
    END IF
    IF(convergenz)THEN
      CALL DataOutput(t, iter)
      CALL FinalizeDataOutput()
      !-----------------------------------------------------------------------------------------------------------------------------
      ! Error Calculation (for exact functions)
      !-----------------------------------------------------------------------------------------------------------------------------
      IF (ExactSolution) THEN
        CALL CalcErrors(t)
      END IF
      EXIT
    END IF
  ENDIF 
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Data Output
  !---------------------------------------------------------------------------------------------------------------------------------
  IF ((printtime - t <= 1.E-15).OR.(iter == printiter)) THEN
    WRITE (*,*) '------------------------------------------------------------'
    WRITE (*,*) ' Data Output at Iteration ', iter
    CALL CPU_TIME(tIOEnd)
    !$ tIOEnd=omp_get_wtime()
    WRITE (*,'(a,F7.2,a)') '  Time since last I/O: ', tIOEnd - tIOStart, ' s'
    tIOStart = tIOEnd
    IF (stationary) THEN
      IF (CalcWing) THEN
        WRITE(*,'(A13,A3,A3,E25.14)') ' Residuals:  ',AbortVarname,' : ' , res_iter(AbortVariable)
        WRITE(*,'(A19,E25.14)') '             cl  : ', res_iter(NVAR+1)
        WRITE(*,'(A19,E25.14)') '             cd  : ', res_iter(NVAR+2)
      ELSE
        WRITE(*,*) ' Residuals:  rho: ', res_iter(RHO)
        WRITE(*,*) '             m1 : ', res_iter(M1)
        WRITE(*,*) '             m2 : ', res_iter(M2)
        WRITE(*,*) '             e  : ', res_iter(E)
      END IF
    ELSE
      WRITE(*,*) ' Time     : ', t
      CALL CalcTimeStep(printtime+1.E150,dt,viscous_timestep)
      WRITE(*,*) ' TimeStep : ', dt
      IF (viscous_timestep) WRITE(*,*) '  Viscous time step dominates!'
    END IF
    IF (ExactSolution) THEN
      CALL CalcErrors(t)
    END IF
  !---------------------------------------------------------------------------------------------------------------------------------
    ! Data Output
    CALL DataOutput(t, iter)
    CALL FinalizeDataOutput()
  !---------------------------------------------------------------------------------------------------------------------------------
    IF (printtime - t <= 1.E-15) THEN
      printTime = printTime + IOTimeInterval
    END IF
    IF (iter == printiter) THEN
      printiter = printiter + IOIterInterval
    END IF
  ENDIF
ENDDO
CALL CPU_Time(tEnd)
!$ tEnd=omp_get_wtime()
!-----------------------------------------------------------------------------------------------------------------------------------
! Error Handler for the case that the maximum iteration number was reached
IF (iter .GT. MaxIter) THEN
  WRITE (*,*) 
  WRITE (*,*) '------------------------------------------------------------'
  WRITE (*,*) ' ERROR - ERROR - ERROR - ERROR - ERROR - ERROR - ERROR'
  WRITE (*,*) '------------------------------------------------------------'
  WRITE (*,*) ' Maximum iteration number reached. Calculation aborted!'
  WRITE (*,*) ' Final time ',Stoptime,' has not been reached.'
  WRITE (*,*) ' Current time: ', t
  WRITE (*,*) ' Final state will be written to disk...'
  WRITE (*,*) '------------------------------------------------------------'
  WRITE (*,*) 
  CALL DataOutput(t, iter-1)
  CALL FinalizeDataOutput()
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Error Calculation (for exact functions)
  !---------------------------------------------------------------------------------------------------------------------------------
  IF (ExactSolution) THEN
    CALL CalcErrors(t)
  END IF
ENDIF
!-----------------------------------------------------------------------------------------------------------------------------------
WRITE (*,*) '------------------------------------------------------------'
WRITE (*,'(a,F15.8,a)') ' Computation Time : ', tEnd - tStart, ' s'
IF(implicit)THEN
  WRITE(*,'(a,i10)') ' Netwon Iterations: ', nNewtonIterGlobal
  WRITE(*,'(a,i10)') ' GMRES Iterations: ', nGMRESIterGlobal
END IF
WRITE (*,*) '------------------------------------------------------------'
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE TimeDisc

END MODULE MOD_TimeDisc
