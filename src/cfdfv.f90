PROGRAM CFDFV
!===================================================================================================================================
! CFD 2D Finite Volume Code
! Control of the program. Performing the initialization.
!===================================================================================================================================
! MODULES
  USE MOD_Analyze           ,ONLY:InitAnalyze
  USE MOD_Boundary          ,ONLY:InitBoundary
  USE MOD_Equation          ,ONLY:InitEquation
  USE MOD_Initialcondition  ,ONLY:InitInitialCondition,InitialCondition
  USE MOD_Mesh              ,ONLY:InitMesh
  USE MOD_FV                ,ONLY:InitFV
  USE MOD_Output            ,ONLY:InitOutput
  USE MOD_Readin            ,ONLY:ucase,getFreeIOUnit
  USE MOD_TimeDisc          ,ONLY:InitTimeDisc,TimeDisc
  USE MOD_TimeDisc_Vars     ,ONLY:IniIterationNumber,Restart,StartTime,RestartTime,stationary
  USE MOD_Output_Vars       ,ONLY:OutputTimes,strOutFile
  USE MOD_Mesh_Vars         ,ONLY:strIniCondFile
  USE MOD_LinearSolver      ,ONLY:InitLinearSolver,FinalizeLinearSolver
  USE MOD_Readintools       ,ONLY:IgnoredStrings,GETLOGICAL
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
  CHARACTER(LEN=200)        :: RestartFile
  CHARACTER(LEN=32)         :: string
  INTEGER                   :: ParUnit,stat,nArgs,iExt1,iExt2
!===================================================================================================================================
  ! Welcome message
   WRITE(*,*) '============================================================'
   WRITE(*,*) '                          C F D F V'
   WRITE(*,*) '============================================================'
   WRITE(*,*)
   WRITE(*,*) '  Solution of the two-dimensional Euler Equations using an'
   WRITE(*,*) '             unstructured finite volume solver'
   WRITE(*,*)
   WRITE(*,*) '  This program may be used under the conditions of the GPL'
   WRITE(*,*)
   WRITE(*,*) '  Created at the Institute of Aerodynamics and Gas Dynamics'
   WRITE(*,*) '           at the University of Stuttgart, Germany'
   WRITE(*,*) '               http://www.iag.uni-stuttgart.de'
   WRITE(*,*)
   WRITE(*,*) '============================================================'
   WRITE(*,*)
!-----------------------------------------------------------------------------------------------------------------------------------
   ! Read Parameter File
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,*) ' Reading Parameter File:'
!-----------------------------------------------------------------------------------------------------------------------------------
   Stationary = GETLOGICAL('Stationary','T')
   ! Check if a restart is requested
   nArgs = COMMAND_ARGUMENT_COUNT()
   IF (nArgs .EQ. 2) THEN
     ! Set Restart Flag
     Restart = .TRUE.
     CALL GET_COMMAND_ARGUMENT(2,RestartFile)
     strIniCondFile=RestartFile
     WRITE(*,'(A,A,A)')' | Restarting from file "',TRIM(RestartFile),'":'
     ! Try to open the file
     OPEN(UNIT=ParUnit, FILE=RestartFile, STATUS = 'OLD', IOSTAT=stat)
     IF (stat /= 0 ) THEN
       WRITE(*,'(a,a,a)') ' Restart Error: No restart file ( ', TRIM(RestartFile),' ).'
       STOP
     END IF
     CLOSE(ParUnit)
     iExt1=INDEX(RestartFile,'_',BACK = .TRUE.)
     iExt2=INDEX(RestartFile,'.',BACK = .TRUE.)
     string = RestartFile(iExt1+1:iExt2-1)
     IF (stationary) THEN
       READ(string,*) IniIterationNumber
       StartTime= 0.
     ELSE
       IniIterationNumber = 0
       READ(string,*) StartTime
       RestartTime=StartTime
     END IF
   ELSE
     Restart = .FALSE.
     IniIterationNumber = 0
     StartTime = 0.
   END IF
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Calling the init routines
   CALL InitOutput()
   CALL InitEquation()
   CALL InitBoundary()
   CALL InitMesh()
   CALL InitInitialCondition()
   CALL InitFV()
   CALL InitTimeDisc()
   CALL InitLinearSolver()

   OutputTimes=> NULL();
!-----------------------------------------------------------------------------------------------------------------------------------
   ! Setting Initial Condition
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,*) ' Setting Initial Condition:'
   CALL InitialCondition()
   WRITE(*,*) ' ...done.'
!-----------------------------------------------------------------------------------------------------------------------------------
   ! Initialize Ca, Cw and Cp calculation as well as record points
   WRITE(*,*) '------------------------------------------------------------'
   WRITE(*,*) ' Initializing Analyse Module:'
   CALL InitAnalyze()
   WRITE(*,*) ' ...done.'
!-----------------------------------------------------------------------------------------------------------------------------------
   ! Print ignored strings
   CALL IgnoredStrings()
   ! Call time stepping routine
   CALL TimeDisc()
END PROGRAM CFDFV
