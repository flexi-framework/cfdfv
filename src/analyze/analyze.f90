MODULE MOD_Analyze
!===================================================================================================================================
! Analyze Module
! Contains calculation of error norms, C_L & C_D, residual
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyze
   MODULE PROCEDURE InitAnalyze
END INTERFACE
INTERFACE Analyze
   MODULE PROCEDURE Analyze
END INTERFACE
INTERFACE Calc_Coef
   MODULE PROCEDURE Calc_Coef
END INTERFACE
INTERFACE CalcErrors
   MODULE PROCEDURE CalcErrors
END INTERFACE
INTERFACE IniWing
   MODULE PROCEDURE IniWing
END INTERFACE
INTERFACE IniRecordPoints
   MODULE PROCEDURE IniRecordPoints
END INTERFACE
INTERFACE EvalRecordPoints
   MODULE PROCEDURE EvalRecordPoints
END INTERFACE
INTERFACE GlobalResidual
   MODULE PROCEDURE GlobalResidual
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: Calc_Coef, CalcErrors, IniWing, IniRecordPoints , EvalRecordPoints, &
             Analyze, InitAnalyze, GlobalResidual
!===================================================================================================================================

CONTAINS

SUBROUTINE InitAnalyze()
!===================================================================================================================================
! Init Analyze
!===================================================================================================================================
! MODULES
USE MOD_Readintools,  ONLY:GETLOGICAL
USE MOD_Analyze_Vars, ONLY: CalcWing,RECORDPOINT,ExactSolution
USE MOD_Readin ,      ONLY: GetFreeIOUnit
USE MOD_TimeDisc_Vars,ONLY: stationary,restart,IniIterationNumber,AbortVarName,AbortVariable
USE MOD_Output_Vars,  ONLY: strOutFile,ResUnit
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                    :: allocStat, i
CHARACTER(LEN=85), POINTER :: ResFile(:)
CHARACTER(LEN=256)         :: ResFileName, demFileName
INTEGER                    :: FileUnit
LOGICAL                    :: FileExist
!===================================================================================================================================

ExactSolution = GETLOGICAL('ExactSolution','F')
CalcWing      = GETLOGICAL('CalcWing','F')

! Open File for residuals
IF (CalcWing.OR.stationary) THEN
  CALL getFreeIOUnit(110,ResUnit)
  ResFileName = TRIM(strOutFile)//'_analysis.csv'
  IF (Restart) THEN
    INQUIRE(FILE=ResFileName,EXIST=FileExist)
    IF(FileExist)THEN
      ALLOCATE(ResFile(1:IniIterationNumber))
      OPEN(UNIT   = ResUnit       , &
           FILE   = ResFileName   , &
           STATUS = 'OLD'         , &
           IOSTAT = allocStat       )
      DO i = 1, IniIterationNumber
        READ(ResUnit,'(a)') ResFile(i)
      END DO
      CLOSE(ResUnit)
      OPEN(UNIT   = ResUnit       , &
           FILE   = ResFileName   , &
           STATUS = 'OLD'         , &
           IOSTAT = allocStat       )
      IF (CalcWing) THEN
        WRITE(ResUnit,'(A7)',ADVANCE='NO') 'iter'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A13)',ADVANCE='NO') 'time'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A9,A3,A3)',ADVANCE='NO') 'residual(',AbortVarName,')'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'C_L'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'C_D'
        WRITE(ResUnit,'(A1)') ' '
      ELSE
        WRITE(ResUnit,'(A7)',ADVANCE='NO') 'iter'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A13)',ADVANCE='NO') 'time'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(rho)'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(v1)'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(v2)'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(E)'
        WRITE(ResUnit,'(A1)') ' '
      END IF
      DO i = 1, IniIterationNumber
        WRITE(ResUnit,'(a)') ResFile(i)
      END DO
      DEALLOCATE(ResFile)
    ELSE ! no residual file exists
      ! create header
      OPEN(UNIT   = ResUnit       , &
           FILE   = ResFileName   , &
           STATUS = 'NEW'     , &
           IOSTAT = allocStat       )
      IF (CalcWing) THEN
        WRITE(ResUnit,'(A7)',ADVANCE='NO') 'iter'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A13)',ADVANCE='NO') 'time'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A9,A3,A3)',ADVANCE='NO') 'residual(',AbortVarName,')'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'C_L'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'C_D'
        WRITE(ResUnit,'(A1)') ' '
      ELSE
        WRITE(ResUnit,'(A7)',ADVANCE='NO') 'iter'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A13)',ADVANCE='NO') 'time'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(rho)'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(v1)'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(v2)'
        WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
        WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(E)'
        WRITE(ResUnit,'(A1)') ' '
      END IF
    END IF
  ELSE ! new calculation
    OPEN(UNIT   = ResUnit       , &
         FILE   = ResFileName   , &
         STATUS = 'REPLACE'     , &
         IOSTAT = allocStat       )
    IF (CalcWing) THEN
      WRITE(ResUnit,'(A7)',ADVANCE='NO') 'iter'
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A13)',ADVANCE='NO') 'time'
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A9,A3,A3)',ADVANCE='NO') 'residual(',AbortVarName,')  '
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A15)',ADVANCE='NO') 'C_L'
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A15)',ADVANCE='NO') 'C_D'
      WRITE(ResUnit,'(A1)') ' '
    ELSE
      WRITE(ResUnit,'(A7)',ADVANCE='NO') 'iter'
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A13)',ADVANCE='NO') 'time'
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(rho)'
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(v1)'
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(v2)'
      WRITE(ResUnit,'(A1)',ADVANCE='NO') ','
      WRITE(ResUnit,'(A15)',ADVANCE='NO') 'residual(E)'
      WRITE(ResUnit,'(A1)') ' '
    END IF
  END IF
END IF
! Gnuplot file for residuals
IF (CalcWing) THEN
  CALL ReadWing()
  CALL getFreeIOUnit(120,FileUnit)
  demFileName = TRIM(strOutFile)//'_residuals.dem'
  OPEN(UNIT   = FileUnit      , &
       FILE   = demFileName   , &
       STATUS = 'REPLACE'     , &
       IOSTAT = allocStat       )
  WRITE(FileUnit, *) 'set title "Residual plot"'
  WRITE(FileUnit, *) 'set logscale y'
  WRITE(FileUnit, *) 'plot "', TRIM(resFileName),'" using 1:3 title "',AbortVarName,'" with lines'
  WRITE(FileUnit, *) 'pause -1'
  CLOSE(FileUnit)
  CALL getFreeIOUnit(120,FileUnit)
  demFileName = TRIM(strOutFile)//'_ClCd.dem'
  OPEN(UNIT   = FileUnit      , &
       FILE   = demFileName   , &
       STATUS = 'REPLACE'     , &
       IOSTAT = allocStat       )
  WRITE(FileUnit, *) 'set title "CL/CD plot"'
  WRITE(FileUnit, *) 'plot "', TRIM(resFileName),'" using 1:4 title "CL" with lines,\'
  WRITE(FileUnit, *) '"', TRIM(resFileName),'" using 1:5 title "CD" with lines'
  WRITE(FileUnit, *) 'pause -1'
  CLOSE(FileUnit)
ELSE
  IF (stationary) THEN
    CALL getFreeIOUnit(120,FileUnit)
    demFileName = TRIM(strOutFile)//'_residuals.dem'
    OPEN(UNIT   = FileUnit      , &
         FILE   = demFileName   , &
         STATUS = 'REPLACE'     , &
         IOSTAT = allocStat       )
    WRITE(FileUnit, *) 'set title "Residual plot"'
    WRITE(FileUnit, *) 'set logscale y'
    WRITE(FileUnit, *) 'plot "', TRIM(resFileName),'" using 1:3 title "rho" with lines,\'
    WRITE(FileUnit, *) '"', TRIM(resFileName),'" using 1:4 title "m1" with lines,\'
    WRITE(FileUnit, *) '"', TRIM(resFileName),'" using 1:5 title "m2" with lines,\'
    WRITE(FileUnit, *) '"', TRIM(resFileName),'" using 1:6 title "e" with lines'
    WRITE(FileUnit, *) 'pause -1'
    CLOSE(FileUnit)
  END IF
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Initialize Ca, Cw and Cp calculation
IF (CalcWing) THEN
  CALL IniWing()
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Initialize Recordpoints
IF (RECORDPOINT%nPoints > 0) THEN
  CALL IniRecordPoints()
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE InitAnalyze


SUBROUTINE Analyze(time, iter, res_iter)
!===================================================================================================================================
! Computes aerodynamic coefficients and extracts values at record points
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars   ,ONLY: stationary
USE MOD_Analyze_Vars    ,ONLY: WING,RECORDPOINT,CalcWing
USE MOD_Mesh_Vars       ,ONLY: FirstElem
USE MOD_Output_Vars     ,ONLY: ResUnit
USE MOD_TimeDisc_Vars   ,ONLY: time_overall,AbortVariable
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: time
INTEGER,INTENT(IN)        :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)        :: res_iter(NVAR+2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

! Record Points
IF (RECORDPOINT%nPoints > 0) THEN
  CALL EvalRecordPoints(time)
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Aerodynamic coeficients
IF (CalcWing) THEN
  res_iter(NVAR+1) = WING%CL
  res_iter(NVAR+2) = WING%CD
  CALL Calc_Coef(iter)
  res_iter(NVAR+1) = ABS(res_iter(NVAR+1) - WING%CL)/FirstElem%dt
  res_iter(NVAR+2) = ABS(res_iter(NVAR+2) - WING%CD)/FirstElem%dt
  ! Write current residuals to disk
  WRITE( ResUnit, '(I7,A1)', ADVANCE='NO') iter,','
  WRITE( ResUnit, '(F13.8,A1,E15.4,A1,F15.10,A1,F15.10)')     &
          time + FirstElem%dt,',', res_iter(AbortVariable),',', &
          WING%CL,',', WING%CD
ELSE
  ! Write current residuals to disk if computation is stationary
  IF (stationary) THEN
    WRITE( ResUnit, '(I7,A1)', ADVANCE='NO') iter,','
    WRITE( ResUnit, '(F13.8,A1,E15.08,A1,E15.08,A1,E15.08,A1,E15.08)') &
            time + FirstElem%dt, ',',res_iter(RHO)        &
            ,',',res_iter(V1),',',res_iter(V2),',',res_Iter(E)
  END IF
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE Analyze



SUBROUTINE CalcErrors(time)
!===================================================================================================================================
! Calc errors of L1, L2, and Linf norm
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ExactFunc
USE MOD_Reconstruction
USE MOD_Mesh_Vars,ONLY:tElem,Elems,nElems
USE MOD_Mesh_Vars,ONLY:totalArea_q,FirstElem
USE MOD_Equation_Vars,ONLY:intExactFunc
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)         :: time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: L1(NVAR), L2(NVAR), Linf(NVAR)
REAL                    :: pvar_ex(NVAR), pvar(NVAR)
REAL                    :: dx(NDIM), error(NVAR)
INTEGER                 :: iGP,iElem
CHARACTER(LEN=40)       :: formatStr
TYPE(tElem), POINTER    :: aElem
!===================================================================================================================================

! null
L1   = 0.
L2   = 0.
Linf = 0.

!-----------------------------------------------------------------------------------------------------------------------------------
! Reconstruction
CALL SpatialReconstruction(time)
!-----------------------------------------------------------------------------------------------------------------------------------
!$omp parallel do private(aElem,iGP,pvar_ex,dx,pvar,error), reduction(max:Linf), reduction(+:L1,L2)
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
! Evaluate pvar and pvar_ex at all element GPs
  DO iGP = 1, aElem%nGP
  ! Exact function
    CALL ExactFunc(intExactFunc, aElem%xGP(iGP,:), time, pvar_ex(:))
  ! Space Expansion to obtain pvar at current GP
    dx = aElem%xGP(iGP,:) - aElem%Bary(:)
    pvar(1:NVAR) = aElem%pVar(1:NVAR) + dx(X_DIR) * aElem%u_x(1:NVAR) &
                                      + dx(Y_DIR) * aElem%u_y(1:NVAR)
  ! Compute errors at GP
    error = ABS(pvar_ex - pvar)
  ! Update Linf error norm
    Linf = MAX(Linf, error)
  ! L1 error
    L1 = L1 + error * aElem%wGP(iGP)
  ! L2 error
    L2 = L2 + error * error * aElem%wGP(iGP)
  END DO
END DO
!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
L1 = L1 * totalArea_q
L2 = sqrt(L2 * totalArea_q)
!-----------------------------------------------------------------------------------------------------------------------------------
! I/O
WRITE(*,*)             '------------------------------------------------------------'
WRITE(*,'(a,F10.5)')   ' Error Analysis at t = ', time
WRITE(*,'(A13,4(A16))') "","rho","v1","v2","p"
WRITE(formatStr,'(A5,I1,A7)')'(A13,',4,'ES16.7)'
WRITE(*,formatStr)' L_1       : ',L1
WRITE(*,formatStr)' L_2       : ',L2
WRITE(*,formatStr)' L_inf     : ',Linf
!WRITE(*,'(a,F10.5)')   ' Error Analysis at t = ', time
!WRITE(*,'(a)')         ' L1 Error:'
!WRITE(*,'(a, F15.10)') '   rho: ', L1(RHO)
!WRITE(*,'(a, F15.10)') '   v1 : ', L1(V1)
!WRITE(*,'(a, F15.10)') '   v2 : ', L1(V2)
!WRITE(*,'(a, F15.10)') '   p  : ', L1(P)
!WRITE(*,'(a)')         ' L2 Error:'
!WRITE(*,'(a, F15.10)') '   rho: ', L2(RHO)
!WRITE(*,'(a, F15.10)') '   v1 : ', L2(V1)
!WRITE(*,'(a, F15.10)') '   v2 : ', L2(V2)
!WRITE(*,'(a, F15.10)') '   p  : ', L2(P)
!WRITE(*,'(a)')         ' Linf Error:'
!WRITE(*,'(a, F15.10)') '   rho: ', Linf(RHO)
!WRITE(*,'(a, F15.10)') '   v1 : ', Linf(V1)
!WRITE(*,'(a, F15.10)') '   v2 : ', Linf(V2)
!WRITE(*,'(a, F15.10)') '   p  : ', Linf(P)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CalcErrors


SUBROUTINE IniWing()
!===================================================================================================================================
! Init required data for calculation of CL CD
!===================================================================================================================================
! MODULES
USE MOD_Boundary_Vars,ONLY: tBoundary
USE MOD_Mesh_Vars    ,ONLY: tPureSidePtr
USE MOD_Mesh_Vars    ,ONLY: FirstBCSide
USE MOD_Boundary_Vars,ONLY: FirstBC
USE MOD_Analyze_Vars ,ONLY: WING
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tBoundary), POINTER     :: aBC, WingBC
TYPE(tPureSidePtr), POINTER  :: aBCSidePtr, aSidePtr, aLastSidePtr, aInsSidePtr
REAL                         :: x(2), xIns(2)
LOGICAL                      :: inserted
INTEGER                      :: iSide
!===================================================================================================================================

! Get Pointer to the wing section's BC
aBC => FirstBC
DO WHILE(ASSOCIATED(aBC))
  IF (aBC%BCType * 100 + aBC%BCID == WING%Wall_ID) THEN
    WingBC => aBC
  END IF
  aBC => aBC%nextBC
END DO
! Build lists for pressure and suction side of the wing section
aBCSidePtr => FirstBCSide 
iSide=0
aBCSidePtr => FirstBCSide
DO WHILE(ASSOCIATED(aBCSidePtr))
  IF (ASSOCIATED(WingBC, aBCSidePtr%Side%BC)) THEN
    inserted = .FALSE.
    ALLOCATE(aInsSidePtr)
    aInsSidePtr%Side => aBCSidePtr%Side%connection
    xIns = aInsSidePtr%Side%GP + aInsSidePtr%Side%Elem%Bary
    NULLIFY(aInsSidePtr%NextSide)
  ! Check the direction of the side's normal vector
    IF (aBCSidePtr%Side%n(2) <= 0.) THEN
      IF (.NOT.ASSOCIATED(WING%FirstSuctionSide)) THEN
        WING%FirstSuctionSide => aInsSidePtr
        inserted = .TRUE.
      END IF
      aSidePtr => WING%FirstSuctionSide
      x = aSidePtr%Side%GP + aSidePtr%Side%Elem%Bary
      IF (xIns(1) < x(1)) THEN
        aInsSidePtr%nextSide => aSidePtr
        WING%FirstSuctionSide => aInsSidePtr
        inserted = .TRUE.
      END IF
      aSidePtr => WING%FirstSuctionSide
    ELSE
      IF (.NOT.ASSOCIATED(WING%FirstPressureSide)) THEN
        WING%FirstPressureSide => aInsSidePtr
        inserted = .TRUE.
      END IF
      aSidePtr => WING%FirstPressureSide
      x = aSidePtr%Side%GP + aSidePtr%Side%Elem%Bary
      IF (xIns(1) < x(1)) THEN
        aInsSidePtr%nextSide => aSidePtr
        WING%FirstPressureSide => aInsSidePtr
        inserted = .TRUE.
      END IF
      aSidePtr => WING%FirstPressureSide
    END IF
    aLastSidePtr => aSidePtr
    aSidePtr     => aSidePtr%nextSide
  ! Sort current side according to its coordinates into the list
    DO WHILE (ASSOCIATED(aSidePtr).AND.(inserted .EQV. .FALSE.))
      x = aSidePtr%Side%GP + aSidePtr%Side%Elem%Bary
      IF (xIns(1) < x(1)) THEN
        aLastSidePtr%nextSide => aInsSidePtr
        aInsSidePtr%nextSide => aSidePtr
        inserted = .TRUE.
      END IF
      aLastSidePtr => aSidePtr
      aSidePtr => aSidePtr%NextSide
    END DO
    IF (.NOT.inserted) THEN
      aLastSidePtr%nextSide => aInsSidePtr
    END IF
  END IF
  aBCSidePtr => aBCSidePtr%nextSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE IniWing


SUBROUTINE ReadWing()
!===================================================================================================================================
!  Get parameter
!===================================================================================================================================
! MODULES
USE MOD_Readintools ,ONLY:GETLOGICAL,GETINT,GETREAL
USE MOD_Analyze_Vars,ONLY:CalcWing,WING
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!===================================================================================================================================
IF (CalcWing) THEN
  WRITE(*,*) '-[Wing]-----------------------------------------------------'
  Wing%Referencelength = GETREAL('ReferenceLength')
  WRITE(*,'(a, F10.3)') '  Reference length of the wing section: ', Wing%Referencelength
  Wing%Wall_ID = GETINT('Wall_ID')
  WRITE(*,'(a, I4)') '  BC Code of the wing: ', Wing%Wall_ID
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ReadWing


SUBROUTINE Calc_Coef(iter)
!===================================================================================================================================
! Calc CL and CD
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Readin                ,ONLY: getFreeIOUnit
USE MOD_Mesh_Vars             ,ONLY: tPureSidePtr
USE MOD_InitialCondition_Vars ,ONLY:RefState,alpha
USE MOD_Analyze_Vars          ,ONLY:WING
USE MOD_Output_Vars           ,ONLY:strOutFile
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                        :: n(2)            ! horizontal/vertical Norm-Vetors
REAL                        :: length          ! length of the edge             
REAL                        :: p_0             ! Absolut Pressure               
REAL                        :: q_inf_q, p_inf
REAL                        :: q_inf_l_q       ! Overall Constant               
REAL                        :: v               ! initial velocity               
REAL                        :: CL              ! Lift Coefficient               
REAL                        :: CD              ! Drag Coefficient               
REAL                        :: CP              ! Pressure Coefficient           
REAL                        :: alpha_loc       ! local variable for alpha
TYPE(tPureSidePtr), POINTER :: aSidePtr        ! actual Side for cl and cd-calculation
INTEGER                     :: CPUnit1,CPUnit2, stat
! gnuplot
INTEGER                     :: demUnit
CHARACTER(LEN=128)          :: demFileName
CHARACTER(LEN=128)          :: suctionFileName
CHARACTER(LEN=128)          :: pressureFileName
LOGICAL                     :: FileExist
!===================================================================================================================================

! Initializing Values                                                       
CL = 0.
CD = 0.
v = RefState(1, V1) / (COS(alpha * ACOS(-1.) / 180.))
q_inf_q = 1. / (RefState(1, RHO)*0.5*(v*v))
q_inf_l_q = q_inf_q / WING%Referencelength
p_inf = RefState(1, P)
!-----------------------------------------------------------------------------------------------------------------------------------
! Open file for writing cp data
CALL getFreeIOUnit(130,CPUnit1)
PressureFileName= TRIM(strOutFile)//'_CP_pressureside.csv'
OPEN(UNIT   = CPUnit1                       , &
     FILE   = PressureFileName              , &
     STATUS = 'REPLACE'                     , &
     IOSTAT = stat                            )
IF (stat.NE.0) THEN
  WRITE(*,*)'ERROR: cannot open outputfile for CP I/O!'
  STOP
END IF

! Pressure side: Ca, Cw, Cp
! pressure side x, y, phi, cp_pressureside
WRITE(CPUnit1,'(A)') 'x, y, Phi, CP_pressureside'
!WRITE(CPUnit,'(A)') '#Cp: Pressure Side'
aSidePtr => WING%FirstPressureSide
DO WHILE (ASSOCIATED(aSidePtr))
  p_0    = aSidePtr%Side%pvar(P)
  n      = aSidePtr%Side%n
  length = aSidePtr%Side%length
  CL     = CL + n(2) * p_0 * length
  CD     = CD + n(1) * p_0 * length
  CP     = (p_0 - p_inf) * q_inf_q
  WRITE(CPUnit1,'(F15.7,A1,F15.7,A1,F15.7,A1,F15.7)') aSidePtr%Side%GP(1) + aSidePtr%Side%Elem%bary(1), ',', &
                                            aSidePtr%Side%GP(2) + aSidePtr%Side%Elem%bary(2),           ',', &
        ATAN2(aSidePtr%Side%GP(2)+aSidePtr%Side%Elem%bary(2),aSidePtr%Side%GP(1)+aSidePtr%Side%Elem%bary(1)),',',CP
  aSidePtr => aSidePtr%nextSide
END DO
! Open file for writing cp data
CALL getFreeIOUnit(132,CPUnit2)
suctionFileName= TRIM(strOutFile)//'_CP_suctionside.csv'
OPEN(UNIT   = CPUnit2                       , &
     FILE   = suctionFileName               , &
     STATUS = 'REPLACE'                     , &
     IOSTAT = stat                            )
IF (stat.NE.0) THEN
  WRITE(*,*)'ERROR: cannot open outputfile for CP I/O!'
  STOP
END IF
!WRITE(CPUnit,'(A)') '#Cp: Suction Side'
WRITE(CPUnit2,'(A)') 'x, y, Phi, CP_suctionside'
aSidePtr => WING%FirstSuctionSide
DO WHILE (ASSOCIATED(aSidePtr))
  p_0    = aSidePtr%Side%pvar(P)
  n      = aSidePtr%Side%n
  length = aSidePtr%Side%length
  CL     = CL + n(2) * p_0 * length
  CD     = CD + n(1) * p_0 * length
  CP     = (p_0 - p_inf) * q_inf_q
  WRITE(CPUnit2,'(F15.7,A1,F15.7,A1,F15.7,A1,F15.7)') aSidePtr%Side%GP(1) + aSidePtr%Side%Elem%bary(1), ',', &
                                            aSidePtr%Side%GP(2) + aSidePtr%Side%Elem%bary(2),           ',', &
        ATAN2(aSidePtr%Side%GP(2)+aSidePtr%Side%Elem%bary(2),aSidePtr%Side%GP(1)+aSidePtr%Side%Elem%bary(1)),',',CP
  aSidePtr => aSidePtr%nextSide
END DO
CL = CL * q_inf_l_q
CD = CD * q_inf_l_q
! Saving as global Variables
alpha_loc=alpha*ATAN(1.)/45.
Wing%CL = CL * COS(alpha_loc) - CD * SIN(alpha_loc)
Wing%CD = CD * COS(alpha_loc) + CL * SIN(alpha_loc)
! Close File
CLOSE(CPUnit1)
CLOSE(CPUnit2)

demFileName = TRIM(strOutFile)//'_CP.dem'
INQUIRE(FILE=demFileName,EXIST=FileExist)
IF(.NOT.FileExist)THEN
  CALL getFreeIOUnit(130,demUnit)
  OPEN(UNIT   = demUnit       , &
       FILE   = demFileName   , &
       STATUS = 'REPLACE'       )
  WRITE(demUnit, *) 'set title "Cp plot"'
  WRITE(demUnit, *) 'set xlabel "theta "'
  WRITE(demUnit, *) 'set ylabel "cp "'
  WRITE(demUnit, *) 'plot "', TRIM(pressureFileName),'" using 3:4 w l lc rgb "black" t "pressure side", "', &
                              TRIM(suctionFileName) ,'" using 3:4 w l lc rgb "blue" t "suction side"'
  WRITE(demUnit, *) 'pause -1'
  CLOSE(demUnit)
END IF

!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE Calc_Coef



SUBROUTINE IniRecordPoints()
!===================================================================================================================================
! Init Record Points
!===================================================================================================================================
! MODULES
USE MOD_Readin,       ONLY: getFreeIOUnit
USE MOD_Mesh_Vars,    ONLY: tElem,tSide
USE MOD_Mesh_Vars,    ONLY: FirstElem
USE MOD_Analyze_Vars, ONLY: Recordpoint
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem), POINTER    :: aElem
INTEGER                 :: iPt, iNode1, iNode2, i, FileStat
REAL                    :: u(2), x(2), projection
LOGICAL                 :: is_inside
!===================================================================================================================================

DO iPt = 1, RECORDPOINT%nPoints
  aElem => FirstElem
  DO WHILE(ASSOCIATED(aElem))
    is_inside = .TRUE.
    DO i = 1, aElem%ElemType
      iNode1 = i
      iNode2 = i+1
      IF (iNode2 > aElem%ElemType) iNode2 = 1
      u(1) = aElem%nodearray(iNode2)%Node%x(2) - aElem%nodearray(iNode1)%Node%x(2)
      u(2) = aElem%nodearray(iNode1)%Node%x(1) - aElem%nodearray(iNode2)%Node%x(1)
      u = u / sqrt(u(1)*u(1)+u(2)*u(2))

      x = RECORDPOINT%x(iPt,:) - aElem%nodearray(iNode1)%Node%x
      projection = DOT_PRODUCT(u, x)
      IF (projection > 0) THEN 
        is_inside = .FALSE.
      END IF
    END DO
    IF (is_inside) THEN
      RECORDPOINT%Elem(iPt)%Elem => aElem
      EXIT
    END IF
    aElem => aElem%nextElem
  END DO
! Error Handler
  IF (.NOT.is_inside) THEN
    WRITE(*,'(a,I2,a)') 'ERROR: Record Point # ', iPt, ' is not in Domain!'
    WRITE(*,*) 'Aborting...'
    STOP
  END IF
! Open File for writing
  CALL getFreeIOUnit(50+iPt,RECORDPOINT%FileUnit(iPt))
  OPEN(UNIT   = RECORDPOINT%FileUnit(iPt)       , &
       FILE   = TRIM(RECORDPOINT%FileName(iPt)) , &
       STATUS = 'REPLACE'                                     , &
       IOSTAT = FileStat                                        )
  WRITE(RECORDPOINT%FileUnit(iPt), '(a)') 'Time, Density, VelocityX, VelocityY, Pressure'
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE IniRecordPoints

SUBROUTINE EvalRecordPoints(time)
!===================================================================================================================================
! Evaluation of RPs
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars   ,ONLY:tElem
USE MOD_Analyze_Vars,ONLY:Recordpoint
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)         :: time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem), POINTER    :: aElem
INTEGER                 :: iPt
!===================================================================================================================================

DO iPt = 1, RECORDPOINT%nPoints
  aElem => RECORDPOINT%Elem(iPt)%Elem
  WRITE(RECORDPOINT%FileUnit(iPt), '(5(F20.12,a))') time + aElem%dt, ';'     , &
                                                         aElem%pvar(RHO), ';', &
                                                         aElem%pvar(V1), ';',  &
                                                         aElem%pvar(V2), ';',  &
                                                         aElem%pvar(P), ';'
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE EvalRecordPoints


SUBROUTINE GlobalResidual(dt,res_iter)
!===================================================================================================================================
! Calculation the global residual
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY:tElem,Elems
USE MOD_Mesh_Vars     ,ONLY:nElems,TotalArea_q
USE MOD_TimeDisc_Vars ,ONLY:implicit,CFL
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)      :: dt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: res_iter(NVAR+2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem), POINTER :: aElem
INTEGER              :: iElem
!===================================================================================================================================
res_iter = 0.0
!$omp parallel do private(aElem) reduction(+:res_iter)
DO iElem = 1, nElems
    aElem => Elems(iElem)%Elem
    ! calculate global residual
    res_iter(1:NVAR) = res_iter(1:NVAR) + aElem%Area &
                       * aElem%u_t(:) * aElem%u_t(:)
END DO
!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
! Compute 2-Norm of the residual

! Make residual independent of dt (still experimental!)
!IF(implicit)THEN
  res_iter(1:NVAR) = SQRT(res_iter(1:NVAR) * TotalArea_q)!*dt/CFL
!ELSE
!  res_iter(1:NVAR) = SQRT(res_iter(1:NVAR) * TotalArea_q)*dt
!END IF

!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE GlobalResidual

END MODULE MOD_Analyze

