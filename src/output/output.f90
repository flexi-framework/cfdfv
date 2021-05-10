MODULE MOD_Output
!===================================================================================================================================
! Output moudels
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE
INTERFACE DataOutput
   MODULE PROCEDURE DataOutput
END INTERFACE
INTERFACE FinalizeDataOutput
   MODULE PROCEDURE FinalizeDataOutput
END INTERFACE
INTERFACE CurveOutput
   MODULE PROCEDURE CurveOutput
END INTERFACE
INTERFACE PartitionArray
   MODULE PROCEDURE PartitionArray
END INTERFACE
INTERFACE SortArray
   MODULE PROCEDURE SortArray
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: InitOutput, DataOutput, FinalizeDataOutput
!===================================================================================================================================

CONTAINS

SUBROUTINE InitOutput()
!===================================================================================================================================
! Init Output
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Readintools  ,ONLY:GETSTR,GETREAL,GETINT
USE MOD_Output_Vars  ,ONLY:strOutFile,IOTimeInterval,IOIterInterval,iVisuProg
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) 'Initializing IO:'
strOutFile     = GETSTR('FileName')
IOTimeInterval = GETREAL('IOTimeInterval')
IOIterInterval = GETREAL('IOIterInterval')
iVisuProg      = GETINT('OutputFormat','1')
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE InitOutput

SUBROUTINE DataOutput(time, iter)
!===================================================================================================================================
! Data output, swith between different formats
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Output_Vars    ,ONLY:tOutputTime,strOutFile
USE MOD_TimeDisc_Vars  ,ONLY:IO,stationary 
USE MOD_Analyze_Vars   ,ONLY: ExactSolution
USE MOD_Output_Vars    ,ONLY: iVisuProg,OutputTimes
USE MOD_Output_cgns    ,ONLY:CGNSOutput
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL   ,INTENT(IN)                 :: time
INTEGER,INTENT(IN)                 :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: time_start, time_end
CHARACTER(LEN=15)       :: cIter
CHARACTER(LEN=255)      :: FileName
TYPE(tOutputTime), POINTER::OutputTime
!$ DOUBLE PRECISION     :: omp_get_wtime
!===================================================================================================================================

! Time Measurement 
CALL CPU_TIME(time_start)
!$ time_start=omp_get_wtime()
!-----------------------------------------------------------------------------------------------------------------------------------
 ! Output Times
ALLOCATE(OutputTime)
OutputTime%Time = time
OutputTime%Iter = iter
OutputTime%Next => OutputTimes
OutputTimes => OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! Write Flow Solution
IF (stationary) THEN 
  WRITE(cIter, '(I9.9)') iter
  FileName =  TRIM(strOutFile) // '_' // TRIM(cIter)
ELSE
  FileName = TIMESTAMP(TRIM(strOutFile),time)
END IF
SELECT CASE(iVisuProg)
CASE(CGNS)
  FileName = TRIM(FileName) // '.cgns'
  CALL CGNSOutput(FileName, time, iter, .FALSE.)
CASE(CURVE)
  FileName = TRIM(FileName) // '.curve'
  CALL CurveOutput(FileName, time, iter, .FALSE.)
CASE(DAT)
  FileName = TRIM(FileName) // '.csv'
  CALL CSVOutput(FileName, time, iter, .FALSE.)
CASE DEFAULT
  WRITE(*,*) 'Error in DataOutput: Output Format unkown.'
  STOP
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
! Write Exact Solution if applicable
IF (ExactSolution) THEN
  IF (stationary) THEN 
    WRITE(cIter, '(I9.9)') iter
    FileName =  TRIM(strOutFile)//'_ex_' // TRIM(cIter)
  ELSE
    FileName = TIMESTAMP(TRIM(strOutFile)//'_ex',time)
  END IF
  SELECT CASE(iVisuProg)
  CASE(CGNS)
    FileName = TRIM(FileName) // '.cgns'
    CALL CGNSOutput(FileName, time, iter, .TRUE.)
  CASE(CURVE)
    FileName = TRIM(FileName) // '.curve'
    CALL CurveOutput(FileName, time, iter, .TRUE.)
  CASE(DAT)
    FileName = TRIM(FileName) // '.csv'
    CALL CSVOutput(FileName, time, iter, .TRUE.)
  END SELECT
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CALL CPU_TIME(time_end)
!$ time_end=omp_get_wtime()
IO = IO + time_end - time_start
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE DataOutput



SUBROUTINE FinalizeDataOutput()
!===================================================================================================================================
! Finalize Data output
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Output_Vars,ONLY:iVisuProg
USE MOD_Output_cgns,ONLY:FinalizeCGNSOutput
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

SELECT CASE(iVisuProg)
CASE(CGNS)
  CALL FinalizeCGNSOutput()
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE FinalizeDataOutput


SUBROUTINE CurveOutput(FileName, time, iter, ExactSolution)
!===================================================================================================================================
! curved data output
! 1D data
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Readin,   ONLY: GetFreeIOUnit
USE MOD_ExactFunc
USE MOD_Mesh_Vars,ONLY: FirstElem,nElems
USE MOD_Mesh_Vars,ONLY: tElem
USE MOD_Equation_Vars,ONLY: intExactFunc
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=200),INTENT(IN)      :: FileName
INTEGER,INTENT(IN)                 :: iter
REAL,INTENT(IN)                    :: time
LOGICAL,INTENT(IN)                 :: ExactSolution
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: FlowData(nElems, NVAR), tmp(NVAR)
TYPE(tElem), POINTER    :: aElem
INTEGER                 :: iElem, iUnit, allocStat, n
LOGICAL                 :: swapped
REAL                    :: pvar(NVAR)
!===================================================================================================================================

! Prepare Data (only equidistant grids are supported)
IF (ExactSolution) THEN
!-----------------------------------------------------------------------------------------------------------------------------------
! store exact solution in array
  iElem = 1
  aElem => FirstElem
  DO WHILE(ASSOCIATED(aElem))
    CALL ExactFunc(intExactFunc, aElem%Bary, time, pvar)
    FlowData(iElem, 1) = aElem%Bary(1)
    FlowData(iElem, 2) = pvar(RHO)
    FlowData(iElem, 3) = pvar(V1)
    FlowData(iElem, 4) = pvar(P)
    iElem = iElem + 1
    aElem => aElem%nextElem
  END DO
!-----------------------------------------------------------------------------------------------------------------------------------
ELSE
!-----------------------------------------------------------------------------------------------------------------------------------
! store numerical solution in array
  iElem = 1
  aElem => FirstElem
  DO WHILE(ASSOCIATED(aElem))
    FlowData(iElem, 1) = aElem%Bary(1)
    FlowData(iElem, 2) = aElem%pvar(RHO)
    FlowData(iElem, 3) = aElem%pvar(V1)
    FlowData(iElem, 4) = aElem%pvar(P)
    iElem = iElem + 1
    aElem => aElem%nextElem
  END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Sort Array
n = nElems
swapped = .TRUE.
DO WHILE(swapped.AND.(n>=1))
  swapped = .FALSE.
  DO iElem = 1, n-1
    IF (FlowData(iElem,1) > FlowData(iElem+1,1)) THEN
      tmp = FlowData(iElem+1,:)
      FlowData(iElem+1,:) = FlowData(iElem,:)
      FlowData(iElem,:) = tmp
      swapped = .TRUE.
    END IF
  END DO
  n = n-1
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Open File
CALL GetFreeIOUnit(137,iUnit)
OPEN(UNIT   = iUnit         , &
     FILE   = TRIM(FileName), &
     STATUS = 'REPLACE'     , &
     IOSTAT = allocStat       )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write Data
WRITE(iUnit,'(a)') '#Density'
DO iElem = 1, nElems
  WRITE(iUnit,'(2F15.9)') FlowData(iElem,1), FlowData(iElem,2)
END DO
WRITE(iUnit,*)
WRITE(iUnit,'(a)') '#Velocity'
DO iElem = 1, nElems
  WRITE(iUnit,'(2F15.9)') FlowData(iElem,1), FlowData(iElem,3)
END DO
WRITE(iUnit,*)
WRITE(iUnit,'(a)') '#Pressure'
DO iElem = 1, nElems
  WRITE(iUnit,'(2F15.9)') FlowData(iElem,1), FlowData(iElem,4)
END DO
WRITE(iUnit,*)
!-----------------------------------------------------------------------------------------------------------------------------------
CLOSE(iUnit)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CurveOutput

SUBROUTINE CSVOutput(FileName, time, iter, ExactSolution)
!===================================================================================================================================
! tabelled CSV output
! 1D data
!===================================================================================================================================
USE MOD_Globals
USE MOD_ExactFunc
USE MOD_Readin       ,ONLY: GetFreeIOUnit
USE MOD_Mesh_Vars    ,ONLY: FirstElem,nElems
USE MOD_Mesh_Vars    ,ONLY: tElem
USE MOD_Equation_Vars,ONLY: intExactFunc,gamma
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=200),INTENT(IN)      :: FileName
INTEGER,INTENT(IN)                 :: iter
REAL,INTENT(IN)                    :: time
LOGICAL,INTENT(IN)                 :: ExactSolution
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: FlowData(nElems, NVAR+1), tmp(NVAR+1)
TYPE(tElem), POINTER    :: aElem
INTEGER                 :: iElem, iUnit, allocStat, n
LOGICAL                 :: swapped
REAL                    :: pvar(NVAR)
!===================================================================================================================================

! Local Variables
! Prepare Data (only equidistant grids are supported)
IF (ExactSolution) THEN
!-----------------------------------------------------------------------------------------------------------------------------------
! store exact solution in array
  iElem = 1
  aElem => FirstElem
  DO WHILE(ASSOCIATED(aElem))
    CALL ExactFunc(intExactFunc, aElem%Bary, time, pvar)
    FlowData(iElem, 1) = aElem%Bary(1)
    FlowData(iElem, 2) = pvar(RHO)
    FlowData(iElem, 3) = pvar(V1)
    FlowData(iElem, 4) = pvar(P)
    FlowData(iElem, 5) = 0.             ! change dummy output for exact solution here
    iElem = iElem + 1
    aElem => aElem%nextElem
  END DO
!-----------------------------------------------------------------------------------------------------------------------------------
ELSE
!-----------------------------------------------------------------------------------------------------------------------------------
! store numerical solution in array
  iElem = 1
  aElem => FirstElem
  DO WHILE(ASSOCIATED(aElem))
    FlowData(iElem, 1) = aElem%Bary(1)
    FlowData(iElem, 2) = aElem%pvar(RHO)
    FlowData(iElem, 3) = aElem%pvar(V1)
    FlowData(iElem, 4) = aElem%pvar(P)
    FlowData(iElem, 5) = 0.             ! change dummy output for numerical solution here
    iElem = iElem + 1
    aElem => aElem%nextElem
  END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Sort Array
n = nElems
swapped = .TRUE.
DO WHILE(swapped.AND.(n>=1))
  swapped = .FALSE.
  DO iElem = 1, n-1
    IF (FlowData(iElem,1) > FlowData(iElem+1,1)) THEN
      tmp = FlowData(iElem+1,:)
      FlowData(iElem+1,:) = FlowData(iElem,:)
      FlowData(iElem,:) = tmp
      swapped = .TRUE.
    END IF
  END DO
  n = n-1
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Open File
CALL GetFreeIOUnit(137,iUnit)
OPEN(UNIT   = iUnit         , &
     FILE   = TRIM(FileName), &
     STATUS = 'REPLACE'     , &
     IOSTAT = allocStat       )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write Data
WRITE(iUnit,'(A12)',ADVANCE='NO') 'CoordinateX'
WRITE(iUnit,'(A1)',ADVANCE='NO') ','
WRITE(iUnit,'(A7)',ADVANCE='NO') 'Density'
WRITE(iUnit,'(A1)',ADVANCE='NO') ','
WRITE(iUnit,'(A8)',ADVANCE='NO') 'Velocity'
WRITE(iUnit,'(A1)',ADVANCE='NO') ','
WRITE(iUnit,'(A8)',ADVANCE='NO') 'Pressure'
WRITE(iUnit,'(A1)',ADVANCE='NO') ','
WRITE(iUnit,'(A11)',ADVANCE='NO') 'DummyOutput' ! remember to change string length A... when changing the string
WRITE(iunit,'(A1)') ' ' 
DO iElem = 1, nElems
  WRITE(iUnit,104) FlowData(iElem,1), ',', FlowData(iElem,2), ',', FlowData(iElem,3), ',', FlowData(iElem,4), ',', FlowData(iElem,5)
END DO
104    FORMAT (4(F15.9,A1),F15.9)
   WRITE(iUnit,*)
!-----------------------------------------------------------------------------------------------------------------------------------
   CLOSE(iUnit)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CSVOutput


RECURSIVE SUBROUTINE SortArray(FlowData)
!===================================================================================================================================
! sort array
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT /OUTPUT VARIABLES
REAL, DIMENSION(:,:) :: FlowData
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
  INTEGER              :: iq
!===================================================================================================================================

WRITE(*,*) 'Size: ', size(FlowData,1)
IF(size(FlowData) > 1) THEN
  CALL PartitionArray(FlowData, iq)
  CALL SortArray (FlowData(:iq-1,:))
  CALL SortArray (FlowData(iq:,:))
ENDIF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE SortArray
                                                                                                  !

SUBROUTINE PartitionArray(A, marker)
!===================================================================================================================================
! Partition Array
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, DIMENSION(:,:) :: A
INTEGER              :: marker
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER   :: i, j
REAL      :: temp(4)
REAL      :: x              ! pivot point
!===================================================================================================================================

  x = A(1,1)

  i = 0
  j = size(A,1) + 1

STOP
  DO
    j = j-1
    DO
WRITE(*,*) j
       IF (A(j,1) <= x) EXIT
       j = j-1
     END DO
     i = i+1
     DO
       IF (A(i,1) >= x) EXIT
       i = i+1
     END DO
     IF (i < j) THEN
       ! exchange A(i) and A(j)
       temp   = A(i,:)
       A(i,:) = A(j,:)
       A(j,:) = temp
     ELSEIF (i == j) THEN
       marker = i+1
       RETURN
     ELSE
       marker = i
       RETURN
     ENDIF
   END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE PartitionArray

FUNCTION TIMESTAMP(Filename,Time)
!===================================================================================================================================
! Creates a timestamp, consistent of a filename (project name + processor) and current time niveau
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: Filename  ! (file)name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=255)             :: TimeStamp ! the complete timestamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                        :: i         ! loop variable
REAL                           :: Time      ! time
!===================================================================================================================================

WRITE(TimeStamp,'(F15.7)')Time
! Replace spaces with 0's
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i).EQ.' ') TimeStamp(i:i)='0'
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)
!-----------------------------------------------------------------------------------------------------------------------------------
END FUNCTION TIMESTAMP

END MODULE MOD_Output
