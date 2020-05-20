MODULE MOD_Readin
!===================================================================================================================================
! Module to read in all parameter
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CGNS_ReadSolution
   MODULE PROCEDURE CGNS_ReadSolution
END INTERFACE
INTERFACE getCmdLine
   MODULE PROCEDURE getCmdLine
END INTERFACE
INTERFACE getFreeIOUnit
   MODULE PROCEDURE getFreeIOUnit
END INTERFACE
INTERFACE ucase
   MODULE PROCEDURE ucase
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: CGNS_ReadSolution, getCmdLine, getFreeIOUnit, ucase
!===================================================================================================================================

CONTAINS


SUBROUTINE CGNS_ReadSolution()
!===================================================================================================================================
! read solution from CGNS file (restart)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_EoS      ,ONLY:PrimCons
USE MOD_Mesh_Vars,ONLY:tElem
USE MOD_Mesh_Vars,ONLY:nElems,FirstElem
USE MOD_Mesh_Vars,ONLY:strIniCondFile
USE MOD_TimeDisc_Vars,ONLY:t,Time_Overall
!-----------------------------------------------------------------------------------------------------------------------------------
! Include CGNS Library:
! (Please note that the CGNS library has to be installed in the computer's
! library and include path (see CGNS documentation for more information:
! www.cgns.org)
USE CGNS
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: ierr
INTEGER              :: CGNSUnit
INTEGER              :: iSize(1,3)
CHARACTER(LEN=32)    :: ZoneName
CHARACTER(LEN=32)    :: DescriptorName
CHARACTER(LEN=256)   :: text
REAL,ALLOCATABLE     :: pressure(:), density(:), x_velocity(:), y_velocity(:)
TYPE(tElem), POINTER :: aElem
!===================================================================================================================================

! Open CGNS file for reading
CALL cg_open_f(TRIM(strIniCondFile), CG_MODE_READ, CGNSUnit, ierr)
! Get number of elements
CALL cg_zone_read_f(CGNSUnit, 1, 1, ZoneName, iSize, ierr)
! Check if the number of elements in the CGNS file corresponds with the
! number if elements in the MESH
IF (nElems /= iSize(1,2)) THEN
  WRITE(*,*) '  ERROR reading CGNS flow solution: wrong number of elements.'
  STOP
END IF
! Allocate array for the flow solution
ALLOCATE(pressure(nElems), density(nElems), x_velocity(nElems), y_velocity(nElems))
! Read Flow solution
CALL cg_field_read_f(CGNSUnit, 1, 1, 1, 'Density', RealDouble, 1, nElems , density(:), ierr)
CALL cg_field_read_f(CGNSUnit, 1, 1, 1, 'VelocityX', RealDouble, 1, nElems , x_velocity(:), ierr)
CALL cg_field_read_f(CGNSUnit, 1, 1, 1, 'VelocityY', RealDouble, 1, nElems , y_velocity(:), ierr)
CALL cg_field_read_f(CGNSUnit, 1, 1, 1, 'Pressure', RealDouble, 1, nElems , pressure(:), ierr)
! Read iteration number, time and wall clock time
CALL cg_goto_f(CGNSUnit, 1, ierr, 'end')
CALL cg_descriptor_read_f(1, DescriptorName, text, ierr)
READ(text, '(F20.12,1X,F20.12)') t, Time_Overall
! Close CGNS File
CALL cg_close_f(CGNSUnit, ierr)
! Save CGNS solution into MESH
aElem => FirstElem
DO WHILE (ASSOCIATED(aElem))
  aElem%pvar(RHO) = density(aElem%ID)
  aElem%pvar(V1)  = x_velocity(aElem%ID)
  aElem%pvar(V2)  = y_velocity(aElem%ID)
  aElem%pvar(P)   = pressure(aElem%ID)
  aElem => aElem%NextElem
END DO
DEALLOCATE(pressure, density, x_velocity, y_velocity)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CGNS_ReadSolution


SUBROUTINE getCmdLine(Unit, Line)
!===================================================================================================================================
! Get next commandline
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: Unit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT)    :: Line
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=256)              :: RawLine
CHARACTER(LEN=256)              :: CmdLine
LOGICAL                         :: Cmd_isEmpty
INTEGER                         :: Cmnt_Start
!===================================================================================================================================

Cmd_isEmpty = .TRUE.
DO WHILE(Cmd_isEmpty) !Skip empty lines
   Read(Unit, '(a)') RawLine
   RawLine = TRIM(ADJUSTL(RawLine))
   Cmnt_Start = INDEX(RawLine, '!')
   SELECT CASE (Cmnt_Start)
      CASE (0)              !No Comment
         CmdLine = TRIM(RawLine)
      CASE (1)              !whole line is Comment
         CmdLine = ''       !empty CmdLine
      CASE DEFAULT
         CmdLine = TRIM(RawLine(:Cmnt_Start-1))
   END SELECT
   Cmd_isEmpty = (LEN_TRIM(ADJUSTL(CmdLine)).EQ.0)
END DO
Line = TRIM(ADJUSTL(CmdLine))
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE getCmdLine

SUBROUTINE ucase(string)
!===================================================================================================================================
! shift string to uppercase
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT & OUTPUT VARIABLES
CHARACTER(LEN=*)    :: string
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i, length
!===================================================================================================================================

length = LEN(string)

DO i = 1, length
  IF (LGE(string(i:i), 'a') .AND. LLE(string(i:i), 'z')) THEN
    string(i:i) = ACHAR(IACHAR(string(i:i)) - 32)
  END IF
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ucase

SUBROUTINE readErrors(Unit)
!===================================================================================================================================
! subroutine to read errors
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars ,ONLY:ExactSolution
USE MOD_Equation_Vars,ONLY:intExactFunc
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: Unit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(256)      :: actLine
!===================================================================================================================================

ExactSolution = .FALSE.
REWIND(UNIT=Unit)
actLine = ''
DO WHILE((TRIM(actLine).ne.'ERRORS:').and.(TRIM(actLine).ne.'END'))
   CALL getCmdLine(Unit, actLine)
END DO
IF (TRIM(actLine).NE.'END') THEN
  ExactSolution = .TRUE.
  WRITE(*,*)
  WRITE(*,*) '-[Error Computation]----------------------------------------'
  WRITE(*,'(a, I4)') '  Errors are computed using exact function # ', intExactFunc
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE readErrors

SUBROUTINE GetFreeIOUnit(Start,IO_Unit)
!===================================================================================================================================
! Open a new file
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: Start
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT) :: IO_Unit
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: fopened
INTEGER             :: TempUnit
!===================================================================================================================================

TempUnit = Start

INQUIRE(UNIT = TempUnit, OPENED = fopened)
IF (fopened) THEN
DO WHILE(fopened)
   PRINT*, TempUnit, ' besetzt!'
   TempUnit = TempUnit + 1
   INQUIRE(UNIT = TempUnit, OPENED = fopened)
END DO
END IF
IO_Unit = TempUnit
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE GetFreeIOUnit

END MODULE MOD_Readin
