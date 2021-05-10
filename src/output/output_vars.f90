MODULE MOD_Output_vars
!===================================================================================================================================
! Output variables
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
! GLOBAL VARIABLES
!===================================================================================================================================
! I/O Data
INTEGER              :: uOutFile                ! unit of outputfile
INTEGER              :: uCaCwFile               ! unit of outputfile
CHARACTER(LEN=60)    :: strOutFile              ! Name of outputfile
REAL                 :: IOTimeInterval          ! time interval of output
INTEGER              :: IOIterInterval          ! time interval of output
INTEGER              :: iVisuProg               ! 
CHARACTER(LEN=256)   :: ParameterFile           ! Name of Parameter File

INTEGER              :: ResUnit
LOGICAL              :: ErrorOutput
!-----------------------------------------------------------------------------------------------------------------------------------
! Time and iterator values for data output
TYPE tOutputTime
  TYPE(tOutputTime), POINTER   :: Next            ! pointer to next list element
  REAL                         :: Time             ! time at output
  INTEGER                      :: Iter             ! iterator value at output
END TYPE tOutputTime
TYPE(tOutputTime), POINTER     ::OutputTimes
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE MOD_Output_vars 
