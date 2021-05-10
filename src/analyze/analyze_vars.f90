MODULE MOD_analyze_vars
!===================================================================================================================================
! Contains the required variables for analyze
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElemPtr,tBoundary,tPureSidePtr
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!===================================================================================================================================

!-----------------------------------------------------------------------------------------------------------------------------------
! Profiling
TYPE tTime
  REAL                 :: overall                 ! overall time
  REAL                 :: IO                      ! Time for I/O
  REAL                 :: Flux                    ! Time for flux calculation
  REAL                 :: SpaceRec                ! Time for spatial reconstruction
  REAL                 :: CK                      ! Time for CK procedure
  REAL                 :: TimeInt                 ! Time for temporal integration
END TYPE tTime
!-----------------------------------------------------------------------------------------------------------------------------------
! Data for aerodynamic coefficient calculation
TYPE tWing
  REAL                 :: Referencelength           ! Referencelenght of foil
  INTEGER              :: Wall_ID                   ! ID of the wing's BC
  REAL                 :: CL                        ! Lift-Coefficient
  REAL                 :: CD                        ! Drag-Coefficient
  TYPE(tBoundary), POINTER    :: WingBC             ! Pointer to the wall bc of the wing
  TYPE(tPureSidePtr), POINTER :: FirstPressureSide  ! Pointer to the first pressure side side
  TYPE(tPureSidePtr), POINTER :: FirstSuctionSide   ! Pointer to the first suction side side
END TYPE tWing
!-----------------------------------------------------------------------------------------------------------------------------------
! Record Points
TYPE tRecordPoint
  INTEGER                     :: nPoints       ! Number of Record Points
  REAL,  POINTER              :: x(:,:)        ! coordinates
  TYPE(tElemPtr), POINTER     :: Elem(:)       ! Pointer to Element containing the coordinates
  CHARACTER(LEN=255), POINTER :: FileName(:)   ! Filename for data output
  INTEGER, POINTER            :: FileUnit(:)   ! units of I/O Files
END TYPE tRecordPoint
!-----------------------------------------------------------------------------------------------------------------------------------
! Analysis of the results
LOGICAL              :: CalcWing
TYPE(tWing)          :: Wing                    ! Foil Analysis
TYPE(tRecordPoint)   :: RecordPoint
LOGICAL              :: ExactSolution
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE MOD_analyze_vars
