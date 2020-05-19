MODULE MOD_Limiter
!===================================================================================================================================
! Module containing the limiter
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE Limiter
   MODULE PROCEDURE Limiter
END INTERFACE
INTERFACE Limiter_BarthJespersen
   MODULE PROCEDURE Limiter_BarthJespersen
END INTERFACE
INTERFACE Limiter_Venkatakrishnan
   MODULE PROCEDURE Limiter_Venkatakrishnan
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC   :: Limiter
PRIVATE  :: Limiter_Venkatakrishnan, Limiter_BarthJespersen
!===================================================================================================================================

CONTAINS

SUBROUTINE Limiter(aElem)
!===================================================================================================================================
! Select Limiter
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Reconstruction_Vars,  ONLY:intLimiter
USE MOD_Mesh_Vars,            ONLY:tElem
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem), POINTER    :: aElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

! Determine uMax and uMin
SELECT CASE(intLimiter)
CASE (BARTHJESPERSEN)
  CALL Limiter_BarthJespersen(aElem)
CASE (VENKATAKRISHNAN)
  CALL Limiter_Venkatakrishnan(aElem)
CASE DEFAULT
  WRITE (*,*) ' ERROR in Limiter.f90: Limiter function unknown.'
  STOP
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE Limiter



SUBROUTINE Limiter_BarthJespersen(aElem)
!===================================================================================================================================
! Limiter after Barth & Jespersen
! 2D, unstructured Limiter
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,  ONLY:tElem,tSide
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem), POINTER    :: aElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: phi(NVAR), phiLoc(NVAR)
REAL                    :: uMax(NVAR), uMin(NVAR)
REAL                    :: MaxDiff(NVAR), MinDiff(NVAR), uDiff(NVAR)
REAL                    :: MaxDiff_sq(NVAR), MinDiff_sq(NVAR), uDiff_sq(NVAR)
INTEGER                 :: iVar
TYPE(tSide), POINTER    :: aSide
!===================================================================================================================================

! this routine computes the limited gradients aElem%u_x and aElem%u_y

        ! Insert your Code here
        
END SUBROUTINE Limiter_BarthJespersen


SUBROUTINE Limiter_Venkatakrishnan(aElem)
!===================================================================================================================================
! Venkatakrishnan limiter, additionally a limiting parameter k
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tElem,tSide
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem), POINTER    :: aElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: phi(NVAR), phiLoc(NVAR)
REAL                    :: uMax(NVAR), uMin(NVAR)
REAL                    :: MaxDiff(NVAR), MinDiff(NVAR), uDiff(NVAR)
REAL                    :: MaxDiff_sq(NVAR), MinDiff_sq(NVAR), uDiff_sq(NVAR)
INTEGER                 :: iVar
TYPE(tSide), POINTER    :: aSide
!===================================================================================================================================

! this routine computes the limited gradients aElem%u_x and aElem%u_y

    ! variables with reqired data:
    ! epsilon^2: aElem%venk_epsilon_sq

        ! Insert your Code here
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE Limiter_Venkatakrishnan

END MODULE MOD_Limiter
