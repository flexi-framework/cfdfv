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

! Determine uMax and uMin
uMax = aElem%pvar
uMin = aElem%pvar
! Loop over all Sides
aSide => aElem%firstside
DO WHILE (ASSOCIATED(aSide))
!-----------------------------------------------------------------------------------------------------------------------------------
  uMax(:) = MAX(uMax, aSide%connection%elem%pvar)
  uMin(:) = MIN(uMin, aSide%connection%elem%pvar)
  aSide => aSide%nextElemSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
maxDiff(:) = uMax(:) - aElem%pvar(:)
minDiff(:) = uMin(:) - aElem%pvar(:)
maxDiff_sq(:) = maxDiff(:) * maxDiff(:)
minDiff_sq(:) = minDiff(:) * minDiff(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! Loop over all Edges: Determine phi
phi(:) = 1.
aSide => aElem%firstside
DO WHILE (ASSOCIATED(aSide))
  phiLoc(:) = 1.
  DO iVar = 1, NVAR
    uDiff(iVar) = aElem%u_x(iVar) * aSide%GP(X_DIR) + &
                  aElem%u_y(iVar) * aSide%GP(Y_DIR)
    uDiff_sq(iVar) = uDiff(iVar) * uDiff(iVar)
    IF (uDiff(iVar) > 0.) THEN
      phiLoc(iVar) = MIN(1., maxDiff(iVar) / uDiff(iVar))
    ELSEIF (uDiff(iVar) < 0.) THEN
      phiLoc(iVar) = MIN(1., minDiff(iVar) / uDiff(iVar))
    END IF
  END DO
  phi = MIN(phi, phiLoc)
  aSide => aSide%nextElemSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Compute limited Gradients
aElem%u_x = aElem%u_x * phi
aElem%u_y = aElem%u_y * phi
!-----------------------------------------------------------------------------------------------------------------------------------
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

! Determine uMax and uMin
uMax = aElem%pvar
uMin = aElem%pvar
! Loop over all Sides
aSide => aElem%firstside
DO WHILE (ASSOCIATED(aSide))
!-----------------------------------------------------------------------------------------------------------------------------------
  uMax(:) = MAX(uMax, aSide%connection%elem%pvar)
  uMin(:) = MIN(uMin, aSide%connection%elem%pvar)
  aSide => aSide%nextElemSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
maxDiff = uMax - aElem%pvar
DO iVar = 1, NVAR
  minDiff(iVar) = uMin(iVar) - aElem%pvar(iVar)
  minDiff(iVar) = SIGN(1.,minDiff(iVar)) * (ABS(minDiff(iVar)) + EPSILON(0.))
END DO
maxDiff_sq = maxDiff * maxDiff
minDiff_sq = minDiff * minDiff
!-----------------------------------------------------------------------------------------------------------------------------------
! Loop over all Edges: Determine phi
phi = 1.
aSide => aElem%firstside
DO WHILE (ASSOCIATED(aSide))
  phiLoc = 1.
  DO iVar = 1, NVAR
    uDiff(iVar) = aElem%u_x(iVar) * aSide%GP(X_DIR) + &
                  aElem%u_y(iVar) * aSide%GP(Y_DIR)
    uDiff_sq(iVar) = uDiff(iVar) * uDiff(iVar)
    IF (uDiff(iVar) > 0.) THEN
      phiLoc(iVar) = 1. / uDiff(iVar) * (((maxDiff_sq(iVar) + aElem%venk_epsilon_sq) * uDiff(iVar) + &
                                          2. * uDiff_sq(iVar) * maxDiff(iVar))                     / &
                                         (maxDiff_sq(iVar) + 2. * uDiff_sq(iVar) + uDiff(iVar)     * &
                                          maxDiff(iVar) + aElem%venk_epsilon_sq))
    ELSEIF (uDiff(iVar) < 0.) THEN
      phiLoc(iVar) = 1. / uDiff(iVar) * (((minDiff_sq(iVar) + aElem%venk_epsilon_sq) * uDiff(iVar) + &
                                          2. * uDiff_sq(iVar) * minDiff(iVar))                     / &
                                         (minDiff_sq(iVar) + 2. * uDiff_sq(iVar) + uDiff(iVar)     * &
                                          minDiff(iVar) + aElem%venk_epsilon_sq))
    END IF
  END DO
  phi = MIN(phi, phiLoc)
  aSide => aSide%nextElemSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Compute limited Gradients
aElem%u_x = aElem%u_x * phi
aElem%u_y = aElem%u_y * phi
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE Limiter_Venkatakrishnan

END MODULE MOD_Limiter
