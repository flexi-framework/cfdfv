MODULE MOD_Reconstruction
!===================================================================================================================================
! Reconstruction
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE SpatialReconstruction
   MODULE PROCEDURE SpatialReconstruction
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC   :: SpatialReconstruction
!===================================================================================================================================

CONTAINS


SUBROUTINE SpatialReconstruction(time)       
!===================================================================================================================================
! Computes the gradients of d U / d x
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Limiter  ,ONLY:Limiter
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:nElems,nSides,Elems,Sides
USE MOD_FV_Vars  ,ONLY:SpatialOrder
USE MOD_Boundary ,ONLY:SetBCatBarys
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)         :: time    ! time needed for BC
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: pDiff(NVAR)
REAL                    :: dx, dy
TYPE(tElem), POINTER    :: aElem
TYPE(tSide), POINTER    :: aSide
INTEGER                 :: iElem, iSide
!===================================================================================================================================

IF (SpatialOrder == 1) THEN
  !$omp parallel do private(aElem,aSide)
  
  ! Set side states to be equal to mean value
  DO iElem = 1, nElems
    aElem => Elems(iElem)%Elem
    aSide => aElem%firstSide
    aElem%u_x(:) = 0
    aElem%u_y(:) = 0
    aElem%u_t(:) = 0
    DO WHILE(ASSOCIATED(aSide))
      aSide%pVar(:) = aElem%pVar(:)
      aSide => aSide%nextElemSide
    END DO
  END DO

  !$omp end parallel do
ELSE
!-----------------------------------------------------------------------------------------------------------------------------------
! Reconstruct values at Side GPs
!-----------------------------------------------------------------------------------------------------------------------------------
! Initialize (zero out) all derivatives
  !$omp parallel do private(aElem)
  DO iElem = 1, nElems
    aElem => Elems(iElem)%Elem
    aElem%u_x(:) = 0
    aElem%u_y(:) = 0
    aElem%u_t(:) = 0
  END DO
  !$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Set Boundary Conditions
  CALL SetBCatBarys(time)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Solve Least Squares Problem for gradient
  ! (Matrix-Vector multiplication)
  DO iSide = 1, nSides
    aSide => Sides(iSide)%Side
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Solve Least Squares Problem for gradient
  ! Compute difference of primitive variables over the side
    pdiff = aSide%connection%Elem%pvar - aSide%Elem%pvar
!-----------------------------------------------------------------------------------------------------------------------------------
  ! least squares reconstruction
  ! main Elem
    aElem     => aSide%Elem
    aElem%u_x = aElem%u_x + aSide%w(1) * pDiff
    aElem%u_y = aElem%u_y + aSide%w(2) * pDiff
  ! connection Elem
    aElem     => aSide%connection%Elem
    aElem%u_x = aElem%u_x - aSide%connection%w(1) * pDiff
    aElem%u_y = aElem%u_y - aSide%connection%w(2) * pDiff
!-----------------------------------------------------------------------------------------------------------------------------------
  END DO

!-----------------------------------------------------------------------------------------------------------------------------------
!! PARALLEL VERSION OF THE ABOVE SIDES LOOP 
!  ! Solve Least Squares Problem for gradient
!  ! (Matrix-Vector multiplication)
!  !$omp parallel do private(aElem,aSide,pdiff)
!  DO iElem = 1, nElems
!    aElem => Elems(iElem)%Elem
!    aSide => aElem%firstSide
!    DO WHILE(ASSOCIATED(aSide))
!!----------------------------------------------------------------------------------------------------------------------------------
!    ! Solve Least Squares Problem for gradient
!    ! Compute difference of primitive variables over the side
!      pdiff = aSide%connection%Elem%pvar - aSide%Elem%pvar
!!----------------------------------------------------------------------------------------------------------------------------------
!    ! least squares reconstruction
!      aElem%u_x = aElem%u_x + aSide%w(1) * pDiff
!      aElem%u_y = aElem%u_y + aSide%w(2) * pDiff
!!----------------------------------------------------------------------------------------------------------------------------------
!      aSide => aSide%nextElemSide
!    END DO
!  END DO
!  !$omp end parallel do

!-----------------------------------------------------------------------------------------------------------------------------------
  !$omp parallel do private(aElem,aSide,dx,dy)

  ! Limit gradient and reconstruct values at side GPs
  DO iElem = 1, nElems
    aElem => Elems(iElem)%Elem
  ! Limit Gradient
    CALL Limiter(aElem)
  ! Reconstruct Values at side GPs
    aSide => aElem%firstSide
    DO WHILE(ASSOCIATED(aSide))
      dx = aSide%GP(X_DIR)
      dy = aSide%GP(Y_DIR)
      aSide%pvar(:) = aElem%pvar(:)     + &
                      dx * aElem%u_x(:) + &
                      dy * aElem%u_y(:)
      aSide => aSide%nextElemSide
    END DO
  END DO

  !$omp end parallel do
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE SpatialReconstruction

END MODULE MOD_Reconstruction

