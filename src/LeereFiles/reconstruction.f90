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

   ! Insert your Code here


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
  
  ! Do this side by side using a loop over all sides

       ! Insert your Code here


!-----------------------------------------------------------------------------------------------------------------------------------
  !$omp parallel do private(aElem,aSide,dx,dy)

  ! Limit gradient and reconstruct values at side GPs
  
  ! Do this element by element using a loop over all elements

       ! Insert your Code here
    


  !$omp end parallel do
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE SpatialReconstruction

END MODULE MOD_Reconstruction

