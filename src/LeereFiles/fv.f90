MODULE MOD_FV
!===================================================================================================================================
! Contains the FV_TimeDerivative
! This subroutine calculates the residual of the FV scheme. It is the spacial discretization of the FV scheme
!===================================================================================================================================
! Modules
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
  INTERFACE InitFV
     MODULE PROCEDURE InitFV
  END INTERFACE
  INTERFACE FV_TimeDerivative
     MODULE PROCEDURE FV_TimeDerivative
  END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
  PUBLIC   :: InitFV, FV_TimeDerivative
!===================================================================================================================================

CONTAINS

SUBROUTINE InitFV()  
!===================================================================================================================================
! Init Finite Volume
!===================================================================================================================================
! Modules
USE MOD_Globals
USE MOD_Mesh_Vars             ,ONLY: tElem
USE MOD_Mesh_Vars             ,ONLY: nElems,Elems
USE MOD_Reconstruction_Vars   ,ONLY: intLimiter,venk_K
USE MOD_Readintools           ,ONLY: GETINT,GETREAL
USE MOD_FV_Vars               ,ONLY: SpatialOrder
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
INTEGER                 :: iElem
!===================================================================================================================================

WRITE(*,*)
WRITE(*,*) '-[SpaceDisc]------------------------------------------------'
! Order of spatial discretization
SpatialOrder = GETINT('SpatialOrder','1')
WRITE(*,'(a,I2)') '   Order of spatial discretization: ', SpatialOrder
IF (SpatialOrder.gt.2) THEN
   PRINT*, 'Error: Spatial discretization order must be 1 or 2'
   STOP
END IF
IF (SpatialOrder.GT. 1) THEN
  ! Limiter
  intLimiter = GETINT('Limiter','1')
  SELECT CASE (intLimiter)
  CASE(BARTHJESPERSEN)
    ! nothing to do
    WRITE(*,*) '  Limiter: Barth & Jespersen'
  CASE(VENKATAKRISHNAN)
    ! read K
    venk_k = GETREAL('venk_K','1')
    WRITE(*,*) '  Limiter: Venkatakrishnan'
    WRITE(*,'(a,F6.2)') '     Constant for Limiter: ', venk_k
    ! Initialize limiter
    DO iElem = 1, nElems
      aElem => Elems(iElem)%Elem
      aElem%venk_epsilon_sq = (venk_K * sqrt(aElem%Area)) ** 3
    END DO
  CASE DEFAULT
     WRITE(*,*) 'Error: Limiter must be 1 or 2, was: ', intLimiter
     STOP
  END SELECT
END IF

! Initialize variables
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  aElem%source = 0.
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE InitFV


SUBROUTINE FV_TimeDerivative(time, iter) 
!===================================================================================================================================
! Calling all subroutines to performe the spacial operator of the FV scheme
! Result: Residual vector
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Reconstruction
USE MOD_Source
USE MOD_EOS
USE MOD_Analyze
USE MOD_FluxCalculation       ,ONLY:FluxCalculation
USE MOD_Mesh_Vars             ,ONLY:tElem,tSide
USE MOD_Mesh_Vars             ,ONLY:nElems,Elems
USE MOD_Boundary              ,ONLY:SetBCatSides 
USE MOD_TimeDisc_Vars         ,ONLY:TimeOrder
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL   ,INTENT(IN)      :: time
INTEGER,INTENT(IN)      :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER              :: iElem
TYPE(tElem), POINTER :: aElem
TYPE(tSide), POINTER :: aSide
!===================================================================================================================================
!$omp parallel do private(aElem)
! Set dt for boundary condition calculation (saved in dt_loc)
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  aElem%dt_loc = 0.5 * aElem%dt * (TimeOrder - 1)
END DO
!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
! Compute values at side GPs
CALL SpatialReconstruction(time)
!-----------------------------------------------------------------------------------------------------------------------------------
! Set Boundary conditions
CALL SetBCatSides(time)
! Calculate the Fluxes
CALL FluxCalculation()
!-----------------------------------------------------------------------------------------------------------------------------------
! Source Terms
CALL Source(time)
!-----------------------------------------------------------------------------------------------------------------------------------
! Timeupdate of the conservative variables
!$omp parallel do private(aElem,aSide)

DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  ! Compute cell update: Loop over all element sides
  
    ! Insert your Code here
    
  ! Source term
  
    ! Insert your Code here
    
END DO

!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE FV_TimeDerivative

END MODULE MOD_FV
