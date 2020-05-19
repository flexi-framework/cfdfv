MODULE MOD_ExplicitTimeStep
!===================================================================================================================================
! Module to performe the explicit time step
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ExplicitTimeStepEuler
   MODULE PROCEDURE ExplicitTimeStepEuler
END INTERFACE
INTERFACE ExplicitTimeStepRK
   MODULE PROCEDURE ExplicitTimeStepRK
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: ExplicitTimeStepEuler, ExplicitTimeStepRK
!===================================================================================================================================


CONTAINS

SUBROUTINE ExplicitTimeStepEuler(time,dt,iter,res_iter)
!===================================================================================================================================
! Performing explicit time step using Euler scheme
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,          ONLY: tElem,tSide
USE MOD_Mesh_Vars,          ONLY: nElems,Elems
USE MOD_FV,                 ONLY: FV_TimeDerivative
USE MOD_EoS,                ONLY: ConsPrim
USE MOD_Analyze,            ONLY: GlobalResidual
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)      :: time,dt
INTEGER,INTENT(IN)      :: iter
REAL,INTENT(OUT)        :: res_iter(NVAR+2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                 :: iElem
TYPE(tElem), POINTER    :: aElem
!===================================================================================================================================

CALL FV_TimeDerivative(time, iter) 
!$omp parallel do private(aElem)

! Timeupdate of the conservative and primitive variables
  
  ! Insert your Code here


!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
CALL GlobalResidual(dt,res_iter)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ExplicitTimeStepEuler


SUBROUTINE ExplicitTimeStepRK(time,dt,iter,res_iter)
!===================================================================================================================================
! Performing explicit time step using a RK nstage scheme
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_TimeDisc_Vars,      ONLY: nRKStages,RKCoeff
USE MOD_Mesh_Vars,          ONLY: tElem,tSide
USE MOD_Mesh_Vars,          ONLY: nElems,Elems
USE MOD_FV,                 ONLY: FV_TimeDerivative
USE MOD_EoS,                ONLY: ConsPrim
USE MOD_Analyze,            ONLY: GlobalResidual
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)      :: time,dt
INTEGER,INTENT(IN)      :: iter
REAL,INTENT(OUT)        :: res_iter(NVAR+2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: dt_Stage
INTEGER                 :: iElem,iStage
TYPE(tElem), POINTER    :: aElem
!===================================================================================================================================
!$omp parallel do private(aElem)

! Save the initial solution as needed for the RK scheme

  ! Insert your Code here


!$omp end parallel do
! Loop over the RK stages
!-----------------------------------------------------------------------------------------------------------------------------------

  ! Insert your Code here
  
!-----------------------------------------------------------------------------------------------------------------------------------
CALL GlobalResidual(dt,res_iter)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ExplicitTimeStepRK

END MODULE MOD_ExplicitTimeStep
