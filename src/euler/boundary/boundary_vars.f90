MODULE MOD_Boundary_Vars
!===================================================================================================================================
! Boundary global vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!===================================================================================================================================

TYPE tBoundary
  INTEGER  :: BCType      ! Boundary Type:
                          ! 1: Slip Wall
                          ! 2: Navier-Stokes Wall
                          ! 3: Supersonic Inflow
                          ! 4: Supersonic Outflow
                          ! 5: Characteristic
                          ! 6: Exact Solution
  INTEGER  :: BCID        ! Sub-ID of the boundary type
! Data for Exact Solution Boundary
  INTEGER  :: ExactFunc
! Data for all BC Type with a prescribed pvar state (inflow and characteristic)
  REAL     :: pvar(NVAR)     ! Inflow state
! Data for Navier-Stokes wall BC
  REAL     :: wallVelocity   ! Wall velocity
  LOGICAL  :: adiabatic      ! Adiabatic Wall ?
  LOGICAL  :: TemperaturePrescribed ! TRUE:  Wall Temperature prescribed
                                    ! FALSE: Heat Flux prescribed
  REAL     :: Temperature    ! Temperature
  REAL     :: Heatflux       ! Heatflux
! Data for periodic BC
  REAL     :: connection(2)
! Pointer to the next BC in List
  TYPE(tBoundary), POINTER :: nextBC
END TYPE tBoundary
! Boundary conditions
INTEGER                  :: nBC
LOGICAL                  :: isPeriodic
TYPE(tBoundary), POINTER :: FirstBC                ! Head of BC list
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE MOD_Boundary_Vars
