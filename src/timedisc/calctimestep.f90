MODULE MOD_CalcTimeStep
!===================================================================================================================================
! Compute the time step
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CalcTimeStep
   MODULE PROCEDURE CalcTimeStep
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: CalcTimeStep
!===================================================================================================================================

CONTAINS

SUBROUTINE CalcTimeStep(printtime,dt_max,viscous_timestep)
!===================================================================================================================================
! Compute the convective and dissipative time step
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_EOS
USE MOD_Equation_Vars,ONLY: gamma,mu,Pr
USE MOD_Mesh_Vars    ,ONLY: tElem
USE MOD_Mesh_Vars    ,ONLY: nElems,Elems
USE MOD_TimeDisc_Vars,ONLY: t,CFL,DFL,StopTime,TimeStep1D
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)         :: printtime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)        :: dt_max
LOGICAL,INTENT(OUT)     :: viscous_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: csound
REAL                    :: dt_conv,dt_visc
REAL                    :: sum_spectral_radii
REAL                    :: gammasPr_max
TYPE(tElem), POINTER    :: aElem
INTEGER                 :: iElem
REAL                    :: dt_visc_max,dt_conv_max
!===================================================================================================================================

! Calculate local timestep for each cell
viscous_timestep=.FALSE.
!-----------------------------------------------------------------------------------------------------------------------------------
IF (TimeStep1D) THEN
  dt_max = 1.E150
  !$omp parallel do private(aElem,csound,dt_conv), reduction(min:dt_max)
  DO iElem = 1, nElems
    aElem => Elems(iElem)%Elem
    csound             = SQRT(gamma*aElem%pvar(P)/aElem%pvar(rho))
    dt_conv            = CFL * aElem%Sy / (abs(aElem%pvar(v1))+csound)
    IF (dt_conv.NE.dt_conv) STOP 'Convective timestep NaN'
    IF (mu.GT.1E-10) STOP 'TimeStep1D not implemented for Navier Stokes. Please set mu = 0 or turn off TimeStep1D'
    dt_max             = MIN(dt_max, dt_conv)
  END DO
  !$omp end parallel do
ELSE ! 2D Case
  gammasPr_max=MAX(4./3.,gamma/Pr)
  dt_conv_max= 1.E150
  !$omp parallel do private(aElem,csound,sum_spectral_radii,dt_conv), reduction(min:dt_conv_max)
  DO iElem = 1, nElems
    aElem => Elems(iElem)%Elem
    ! Convective time step
    csound             = SQRT(gamma*aElem%pvar(P)/aElem%pvar(rho))
    sum_spectral_radii = (abs(aElem%pvar(v1))+csound) * aElem%sx  + &
                         (abs(aElem%pvar(v2))+csound) * aElem%sy
    dt_conv            = CFL * aElem%area / sum_spectral_radii
    IF (dt_conv.NE.dt_conv) STOP 'Convective timestep NaN'
    dt_conv_max        = MIN(dt_conv_max,dt_conv)
  END DO
  !$omp end parallel do
  ! Viscous time step
  dt_visc_max= 1.E150
  IF (mu.GT.1E-10) THEN
    !$omp parallel do private(aElem,csound,sum_spectral_radii,dt_visc), reduction(min:dt_visc_max)
    DO iElem = 1, nElems
      aElem => Elems(iElem)%Elem
      sum_spectral_radii = gammasPr_max*mu*aElem%sx**2. + gammasPr_max*mu*aElem%sy**2.
      dt_visc          = DFL * (aElem%pvar(RHO))**2. * aElem%area**2./(4.*sum_spectral_radii) 
      IF (dt_visc.NE.dt_visc) STOP 'Viscous timestep NaN'
      dt_visc_max        = MIN(dt_visc_max,dt_visc)
    END DO
    !$omp end parallel do
  END IF
  dt_max = MIN(dt_conv_max, dt_visc_max)
  IF (dt_visc_max .LT. dt_conv_max) viscous_timestep=.TRUE.
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Special treatment for data output and stoptime
IF ((t+dt_max > printtime).OR.(t+dt_max > StopTime)) THEN
  dt_max = MIN(printtime, StopTime) - t
ELSEIF ((t+1.5*dt_max > printtime).OR.(t+1.5*dt_max > StopTime)) THEN
  dt_max = 0.5*(MIN(printtime, StopTime) - t)
ENDIF
!-----------------------------------------------------------------------------------------------------------------------------------
!$omp parallel do private(aElem)
! Set local timestep for each cell to the global timestep
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  aElem%dt = dt_max
END DO
! WRITE(*,*) 'dt_visc: ',dt_visc_max, ' dt_conv: ',dt_conv_max,dt_max
!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CalcTimeStep

END MODULE MOD_CalcTimeStep
