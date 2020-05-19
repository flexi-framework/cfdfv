MODULE MOD_FluxCalculation
!===================================================================================================================================
! FluxCalculation main routine
! select the used numerical flux
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ConvectiveFlux
   MODULE PROCEDURE ConvectiveFlux
END INTERFACE
INTERFACE FluxCalculation
   MODULE PROCEDURE FluxCalculation
END INTERFACE
INTERFACE FluxJacobianFD
   MODULE PROCEDURE FluxJacobianFD
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PRIVATE  :: ConvectiveFlux
PUBLIC   :: FluxCalculation
PUBLIC  :: FluxJacobianFD
!===================================================================================================================================

CONTAINS

SUBROUTINE FluxCalculation()
!===================================================================================================================================
! Calculation of left and right state. The velocity vector is transforemd into the normal system of the cell interfaces.
! finishes with backrotation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tSide
USE MOD_Mesh_Vars,ONLY:nSides,Sides
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: flux_local(NVAR)
REAL                    :: pvar(NVAR)
REAL                    :: pvar_l(NVAR), pvar_r(NVAR)
INTEGER                 :: iSide
TYPE(tSide), POINTER    :: aSide
REAL                    :: flux_conv(NVAR),flux_diffX(NVAR),flux_diffY(NVAR)
REAL                    :: state_mean(NVAR)
REAL                    :: grad_Uxmean(NVAR),grad_Uymean(NVAR)
REAL                    :: grad_Ux(NVAR),grad_Uy(NVAR)
REAL                    :: correction(NVAR)
REAL                    :: BaryBary(2) ! scaled vector of connection from 
                                       ! cell barycenter to cell barycenter
!===================================================================================================================================
!$omp parallel do private(aSide,pvar_l,pvar_r,pvar,&
!$omp& state_mean,grad_Uxmean,grad_Uymean,BaryBary,correction,grad_Ux,grad_Uy,&
!$omp& flux_conv,flux_diffX,flux_diffY,flux_local) 

! Loop over all sides
DO iSide = 1, nSides
  aSide => Sides(iSide)%Side
! Extract left state 
  pvar(:) = aSide%pvar(:)
! rotate it into normal direction
  pvar_l(RHO) = pvar(RHO)
  pvar_l(V1)  = aSide%n(1) * pvar(V1) + aSide%n(2) * pvar(V2)
  pvar_l(V2)  = - aSide%n(2) * pvar(V1) + aSide%n(1) * pvar(V2)
  pvar_l(P)   = pvar(P)
! Extract right state:
  pvar(:) = aSide%connection%pvar(:)
! rotate state into normal direction
  pvar_r(RHO) = pvar(RHO)
  pvar_r(V1)  = aSide%n(1) * pvar(V1) + aSide%n(2) * pvar(V2)
  pvar_r(V2)  = - aSide%n(2) * pvar(V1) + aSide%n(1) * pvar(V2)
  pvar_r(P)   = pvar(P)
!-----------------------------------------------------------------------------------------------------------------------------------
! Extract left and right gradient

  ! Insert your Code here
  
!-----------------------------------------------------------------------------------------------------------------------------------
! Flux calculation
  ! advection flux
  CALL ConvectiveFlux(pvar_l(RHO), pvar_r(RHO), &
                      pvar_l(V1),  pvar_r(V1),  &
                      pvar_l(V2),  pvar_r(V2),  &
                      pvar_l(P),   pvar_r(P),   &
                      flux_conv                )
  ! diffusion flux (in different axis direction)
  CALL DiffusionFlux(state_mean,                  &
                     grad_Ux,grad_Uy,             &
                     flux_diffX,flux_diffY        )
  flux_local(:) = flux_conv
!-----------------------------------------------------------------------------------------------------------------------------------
! Rotate Flux into global coordinate system and update residual
  aSide%flux(RHO) = flux_local(RHO)
  aSide%flux(M1)  = aSide%n(1) *flux_local(M1) - &
                    aSide%n(2) *flux_local(M2)
  aSide%flux(M2)  = aSide%n(2) *flux_local(M1) + &
                    aSide%n(1) *flux_local(M2)
  aSide%flux(E)   = flux_local(E)
! Sum up diffusion part of the flux
  aSide%flux(:)  = aSide%flux(:) - flux_diffX(:)*aSide%n(1) - flux_diffY(:)*aSide%n(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! Integrate flux over edge using the midpoint rule
  aSide%flux = aSide%flux * aSide%length
!-----------------------------------------------------------------------------------------------------------------------------------
! Set flux of connection cell
  aSide%connection%flux = - aSide%flux
!-----------------------------------------------------------------------------------------------------------------------------------
END DO

!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE FluxCalculation


SUBROUTINE FluxJacobianFD(t,ElemID,iVar)
!===================================================================================================================================
! Calculation of left and right state. The velocity vector is transformed into the normal system of the cell interfaces.
! Finishes with backrotation
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY:nSides,Sides
USE MOD_Mesh_Vars,              ONLY:Elems,tElem,tSide,nElems
USE MOD_LinearSolver_Vars,      ONLY:R_XK,rEps0,srEps0
USE MOD_Boundary,               ONLY:Boundary
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: ElemID,iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)         :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem), POINTER    :: aElem,gElem
TYPE(tSide), POINTER    :: aSide,gSide
REAL                    :: flux_local(NVAR),flux_diff(NVAR)
REAL                    :: pvar(NVAR)
REAL                    :: pvar_l(NVAR), pvar_r(NVAR),pvar_dummy(NVAR),n(2)
INTEGER                 :: NBElemID
!===================================================================================================================================


! get element
aElem => Elems(ElemID)%Elem
aElem%u_t(:) = 0.
! loop over all sides
aSide => aElem%FirstSide
DO WHILE(ASSOCIATED(aSide))
  ! Extract left state 
  pvar(:) = aElem%pvar(:)
  pvar(iVar) = aElem%pvar(iVar)+rEps0
  ! Extract right state:
  IF(ASSOCIATED(aSide%connection%BC))THEN
    ! boundary condition
    gSide =>aSide%connection
    gElem =>gSide%Elem
    CALL Boundary(gSide,                  &
                  t,                      &
                  pvar,                   &
                  pvar_dummy,             &
                  gElem%Bary              )
    n=aSide%n
  ELSE ! no BC SIDE
    pvar_dummy(:) = aSide%connection%Elem%pvar(:)
    n=aSide%n
  END IF
  ! rotate elment state into normal direction
  pvar_l(RHO) = pvar(RHO)
  pvar_l(V1)  = n(1) * pvar(V1) + n(2) * pvar(V2)
  pvar_l(V2)  = - n(2) * pvar(V1) + n(1) * pvar(V2)
  pvar_l(P)   = pvar(P)
  ! rotate neighbor state into normal direction
  pvar_r(RHO) = pvar_dummy(RHO)
  pvar_r(V1)  = n(1) * pvar_dummy(V1) + n(2) * pvar_dummy(V2)
  pvar_r(V2)  = - n(2) * pvar_dummy(V1) + n(1) * pvar_dummy(V2)
  pvar_r(P)   = pvar_dummy(P)
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Flux calculation
  CALL ConvectiveFlux(pvar_l(RHO), pvar_r(RHO), &
                      pvar_l(V1),  pvar_r(V1),  &
                      pvar_l(V2),  pvar_r(V2),  &
                      pvar_l(P),   pvar_r(P),   &
                      flux_local                )
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Rotate Flux into global coordinate system, build differenz to old flux 
  ! and Integrate flux over edge using the midpoint rule
  flux_diff(RHO) = -aSide%flux(RHO) + aSide%length*flux_local(RHO)
  flux_diff(M1)  = -aSide%flux(M1 ) + aSide%length*(n(1) *flux_local(M1) - n(2) *flux_local(M2))
  flux_diff(M2)  = -aSide%flux(M2 ) + aSide%length*(n(2) *flux_local(M1) + n(1) *flux_local(M2))
  flux_diff(E)   = -aSide%flux(E  ) + aSide%length*flux_local(E)

  aElem%u_t = aElem%u_t + flux_diff
!-----------------------------------------------------------------------------------------------------------------------------------
  ! add contribution to neibor element. caution change sign for flux AND  for u_t, therefore, sign is not changed
  NBElemID=aSide%connection%Elem%ID
  IF((NBElemID.GT.0).AND.(NBElemID).LE.nElems)THEN
    aSide%connection%Elem%u_t =aSide%connection%Elem%Area_q *flux_diff*srEps0
  END IF
  aSide => aSide%nextElemSide
END DO

aElem%u_t(:) = - aElem%u_t(:)*aElem%Area_q*srEps0

END SUBROUTINE FluxJacobianFD


SUBROUTINE ConvectiveFlux(rho_l, rho_r, v1_l, v1_r,           & 
                          v2_l, v2_r, p_l, p_r, flux_local    )
!===================================================================================================================================
! Select the convective flux
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Equation_Vars, ONLY: iFlux
USE MOD_flux_godunov
USE MOD_Flux_Roe
USE MOD_Flux_HLL
USE MOD_Flux_HLLE
USE MOD_Flux_HLLC
USE MOD_Flux_lax_friedrichs
USE MOD_Flux_Steger_Warming
USE MOD_Flux_Central
USE MOD_Flux_AUSMD
USE MOD_Flux_AUSMDV
USE MOD_Flux_van_leer
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                     :: rho_l, rho_r, v1_l, v1_r
REAL,INTENT(IN)                     :: v2_l, v2_r, p_l, p_r
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                    :: flux_local(NVAR)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

SELECT CASE (iFlux)
  CASE (GOD)
    CALL flux_godunov( rho_l, rho_r, &
                       v1_l, v1_r,   &
                       v2_l, v2_r,   &
                       p_l, p_r,     &
                       flux_local    )
  CASE (ROE)
    CALL flux_roe( rho_l, rho_r, &
                   v1_l, v1_r,   &
                   v2_l, v2_r,   &
                   p_l, p_r,     &
                   flux_local    )
  CASE (HLL)
    CALL flux_hll(  rho_l, rho_r, &
                    v1_l, v1_r,   &
                    v2_l, v2_r,   &
                    p_l, p_r,     &
                    flux_local    )
  CASE (HLLE)
    CALL flux_hlle( rho_l, rho_r, &
                    v1_l, v1_r,   &
                    v2_l, v2_r,   &
                    p_l, p_r,     &
                    flux_local    )
  CASE (HLLC)
    CALL flux_hllc( rho_l, rho_r, &
                    v1_l, v1_r,   &
                    v2_l, v2_r,   &
                    p_l, p_r,     &
                    flux_local    )
  CASE (LXF)
    CALL flux_lax_friedrichs( rho_l, rho_r, &
                              v1_l, v1_r,   &
                              v2_l, v2_r,   &
                              p_l, p_r,     &
                              flux_local    )
  CASE (STW)
    CALL flux_steger_warming( rho_l, rho_r, &
                              v1_l, v1_r,   &
                              v2_l, v2_r,   &
                              p_l, p_r,     &
                              flux_local    )
  CASE (CEN)
    CALL flux_central( rho_l, rho_r, &
                       v1_l, v1_r,   &
                       v2_l, v2_r,   &
                       p_l, p_r,     &
                       flux_local    )
  CASE (AUSMD)
    CALL flux_ausmd( rho_l, rho_r, &
                     v1_l, v1_r,   &
                     v2_l, v2_r,   &
                     p_l, p_r,     &
                     flux_local    )
  CASE (AUSMDV)
    CALL flux_ausmdv( rho_l, rho_r, &
                     v1_l, v1_r,   &
                     v2_l, v2_r,   &
                     p_l, p_r,     &
                     flux_local    )
  CASE (VANLEER)
    CALL flux_van_leer( rho_l, rho_r, &
                        v1_l, v1_r,   &
                        v2_l, v2_r,   &
                        p_l, p_r,     &
                        flux_local    )
  CASE DEFAULT
    write(*,*) 'ERROR: no other flux function known'
    write(*,*) '       choose  1 for Godunov scheme'
    write(*,*) '       choose  2 for Roe scheme'
    write(*,*) '       choose  3 for HLL scheme'
    write(*,*) '       choose  4 for HLLE scheme'
    write(*,*) '       choose  5 for HLLC scheme'
    write(*,*) '       choose  6 for Lax-Friedrichs scheme'
    write(*,*) '       choose  7 for Steger-Warming scheme'
    write(*,*) '       choose  8 for central scheme'
    write(*,*) '       choose  9 for AUSMD scheme'
    write(*,*) '       choose 10 for AUSMDFV scheme'
    write(*,*) '       choose 11 for van Leer scheme'
    STOP
ENDSELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ConvectiveFlux



SUBROUTINE DiffusionFlux(state,gradx,grady,f,g)
!===================================================================================================================================
! Calculate the diffusive flux
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Equation_Vars,ONLY:gamma,mu,Pr,gamma1
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)           :: state(NVAR)
REAL,INTENT(IN)           :: gradx(NVAR),grady(NVAR)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: f(NVAR),g(NVAR)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

! Insert your Code here

END SUBROUTINE DiffusionFlux

END MODULE MOD_FluxCalculation
