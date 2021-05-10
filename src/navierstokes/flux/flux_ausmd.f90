MODULE MOD_flux_ausmd
!===================================================================================================================================
! AUSMD Flux splitting
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE flux_ausmd
   MODULE PROCEDURE flux_ausmd
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: flux_ausmd
!===================================================================================================================================

CONTAINS

SUBROUTINE flux_ausmd( rho_l, rho_r, v1_l , v1_r, &
                       v2_l , v2_r , p_l  , p_r , &
                       flux_side)
!===================================================================================================================================
! AUSMD flux
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars,ONLY:gamma,gamma1q
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: rho_l, rho_r, v1_l , v1_r,v2_l , v2_r , p_l  , p_r 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: flux_side(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL             :: e_r,e_l,H_r, H_l
REAL             :: u_plus,u_minus,p_plus,p_minus
REAL             :: rho_u,cm
REAL             :: alpha_l,alpha_r
!===================================================================================================================================

! Calculate left/right energy and enthalpy
e_r=gamma1q*p_r+0.5*rho_r*(v1_r*v1_r+v2_r*v2_r)
e_l=gamma1q*p_l+0.5*rho_l*(v1_l*v1_l+v2_l*v2_l)

H_r=(e_r + p_r)/rho_r
H_l=(e_l + p_l)/rho_l

! Maximum speed of sound
cm = max(sqrt(gamma*p_r/rho_r),sqrt(gamma*p_l/rho_l))

alpha_r = 2.*p_r/rho_r/(p_l/rho_l+p_r/rho_r)
alpha_l = 2.*p_l/rho_l/(p_l/rho_l+p_r/rho_r)

IF(abs(v1_l).le.cm)THEN
    u_plus = 0.25*alpha_l*(v1_l+cm)*(v1_l+cm)/cm + &
                0.5*(1-alpha_l)*(v1_l+abs(v1_l))
    p_plus = 0.25*p_l*(v1_l+cm)*(v1_l+cm)/(cm*cm)*(2.-v1_l/cm)
ELSE
    u_plus = 0.5*(v1_l+abs(v1_l))
    p_plus = 0.5*p_l*(v1_l+abs(v1_l))/v1_l
ENDIF

IF(abs(v1_r).le.cm)THEN
    u_minus = -0.25*alpha_r*(v1_r-cm)*(v1_r-cm)/cm + &
                0.5*(1-alpha_r)*(v1_r-abs(v1_r))
    p_minus = 0.25*p_r*(v1_r-cm)*(v1_r-cm)/(cm*cm)*(2.+v1_r/cm)
ELSE
    u_minus = 0.5*(v1_r-abs(v1_r))
    p_minus = 0.5*p_r*(v1_r-abs(v1_r))/v1_r
ENDIF

! Mass flux
rho_u = u_plus*rho_l + u_minus*rho_r

flux_side(1) = rho_u
flux_side(2) = 0.5*(rho_u*(v1_r+v1_l)-abs(rho_u)*(v1_r-v1_l))+ &
                    (p_plus+p_minus)
flux_side(3) = 0.5*(rho_u*(v2_r+v2_l)-abs(rho_u)*(v2_r-v2_l))
flux_side(4) = 0.5*(rho_u*(H_r+H_l) - abs(rho_u)*(H_r-H_l))
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE flux_ausmd

END MODULE MOD_flux_ausmd