MODULE MOD_flux_ausmdv
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
INTERFACE flux_ausmdv
   MODULE PROCEDURE flux_ausmdv
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: flux_ausmdv
!===================================================================================================================================

CONTAINS

SUBROUTINE flux_ausmdv( rho_l, rho_r, v1_l , v1_r, &
                         v2_l , v2_r , p_l  , p_r , &
                         flux_side)
!===================================================================================================================================
! AUSMDV flux
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
REAL             :: rho_u,rho_u_sq,c_l,c_r,cm,s
REAL             :: alpha_l,alpha_r
LOGICAL          :: tmpa=.FALSE.,tmpb=.FALSE.
!===================================================================================================================================

! Calculate left/right energy and enthalpy
e_r=gamma1q*p_r+0.5*rho_r*(v1_r*v1_r+v2_r*v2_r)
e_l=gamma1q*p_l+0.5*rho_l*(v1_l*v1_l+v2_l*v2_l)

H_r=(e_r + p_r)/rho_r
H_l=(e_l + p_l)/rho_l

! Maximum speed of sound
c_l = sqrt(gamma*p_l/rho_l)
c_r = sqrt(gamma*p_r/rho_r)
cm = max(c_r,c_l)

alpha_r = 2.*p_r/rho_r/(p_l/rho_l+p_r/rho_r)
alpha_l = 2.*p_l/rho_l/(p_l/rho_l+p_r/rho_r)

IF(abs(v1_l).le.cm)THEN
    p_plus = 0.25*p_l*(v1_l+cm)*(v1_l+cm)/(cm*cm)*(2.-v1_l/cm)
    IF (v1_l.GE.0.) then
        u_plus = v1_l + alpha_l*(v1_l-cm)*(v1_l-cm)
    ELSE
        u_plus =        alpha_l*(v1_l+cm)*(v1_l+cm)
    END IF
ELSE
    IF(v1_l.GE.0) THEN
        u_plus = v1_l
        p_plus = p_l
    ELSE
        u_plus = 0.
        p_plus = 0.
    END IF
ENDIF

IF(abs(v1_r).le.cm)THEN
    p_minus = 0.25*p_r*(v1_r-cm)*(v1_r-cm)/(cm*cm)*(2.+v1_r/cm)
    IF (v1_r.GE.0.) then
        u_minus =      - alpha_r*(v1_l-cm)*(v1_l-cm)
    ELSE
        u_minus = v1_r - alpha_r*(v1_r+cm)*(v1_r+cm)
    END IF
ELSE
    IF(v1_r.GE.0) THEN
        u_minus = 0.0
        p_minus = 0.0
    ELSE
        u_minus = v1_r
        p_minus = p_r
    END IF
ENDIF

! Mass flux
rho_u = u_plus*rho_l + u_minus*rho_r

! momentum flux (blending between AUSMD and AUSMV)
s = MIN(1.,10.*ABS(p_r-p_l)/MIN(p_r,p_l))
rho_u_sq =             0.5*(1.+s)*(rho_l*v1_l*u_plus + rho_r*v1_r*u_minus)   ! AUSMV flux
rho_u_sq = rho_u_sq + 0.25*(1.-s)*(rho_u*(v1_r+v1_l)-ABS(rho_u)*(v1_r-v1_l)) ! AUSMD flux

flux_side(1) = rho_u
flux_side(2) = rho_u_sq + (p_plus+p_minus)
flux_side(3) = 0.5*(rho_u*(v2_r+v2_l)-ABS(rho_u)*(v2_r-v2_l))
flux_side(4) = 0.5*(rho_u*(H_r+H_l) - ABS(rho_u)*(H_r-H_l))

! entropy fix
tmpa = (v1_l-c_l < 0.) .AND. (v1_r-c_r > 0.)
tmpb = (v1_l+c_l < 0.) .AND. (v1_r+c_r > 0.)
IF(tmpa .AND. .NOT.(tmpb))THEN
  flux_side(:) = flux_side(:) - 0.125*((v1_r-c_r)-(v1_l-c_l))*(rho_r*(/1.,v1_r,v2_r,h_r/)-rho_l*(/1.,v1_l,v2_l,h_l/))
END IF
IF(.NOT. tmpa .AND. tmpb)THEN
  flux_side(:) = flux_side(:) - 0.125*((v1_r+c_r)-(v1_l+c_l))*(rho_r*(/1.,v1_r,v2_r,h_r/)-rho_l*(/1.,v1_l,v2_l,h_l/))
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE flux_ausmdv

END MODULE MOD_flux_ausmdv
