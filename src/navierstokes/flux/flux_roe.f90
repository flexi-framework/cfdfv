MODULE MOD_flux_roe
!===================================================================================================================================
! Roe flux
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE flux_roe
   MODULE PROCEDURE flux_roe
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: flux_roe
!===================================================================================================================================

CONTAINS

SUBROUTINE flux_roe( rho_l, rho_r, &
                     v1_l, v1_r,   &
                     v2_l, v2_r,   &
                     p_l, p_r,     &
                     flux_side     )
!===================================================================================================================================
! Calculation of Roe flux
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars, ONLY: gamma,gamma1q
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)             :: rho_l, rho_r
REAL,INTENT(IN)             :: v1_l , v1_r
REAL,INTENT(IN)             :: v2_l , v2_r
REAL,INTENT(IN)             :: p_l  , p_r
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)            :: flux_side(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                         :: m1_l , m1_r
REAL                         :: m2_l , m2_r
REAL                         :: e_l  , e_r
REAL                         :: H_r, H_l
REAL                         :: sqrt_rho_r, sqrt_rho_l, sum_sqrt_rho_q
REAL                         :: v1_quer,v2_quer,H_quer,c_quer
REAL                         :: a1,a2,a3,a4
REAL                         :: r1(4),r2(4),r3(4),r4(4)
REAL                         :: delta_rho,delta_m1,delta_m2,delta_e,delta_eq
REAL                         :: c_quer_q
REAL                         :: gamma1,gamma2,gamma3,gamma4
REAL                         :: f1_r(4),f1_l(4)
!===================================================================================================================================

! Calculate left/right conservative variables
m1_l = rho_l * v1_l
m1_r = rho_r * v1_r
m2_l = rho_l * v2_l
m2_r = rho_r * v2_r
e_r  = gamma1q * p_r + 0.5 * rho_r * &
       (v1_r * v1_r + v2_r * v2_r)
e_l  = gamma1q * p_l + 0.5 * rho_l * &
       (v1_l * v1_l + v2_l * v2_l)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate left/right enthalpy
H_r=(e_r + p_r)/rho_r
H_l=(e_l + p_l)/rho_l
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate sqrt(rho)
sqrt_rho_r = SQRT(rho_r)
sqrt_rho_l = SQRT(rho_l)
sum_sqrt_rho_q = 1. / (sqrt_rho_l + sqrt_rho_r)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate Roe mean values
v1_quer=(sqrt_rho_r*v1_r+sqrt_rho_l*v1_l) * sum_sqrt_rho_q
v2_quer=(sqrt_rho_r*v2_r+sqrt_rho_l*v2_l) * sum_sqrt_rho_q
H_quer=(sqrt_rho_r* H_r+sqrt_rho_l* H_l) * sum_sqrt_rho_q
c_quer=SQRT((gamma-1.)*(H_quer-0.5*(v1_quer*v1_quer + v2_quer*v2_quer)))
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate mean eigenvalues
a1=v1_quer-c_quer
a2=v1_quer
a3=v1_quer
a4=v1_quer+c_quer
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate mean eigenvectors
r1(1)=1.
r1(2)=a1
r1(3)=v2_quer
r1(4)=H_quer-v1_quer*c_quer

r2(1)=1.
r2(2)=v1_quer
r2(3)=v2_quer
r2(4)=0.5*(v1_quer*v1_quer + v2_quer*v2_quer)

r3(1)=0.
r3(2)=0.
r3(3)=1.
r3(4)=v2_quer

r4(1)=1.
r4(2)=a4
r4(3)=v2_quer
r4(4)=H_quer+v1_quer*c_quer
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate differences
delta_rho = rho_r - rho_l
delta_m1  = m1_r  - m1_l
delta_m2  = m2_r  - m2_l
delta_e   = e_r   - e_l
delta_eq  = delta_e-(delta_m2-v2_quer*delta_rho)*v2_quer
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate wave strengths gamma1,...,gamma4
c_quer_q = 1./c_quer
gamma2   = -(gamma-1.)*c_quer_q*c_quer_q*          &
           (delta_rho*(v1_quer*v1_quer-H_quer) &
           + delta_eq &
           - delta_m1*v1_quer )
gamma1   = -0.5*c_quer_q*(delta_m1-delta_rho*(v1_quer+c_quer))-0.5*gamma2
gamma4   = delta_rho - gamma1 - gamma2 
gamma3   = delta_m2 - v2_quer*delta_rho
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculalte the physical fluxes
f1_r(1) = m1_r
f1_r(2) = m1_r*v1_r + p_r
f1_r(3) = m1_r*v2_r
f1_r(4) = v1_r*(e_r + p_r)

f1_l(1) = m1_l
f1_l(2) = m1_l*v1_l + p_l
f1_l(3) = m1_l*v2_l
f1_l(4) = v1_l*(e_l + p_l)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate Roe flux
flux_side(:) = 0.5*( f1_r(:)+ f1_l(:) &
               - gamma1*ABS(a1)*r1(:) &
               - gamma2*ABS(a2)*r2(:) &
               - gamma3*ABS(a3)*r3(:) &
               - gamma4*ABS(a4)*r4(:) )
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE flux_roe

END MODULE MOD_flux_roe
