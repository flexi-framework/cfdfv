MODULE MOD_flux_hllc
!-----------------------------------------------------------------------------------------------------------------------------------
! HLLC flux 
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE flux_hllc
   MODULE PROCEDURE flux_hllc
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: flux_hllc
!===================================================================================================================================

CONTAINS

SUBROUTINE flux_hllc( rho_l, rho_r, &
                      v1_l, v1_r,   &
                      v2_l, v2_r,   &
                      p_l, p_r,     &
                      flux_side     )
!===================================================================================================================================
! Computation of the HLLC flux
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Equation_Vars, ONLY: gamma,gamma1,gamma1q
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
REAL                        :: e_l, e_r
REAL                        :: fl(4),fr(4),ul(4),ur(4),us(4)
REAL                        :: c_l,c_r,H_l,H_r
REAL                        :: um, vm, Hm, cm
REAL                        :: arp, alm, as
REAL                        :: sqrt_rho_r, sqrt_rho_l, sum_sqrt_rho_q
REAL                        :: rho_lq, rho_rq, factor
!===================================================================================================================================

! Calculation of auxiliary values
rho_lq         = 1. / rho_l
rho_rq         = 1. / rho_r
sqrt_rho_r     = SQRT(rho_r)
sqrt_rho_l     = SQRT(rho_l)
sum_sqrt_rho_q = 1. / (sqrt_rho_r + sqrt_rho_l)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculate Energies
e_r = gamma1q * p_r + 0.5 * rho_r * &
      (v1_r * v1_r + v2_r * v2_r)
e_l = gamma1q * p_l + 0.5 * rho_l * &
      (v1_l * v1_l + v2_l * v2_l)
!-----------------------------------------------------------------------------------------------------------------------------------
! Left conservative state vector
ul(RHO) = rho_l
ul(M1) = rho_l * v1_l
ul(M2) = rho_l * v2_l
ul(E) = e_l
! Right conservative state vector
ur(RHO) = rho_r
ur(M1) = rho_r * v1_r
ur(M2) = rho_r * v2_r
ur(E) = e_r
!-----------------------------------------------------------------------------------------------------------------------------------
! Flux in left cell
fl(RHO) = ul(2)
fl(M1) = fl(1)*v1_l+p_l
fl(M2) = fl(1)*v2_l
fl(E) = v1_l*(e_l+p_l)
! Flux in right cell
fr(RHO) = ur(2)
fr(M1) = fr(1)*v1_r+p_r
fr(M2) = fr(1)*v2_r
fr(E) = v1_r*(e_r+p_r)
!-----------------------------------------------------------------------------------------------------------------------------------
! Calculation of sound speeds
c_l=SQRT(gamma * p_l * rho_lq)
c_r=SQRT(gamma * p_r * rho_rq)
! Left / right enthalpy
H_l = (e_l + p_l) * rho_lq
H_r = (e_r + p_r) * rho_rq
!-----------------------------------------------------------------------------------------------------------------------------------
! Roe mean values
um=(sqrt_rho_r*v1_r+sqrt_rho_l*v1_l) * sum_sqrt_rho_q
vm=(sqrt_rho_r*v2_r+sqrt_rho_l*v2_l) * sum_sqrt_rho_q
Hm=(sqrt_rho_r*H_r +sqrt_rho_l*H_l ) * sum_sqrt_rho_q
cm=SQRT(gamma1*(Hm-0.5*(um*um+vm*vm)))
!-----------------------------------------------------------------------------------------------------------------------------------
! Signal speeds
arp=MAX(v1_r+c_r,um+cm)
alm=MIN(v1_l-c_l,um-cm)
!-----------------------------------------------------------------------------------------------------------------------------------
! HLLC-Flux
IF (alm > 0.) THEN
  flux_side = fl
ELSEIF (arp < 0.) THEN
  flux_side = fr
ELSE
  as = (p_r-p_l+ul(M1)*(alm-v1_l)-ur(M1)*(arp-v1_r))/(rho_l*(alm-v1_l)-rho_r*(arp-v1_r))
  IF ((alm <= 0.).AND.(as >= 0.)) THEN
    factor = rho_l * (alm-v1_l)/(alm-as)
    us(RHO) = factor
    us(M1)  = as * factor
    us(M2)  = v2_l * factor
    us(E)   = factor * (e_l/rho_l+(as-v1_l)*(as+p_l/(rho_l*(alm-v1_l))))

    flux_side = fl + alm * (us - ul)
  ELSE
    factor = rho_r * (arp-v1_r)/(arp-as)
    us(RHO) = factor
    us(M1)  = as * factor
    us(M2)  = v2_r * factor
    us(E)   = factor * (e_r/rho_r+(as-v1_r)*(as+p_r/(rho_r*(arp-v1_r))))

    flux_side = fr + arp * (us - ur)
  END IF
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE flux_hllc

END MODULE MOD_flux_hllc
