MODULE MOD_Equation_vars
!===================================================================================================================================
! Variables for the equation system
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!===================================================================================================================================

! General constants
REAL                 :: Pi                      ! Pi
! Constants regarding the equations to solve
LOGICAL              :: CalcSource              ! calc source
REAL                 :: R                       ! gas constant
REAL                 :: gamma                   ! Ratio of spec. heats
REAL                 :: gamma1                  ! gamma - 1
REAL                 :: gamma2                  ! gamma - 2
REAL                 :: gamma1q                 ! 1. / (gamma - 1.)
REAL                 :: CFL                     ! CFL number
!REAL                 :: cp                      ! specific heat capacity
REAL                 :: Pr                      ! Prandtl number
REAL                 :: mu                      ! viscosity
! Constants controlling the numerical scheme.
INTEGER              :: iFlux                   ! Numerical scheme
! Constants regarding the abort criterion 
! other Constants 
INTEGER              :: intExactFunc            ! Exact Function
INTEGER              :: SourceFunc              ! Source Term
REAL                 :: sqrt3_q, sqrt3, sqrt2

END MODULE MOD_Equation_vars 
