MODULE MOD_Initialcondition_Vars
!===================================================================================================================================
! Global vars
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! GLOBAL VARIABLES
!===================================================================================================================================
! Initial Conditions
INTEGER              :: icType                   ! Type of initial condition    
! Variables for Multi-Domain Initial Condition
INTEGER              :: ndomains                 ! number of domains            
INTEGER,ALLOCATABLE  :: domain_ID(:)             ! goes from 1 to ndomains      
REAL                 :: RP_1D_interface          ! Interface                    
REAL                 :: alpha
REAL,ALLOCATABLE     :: RefState(:, :)           ! Primitive variable state in domain
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE MOD_Initialcondition_Vars
