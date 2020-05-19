MODULE MOD_Globals
!===================================================================================================================================
! Module containing the global variables
! These variables are used in almost all subroutines
!
! parameters are helpful for replacing index numbers by their particular meaning. used in pvar, cvar, and xy direction           
! they are globally defined, used instead of an index! CAPITAL LETTERS!           
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!===================================================================================================================================
  ! parameters for cvar
  INTEGER, PARAMETER      :: RHO  = 1               ! density (pvar and cvar)      
  INTEGER, PARAMETER      :: M1   = 2               ! x momentum                   
  INTEGER, PARAMETER      :: M2   = 3               ! y momentum                   
  INTEGER, PARAMETER      :: E    = 4               ! Energy                       
  ! parameters for pvar
  INTEGER, PARAMETER      :: V1   = 2               ! x velocity                   
  INTEGER, PARAMETER      :: V2   = 3               ! y momentum                   
  INTEGER, PARAMETER      :: P    = 4               ! pressure                     
  ! x,y directions
  INTEGER, PARAMETER      :: X_DIR = 1              ! x comp of a vector           
  INTEGER, PARAMETER      :: Y_DIR = 2              ! y comp of a vector           
  !parameter for BCs
  INTEGER, PARAMETER      :: SLIPWALL = 1           ! index of  slip wall BC edges
  INTEGER, PARAMETER      :: WALL     = 2           ! index of wall BC edges
  INTEGER, PARAMETER      :: INFLOW   = 3           ! index of supersonic inflow BC edges
  INTEGER, PARAMETER      :: OUTFLOW  = 4           ! index of supersonic outflow BC edges
  INTEGER, PARAMETER      :: CHARACTERISTIC = 5     ! index of characteristic BC edges (subsonic inflow)
  INTEGER, PARAMETER      :: EXACTSOL = 6           ! index of exact solution BC edges
  INTEGER, PARAMETER      :: PERIODIC = 7           ! index of periodic BC edges
  INTEGER, PARAMETER      :: PRESSURE_OUT = 8       ! index of subsonic pressure outflow BC
  ! parameter for BC of a Cartesian mesh
  INTEGER, PARAMETER      :: BOTTOM = 1             ! index of bottom boundary      
  INTEGER, PARAMETER      :: RIGHT  = 2             ! index of right boundary       
  INTEGER, PARAMETER      :: TOP    = 3             ! index of top boundary         
  INTEGER, PARAMETER      :: LEFT   = 4             ! index of left boundary        
  ! parameter for mesh type
  INTEGER, PARAMETER      :: UNST    = 0            ! index for unstructured mesh  
  INTEGER, PARAMETER      :: CART    = 1            ! index for cartesian mesh     
  ! parameter for variable conversion
  INTEGER, PARAMETER      :: PVAR_TO_CVAR = 1       ! primitive to conservative    
  INTEGER, PARAMETER      :: CVAR_TO_PVAR = 2       ! vice versa                   
  ! parameter for flux functions
  INTEGER, PARAMETER      :: GOD   = 1              ! Godunov        flux function 
  INTEGER, PARAMETER      :: ROE   = 2              ! Roe            flux function 
  INTEGER, PARAMETER      :: HLL   = 3              ! HLLE           flux function 
  INTEGER, PARAMETER      :: HLLE  = 4              ! HLLE           flux function 
  INTEGER, PARAMETER      :: HLLC  = 5              ! HLLC           flux function 
  INTEGER, PARAMETER      :: LXF   = 6              ! Lax-Friedrichs flux function
  INTEGER, PARAMETER      :: STW   = 7              ! Steger-Warming flux function
  INTEGER, PARAMETER      :: CEN   = 8              ! central        flux function
  INTEGER, PARAMETER      :: AUSMD = 9              ! AUSMD          flux function
  INTEGER, PARAMETER      :: AUSMDV=10              ! AUSMDV         flux function
  INTEGER, PARAMETER      :: VANLEER=11             ! VAN LEER       flux function
  ! parameter for timestepping information
  INTEGER, PARAMETER      :: GLOBAL_STEP = 0        ! global timesteps      
  INTEGER, PARAMETER      :: LOCAL_STEP  = 1        ! local timesteps       
  ! parameter for limiter information
  INTEGER, PARAMETER      :: BARTHJESPERSEN  = 1    
  INTEGER, PARAMETER      :: VENKATAKRISHNAN = 2 
  ! parameters for different equation types
  INTEGER, PARAMETER      :: EULER_EQ  = 1            
  INTEGER, PARAMETER      :: NAVIER_STOKES_EQ = 2  
  ! parameters for boundary conditions
  INTEGER, PARAMETER      :: MEAN_VALUES = 1 
  INTEGER, PARAMETER      :: EDGE_VALUES = 2 
  ! global parameters
  INTEGER, PARAMETER      :: NVAR = 4
  INTEGER, PARAMETER      :: NDIM = 2
  ! I/O Formats
  INTEGER, PARAMETER      :: CGNS  = 1
  INTEGER, PARAMETER      :: CURVE = 2
  INTEGER, PARAMETER      :: DAT   = 3
  ! ENTROPY FIX
  REAL, PARAMETER         :: LAMBDA = 0.1
  ! Real tolerance
  REAL, PARAMETER         :: REALTOL = 1.E-7
  ! default unit for std_out
  INTEGER,PARAMETER       ::UNIT_stdOut=6
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE MOD_Globals
