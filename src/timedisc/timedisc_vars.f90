MODULE MOD_TimeDisc_Vars
!===================================================================================================================================
! Module with the timedisc variables
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                 :: CFL                     ! CFL number
REAL                 :: DFL                     ! DFL number
REAL                 :: t
! Profiling
REAL                 :: time_overall            ! overall time
REAL                 :: IO                      ! Time for I/O
! Constants controlling the numerical scheme
INTEGER              :: TimeOrder               ! Order of time integration
LOGICAL              :: TimeStep1D              ! 1D Problem? (for time step calculation)
! Constants regarding the abort criterion 
LOGICAL              :: stationary              ! stationary solution expected?
INTEGER              :: MaxIter                 ! No. of max. iterations
REAL                 :: Stoptime                ! Simulation endtime
INTEGER              :: IniIterationNumber      ! iteration number to start with
REAL                 :: StartTime
REAL                 :: AbortResidual
INTEGER              :: AbortVariable
CHARACTER(LEN=3)     :: AbortVarName
REAL                 :: Cl_AbortResidual
REAL                 :: Cd_AbortResidual
LOGICAL              :: AbortResidualCl,AbortResidualCd
LOGICAL              :: Restart
REAL                 :: RestartTime
! io
INTEGER              :: printiter
REAL                 :: printtime
! Constants for RK time stepping
INTEGER              :: RK                      ! switch for RK time stepping
INTEGER              :: nRKstages               ! number of RK stages used
REAL                 :: RKcoeff(0:5)            ! RK coefficient for stage
LOGICAL              :: implicit                ! logical for explict/implicit switch
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE MOD_TimeDisc_Vars
