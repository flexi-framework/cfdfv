MODULE MOD_LinearSolver_Vars
!===================================================================================================================================
! Contains global variables used by the DG modules.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! number of DOF
INTEGER                               :: nKDim                      ! max Krylov space
INTEGER                               :: nNewtonIter                ! max Newton iterations
INTEGER                               :: nNewtonIterGlobal          ! max Newton iterations global for the calculation
INTEGER                               :: nGMRESIterGlobal           ! max GMRES iterations global for the calculation
INTEGER                               :: nInnerNewton               ! max Newton iterations for actual Newton stage
INTEGER                               :: nInnerGMRES                ! max GMRES iterations for actual GMRES stage 
INTEGER                               :: iterGlobal
LOGICAL                               :: Precond                    ! use LU-SGS as a preconditioner
! epsilons
REAL                                  :: rEps0,srEps0               ! SQRT(EPSILON(0.0)), used for EpsFD
REAL                                  :: Eps2Newton                 ! square of newton relative epsilon
REAL                                  :: EpsGMRES
REAL                                  :: gammaEW                    ! gamma parameter for Eisenstat Walker
REAL,ALLOCATABLE,DIMENSION(:,:)       :: XK,R_XK
REAL,ALLOCATABLE,DIMENSION(:,:,:)     :: Dinv
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)   :: LowerUpper
INTEGER,ALLOCATABLE,DIMENSION(:,:,:)  :: ElemToElem
! 
LOGICAL                               :: ImplicitInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MOD_LinearSolver_Vars

