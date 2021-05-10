MODULE MOD_InitialCondition
!===================================================================================================================================
! Initialization of stream
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitInitialCondition
   MODULE PROCEDURE InitInitialCondition
END INTERFACE
INTERFACE InitialCondition
   MODULE PROCEDURE InitialCondition
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
  PUBLIC  :: InitInitialCondition,InitialCondition
!===================================================================================================================================

CONTAINS

SUBROUTINE InitInitialCondition()
!===================================================================================================================================
! Get initial flow conditions from init file
!===================================================================================================================================
! MODULES
USE MOD_Readintools
USE MOD_Initialcondition_Vars
!USE MOD_TimeDisc_Vars                ,ONLY:Restart
USE MOD_Globals                      ,ONLY:NVAR,RHO,P,V1,V2
USE MOD_Equation_Vars                ,ONLY:intExactFunc,gamma
USE MOD_Initialcondition_Vars        ,ONLY:icType,Refstate
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL           :: MachNumber, AngleOfAttack,c,v
INTEGER        :: i
!===================================================================================================================================

WRITE(*,*) "-[Initial Conditions]---------------------------------------"
icType = GETINT('ICtype')
SELECT CASE(icType)
CASE(1)
  nDomains = GETINT('nDomains','1')
  !WRITE(*,'(a, I2, a)') '  Homogeneous Initial Conditions with ', ndomains, ' domain(s)'
  ALLOCATE ( domain_ID(1:ndomains))
  ALLOCATE (Refstate(1:nDomains,NVAR))
  DO i = 1,INT(ndomains)
    WRITE(*,'(a, I2, a)') '    Domain ', i, ':'
    ! Get the Domain_ID of first Domain
    domain_ID(i)=GETINT('DomainID')
    WRITE(*,'(a, I2)') '      Domain ID: ', i
    ! Get density of first Domain
    Refstate(i, RHO) = GETREAL('Rho')
    ! Get Mach-Number
    MachNumber = GETREAL('Mach')
    WRITE(*,'(a, F10.5)') '      Mach Number: ', MachNumber
    ! Get Angle of Attack
    AngleOfAttack = GETREAL('Alpha')
    alpha= AngleOfAttack
    WRITE(*,'(a, F10.5)') '      Angle of Attack [deg]: ', AngleOfAttack
    ! Get pressure first Domain
    Refstate(i, P) = GETREAL('Pressure')
    WRITE(*,'(a, F10.5)') '      Pressure: ', Refstate(i, P)
    c = SQRT(gamma * Refstate(i, P) / Refstate(i, RHO))
    v = MachNumber * c
    Refstate(i, V1) = v * COS(AngleOfAttack * ACOS(-1.) / 180.)
    Refstate(i, V2) = v * SIN(AngleOfAttack * ACOS(-1.) / 180.)
  END DO
CASE(2)
  intExactFunc = GETINT('ExactFunc')
  WRITE(*,'(a, I2)') '  Initial Condition with exact function # ', intExactFunc
  IF(intExactFunc == 5) THEN
    ALLOCATE(Refstate(2, 1:NVAR))
    Refstate = 0.
  ! 1D Riemann Problem: Get Location of Interface
    RP_1D_interface = GETREAL('RP_1D_interface','0.5')
    WRITE(*,'(a, F10.5)') '  Position of interface for 1D Riemann Problem: ', RP_1D_interface
    Refstate(1,1:NVAR) = GETREALARRAY('StateLeft',4)
    Refstate(2,1:NVAR) = GETREALARRAY('StateRight',4)
  END IF
CASE DEFAULT
  WRITE(*,*) 'ERROR InitialCondition not known!'
  WRITE(*,*) icType
  STOP
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE InitInitialCondition

SUBROUTINE InitialCondition()
!===================================================================================================================================
! Set init. flow to all cells
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ExactFunc
USE MOD_EoS
USE MOD_Readin               ,ONLY:CGNS_ReadSolution
USE MOD_Initialcondition_Vars,ONLY:icType,Refstate,nDomains
USE MOD_Equation_Vars        ,ONLY:intExactFunc
USE MOD_TimeDisc_Vars        ,ONLY:Restart
USE MOD_Mesh_Vars            ,ONLY:FirstElem
USE MOD_Mesh_Vars            ,ONLY:tElem
USE MOD_Readintools          ,ONLY:GETINT
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                    :: x(2)
TYPE(tElem), POINTER    :: aElem
!===================================================================================================================================

IF (Restart) THEN
! Read IC from previous flow solution
  CALL CGNS_ReadSolution()
ELSE
  ! Set initial IC
  SELECT CASE(icType)
  CASE(0)
    ! Cell Test
    aElem => FirstElem
    DO WHILE(ASSOCIATED(aElem))
      aElem%pvar(V1)  = 0.
      aElem%pvar(V2)  = 0.
      aElem%pvar(RHO) = 1. * aElem%ID
      aElem%pvar(P)   = 1.
   ! Next Element
      aElem => aElem%NextElem
    END DO
  CASE(1)
    ! 1: Multi-Domain homogeneous initial condition
    ! 2: homogeneous initial condition over all domains
    aElem => FirstElem
    DO WHILE(ASSOCIATED(aElem))
      IF (nDomains == 1) THEN
        aElem%pvar(RHO:P) = Refstate(1, RHO:P)
      ELSE
        aElem%pvar(RHO:P) = Refstate(aElem%domain, RHO:P)
      END IF
   ! Next Element
      aElem=>aElem%nextElem
    END DO
  CASE(2)
    ! Exact Function
    aElem => firstElem
    DO WHILE(ASSOCIATED(aElem))
      x = aElem%Bary
      CALL ExactFunc(intExactFunc, x, 0., aElem%pVar(1:4))
      aElem => aElem%nextElem
    END DO
  CASE DEFAULT
    WRITE(*,*) 'ERROR InitialCondition not known!'
    WRITE(*,*) icType
    STOP
  END SELECT
END IF
!-----------------------------------------------------------------------------------------------------------------------------------

! Compute cvar
aElem => FirstElem
DO WHILE(ASSOCIATED(aElem))
  CALL PrimCons(aElem%pvar, aElem%cvar)
  aElem => aElem%nextElem
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE InitialCondition

END MODULE MOD_InitialCondition
