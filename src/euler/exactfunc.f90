MODULE MOD_ExactFunc
!===================================================================================================================================
! Exact function
! Used for initialization, reference state and BC
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: ExactFunc
!===================================================================================================================================

CONTAINS

SUBROUTINE ExactFunc(iExactFunc, x, time, pvar)
!===================================================================================================================================
! Exact function

!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars            ,ONLY:xMin,xMax,yMin,yMax,dxref
USE MOD_Equation_Vars        ,ONLY:gamma,gamma1
USE MOD_Equation_Vars        ,ONLY:Pi,sqrt3_q
USE MOD_InitialCondition_Vars,ONLY:RefState,RP_1D_interface
USE MOD_EOS                  ,ONLY:ConsPrim
USE MOD_ExactRiemann
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)          :: iExactFunc
REAL,INTENT(IN)             :: time           ! time
REAL,INTENT(IN)             :: x(2)           ! location
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)            :: pvar(NVAR)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
! RMI
REAL             :: xLen, yLen
! Gaussian Pulse:
REAL             :: amplitude, hwidth
REAL             :: xPeak, yPeak, explog, hwidth2, dr2
! Double Mach Reflection
REAL             :: cvar(NVAR), cSum(NVAR), cDiff(NVAR), factor
! 1D Riemann Problem
REAL             :: rho_l, rho_r, u_l, u_r, p_l, p_r, c_l, c_r
REAL             :: s
! viscous convtest
REAL             :: resu(NVAR),Omega,Frequency,a
! sine wave v+c 
REAL             :: cons0(NVAR),prims0(NVAR),H0,c0,R(4,4)
!===================================================================================================================================

SELECT CASE (iExactFunc)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(1)
!-----------------------------------------------------------------------------------------------------------------------------------
! Richtmyer-Meshkov Instability (independent of time, only initial condition)
  pvar(V1)  = 0.
  pvar(V2)  = 0.
  pvar(RHO) = 1.
  pvar(P)   = 1.
  xLen      = xMax - xMin
  yLen      = yMax - yMin
  IF (x(1) >= 0.3 * xLen + (1./30.) * xLen * COS(2.*Pi*3./yLen*x(2))) THEN
    pvar(RHO) = 0.25
  END IF
  IF ((x(1) <= 0.1 * xLen).AND.       &
      (x(1) >= (1./30.) * xLen)) THEN
    pvar(P)   = 4.9
    pvar(RHO) = 4.22
  END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! Gaussian Pressure Pulse (independent of time, only initial condition)
  pvar(RHO) = 1.
  pvar(V1)  = 0.
  pvar(V2)  = 0.

  hwidth    = MIN((xMax - xMin), (yMax - yMin)) / 100. * 6.
  amplitude = 1.

  xPeak     = xmin + 0.5 * (xmax - xmin)
  yPeak     = ymin + 0.5 * (ymax - ymin)
  explog    = EXP(-LOG(2.0))
  hwidth2   = hwidth * hwidth
  dr2       =  (x(1) - xPeak)*(x(1) - xPeak) + &
               (x(2) - yPeak)*(x(2) - yPeak)

  pvar(P)   = 1. + amplitude *(explog) ** (dr2 / hwidth2)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! Sinewave: Test case for convergence tests
   Frequency=1.
   Amplitude=0.1
   Omega=Pi*Frequency
   a = 2.*Pi
   resu(1:3) = 2.+ Amplitude*SIN(Omega*SUM(x(:)) - a*time)
   resu(4)   = Resu(1)*Resu(1)
   CALL ConsPrim(pvar,resu)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! Double Mach Reflection: Undisturbed Shock for initial condition and upper boundary
  factor = Pi / (0.2*dxRef)

  cSum(1) = 9.4
  cSum(2) = 57.157676649772950686
  cSum(3) = -33.
  cSum(4) = 566.

  cDiff(1) = 6.6
  cDiff(2) = 57.157676649772950686
  cDiff(3) = - 33.
  cDiff(4) = 561.

  cvar(:) = 0.5*(cSum-cDiff*TANH((x(1)-(1./6. + (20.*time + x(2))*sqrt3_q))*factor));

  pvar(1) = cvar(1)
  pvar(2) = cvar(2) / cvar(1)
  pvar(3) = cvar(3) / cvar(1)
  pvar(4) = gamma1*(cvar(4) - 0.5*(cvar(2)*cvar(2) + cvar(3)*cvar(3))/cvar(1))
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(5)
!-----------------------------------------------------------------------------------------------------------------------------------
! 1D Riemann Problem in x direction
  rho_l = RefState(1,RHO)
  u_l   = RefState(1,V1)
  p_l   = RefState(1,P)
  c_l   = SQRT(gamma*p_l/rho_l)

  rho_r = RefState(2,RHO)
  u_r   = RefState(2,V1)
  p_r   = RefState(2,P)
  c_r   = SQRT(gamma*p_r/rho_r)

  IF (time == 0.) THEN
    IF (x(1) <= RP_1D_interface) THEN
      pvar(RHO) = rho_l
      pvar(V1)  = u_l
      pvar(P)   = p_l
    ELSE
      pvar(RHO) = rho_r
      pvar(V1)  = u_r
      pvar(P)   = p_r
    END IF
  ELSE
    s     = (x(1) - RP_1D_interface) / time
    CALL exact_riemann(gamma, rho_l, rho_r, pvar(RHO), &
                              u_l,   u_r,   pvar(V1),  &
                              p_l,   p_r,   pvar(P),   &
                              c_l,   c_r,   s          )
  END IF

  pvar(V2) = 0.
!-----------------------------------------------------------------------------------------------------------------------------------
CASE(6)
!-----------------------------------------------------------------------------------------------------------------------------------
! 1D sine wave
   Frequency=1.
   Amplitude=0.00001
   Omega=Pi*Frequency
   cons0=(/1., .1, 0., 1./) !define base state
   CALL ConsPrim(prims0,cons0)
   c0=SQRT(gamma*prims0(P)/prims0(RHO))
   H0=(cons0(E)+prims0(P))/prims0(RHO)
   !rotation matrix
   R(:,1)=(/1.              , 1.                                              , 0.        , 1.               /)
   R(:,2)=(/prims0(V1)-c0   , prims0(V1)                                      , 0.        , prims0(V1)+c0    /)
   R(:,3)=(/prims0(V2)      , prims0(V2)                                      , 1.        , prims0(V2)       /)
   R(:,4)=(/H0-prims0(V1)*c0, (prims0(V1)*prims0(V1)+prims0(V2)*prims0(V2))/2., prims0(V2), H0+prims0(V1)*c0 /)
   !resu = base state + sin(x) in invariant variables (diagonalized system). This sin(x) has to be transformed back to cons.
   resu=cons0+MATMUL(TRANSPOSE(R),(/0., 0., 0., Amplitude*SIN(Omega*(x(1)-(c0+prims0(V1))*time)) /) )
   CALL ConsPrim(pvar,resu)
!-----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
    WRITE(*,*) ' Error in ExactFunc.f90. Exact Function unknown.'
    STOP
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ExactFunc

END MODULE MOD_ExactFunc
