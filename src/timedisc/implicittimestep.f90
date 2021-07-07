MODULE MOD_ImplicitTimeStep
!===================================================================================================================================
! Module containing the implicit time stepping scheme
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ImplicitTimeStep
   MODULE PROCEDURE ImplicitTimeStep
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: ImplicitTimeStep
!===================================================================================================================================


CONTAINS

SUBROUTINE ImplicitTimeStep(t,dt,iter,res_iter)
!===================================================================================================================================
! Euler Implicit Time Integration
! Non-linear equations requires to use a Newton Method with internal subiteration
! using a GMRES method
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,          ONLY: tElem,tSide
USE MOD_Mesh_Vars,          ONLY: nElems,Elems
USE MOD_FV,                 ONLY: FV_TimeDerivative
USE MOD_EoS,                ONLY: ConsPrim
USE MOD_LinearSolver_Vars,  ONLY: nNewtonIter,XK,R_XK,iterGlobal,Precond
USE MOD_LinearSolver_Vars,  ONLY: gammaEW,nNewtonIterGlobal,Eps2Newton
USE MOD_LinearSolver,       ONLY: VectordotProduct, GMRES_M, BuildMatrix
USE MOD_Analyze,            ONLY: GlobalResidual
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)      :: iter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)      :: t,dt
REAL,INTENT(OUT)        :: res_iter(NVAR+2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem), POINTER    :: aElem
REAL                    :: alpha, beta,time
! mapping from Pointer - to - Array list
REAL             :: Norm2_F_X0,Norm2_F_Xk,Norm2_F_Xk_old
REAL             :: F_X0  (1:NVAR,1:nElems)
REAL             :: Q     (1:NVAR,1:nElems)
REAL             :: F_Xk  (1:NVAR,1:nElems)
!REAL             :: Xk    (1:NVAR,1:nElems)
REAL             :: DeltaX(1:NVAR,1:nElems)
REAL             :: AbortCritGMRES,gammaA,gammaB
INTEGER          :: nInnerNewton
INTEGER          :: iElem,iVar
!===================================================================================================================================

! input parameters for Netwon
iterGlobal = iter
alpha = 1.
beta  = 1.
! not required for Euler time integration
!Q     = U
!$omp parallel do private(aElem)
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  Q(:,iElem) = aElem%cvar(:)
  CALL ConsPrim(aElem%pvar, aElem%cvar)
END DO
!$omp end parallel do

!-----------------------------------------------------------------------------------------------------------------------------------
! Newton
time  = t + beta*dt
!-----------------------------------------------------------------------------------------------------------------------------------
CALL FV_TimeDerivative(time, iter )
!$omp parallel do private(aElem)
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  !CALL ConsPrim(aElem%pvar, aElem%cvar)
  F_X0(:,iElem) = aElem%cvar(:) - Q(:,iElem) -Alpha*dt*aElem%u_t(:)
  Xk  (:,iElem) = aElem%cvar(:)
  R_Xk(:,iElem) = aElem%u_t(:)
  F_Xk(:,iElem) = F_X0(:,iElem)
END DO
!$omp end parallel do
CALL VectorDotProduct(F_X0,F_X0,Norm2_F_X0)
! Prepare stuff for matrix vector multiplication
IF (Norm2_F_X0.LE.(1.E-12)**2*nElems) THEN ! do not iterate, as U is already the implicit solution
  Norm2_F_Xk=TINY(1.)
ELSE ! we need iterations
 Norm2_F_Xk=Norm2_F_X0
END IF
nInnerNewton=0
IF(Precond) CALL BuildMatrix(t,dt)
!-----------------------------------------------------------------------------------------------------------------------------------
! Newton iterations
!-----------------------------------------------------------------------------------------------------------------------------------
DO WHILE((Norm2_F_Xk.GT.Eps2Newton*Norm2_F_X0).AND. (nInnerNewton.LT.nNewtonIter))
  IF (nInnerNewton.EQ.0) THEN
    AbortCritGMRES=0.999
    Norm2_F_Xk_old=Norm2_F_Xk
  ELSE
    gammaA = gammaEW*(Norm2_F_Xk)/(Norm2_F_Xk_old)
    IF (gammaEW*AbortCritGMRES*AbortCritGMRES < 0.1) THEN
      gammaB = min(0.999,gammaA)
    ELSE
      gammaB = min(0.999, max(gammaA,gammaEW*AbortCritGMRES*AbortCritGMRES))
    ENDIF
    AbortCritGMRES = min(0.999,max(gammaB,0.5*SQRT(Eps2Newton)/SQRT(Norm2_F_Xk)))
    Norm2_F_Xk_old=Norm2_F_Xk
  END IF
  nInnerNewton=nInnerNewton+1
  CALL GMRES_M(time,dt,Alpha,Beta,-F_Xk,SQRT(Norm2_F_Xk),AbortCritGMRES,DeltaX)
  !$omp parallel do private(aElem)
  DO iElem=1,nElems
    aElem => Elems(iElem)%Elem
    DO iVar=1,NVAR
      Xk(iVar,iElem)=Xk(iVar,iElem)+DeltaX(iVar,iElem)
      aElem%cvar(iVar) = XK(iVar,iElem)
    END DO
    CALL ConsPrim(aElem%pvar, aElem%cvar)
  END DO
  !$omp end parallel do
  CALL FV_TimeDerivative(time, iter )
  !$omp parallel do private(aElem)
  DO iElem=1,nElems
    aElem => Elems(iElem)%Elem
    DO iVar=1,NVAR
      R_Xk(iVar,iElem)=aElem%u_t(iVar)
      F_Xk(iVar,iElem)=aElem%cvar(iVar)-Q(iVar,iElem)-alpha*dt*aElem%u_t(iVar)
    END DO
    !CALL ConsPrim(aElem%pvar, aElem%cvar)
  END DO
  !$omp end parallel do
  CALL VectorDotProduct(F_Xk,F_Xk,Norm2_F_Xk)
END DO
nNewtonIterGlobal=nNewtonIterGlobal+nInnerNewton
IF (nInnerNewton.EQ.nNewtonIter) THEN
  WRITE(*,*) Eps2Newton
  WRITE(*,*) ' Newton NOT converged with NEWTON ITERATIONS :', nInnerNewton
  WRITE(*,'(A22,E16.8)')   ' Norm / Norm_R0     : ',Norm2_F_XK / Norm2_F_X0
  STOP
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
CALL GlobalResidual(dt,res_iter)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ImplicitTimeStep

END MODULE MOD_ImplicitTimeStep
