MODULE MOD_LinearSolver
!===================================================================================================================================
! Contains the initialization of the Implicit global variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitLinearSolver
  MODULE PROCEDURE InitLinearSolver
END INTERFACE
INTERFACE FinalizeLinearSolver
  MODULE PROCEDURE FinalizeLinearSolver
END INTERFACE
INTERFACE VectorDotProduct
  MODULE PROCEDURE VectorDotProduct
END INTERFACE
INTERFACE GMRES_M
  MODULE PROCEDURE GMRES_M
END INTERFACE
INTERFACE BuildMatrix
  MODULE PROCEDURE BuildMatrix
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC::InitLinearSolver,FinalizeLinearSolver,VectorDotProduct,GMRES_M, BuildMatrix
!===================================================================================================================================

CONTAINS

SUBROUTINE InitLinearSolver()
!===================================================================================================================================
! Allocate global variable 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_LinearSolver_Vars
USE MOD_Readintools
USE MOD_Mesh_Vars,          ONLY: nElems,tElem,tSide,Elems
USE MOD_TimeDisc_Vars,      ONLY: implicit
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iElem,iSide,iSide2,NBElemID,ElemID
TYPE(tElem), POINTER:: aElem
TYPE(tSide), POINTER:: aSide,aSide2
!===================================================================================================================================
IF (implicit) THEN
  nKDim            = GETINT('nKDim','5')
  Eps2Newton       = GETREAL('EpsNewton','0.001')
  Eps2Newton       = Eps2Newton**2
  EpsGMRES         = GETREAL('EpsGMRES','0.001')
  rEps0            = SQRT(EPSILON(0.0))
  srEps0           = 1./rEps0
  nNewtonIter      = GETINT('nNewtonIter','20')
  nNewtonIterGlobal= 0
  nGMRESIterGlobal = 0
  nInnerNewton     = 0
  nInnerGMRES      = 0
  gammaEW          = GETREAL('gammaEW','0.9')
  ALLOCATE(  XK(1:NVAR,1:nElems) &
          ,R_XK(1:NVAR,1:nElems) )
  Precond          = GETLOGICAL('Precond','.FALSE.')
  IF(Precond)THEN
    ALLOCATE(Dinv(1:NVAR,1:NVAR,1:nElems))
    ALLOCATE(LowerUpper(1:NVAR,1:NVAR,1:4,1:nElems))
    ALLOCATE(ElemToElem(1:2,1:4,1:nElems))
    ElemToElem=0
    DO iElem = 1, nElems
      aElem => Elems(iElem)%Elem
      aSide => aElem%FirstSide
      iSide=1
      DO WHILE(ASSOCIATED(aSide))
        NBElemID=aSide%connection%Elem%ID
        IF (NBElemID .GT. 0 .AND. NBElemID .LE. nElems) THEN
          aSide2 => aSide%Connection%Elem%FirstSide
          iSide2=0
          DO WHILE(ASSOCIATED(aSide2))
            iSide2=iSide2+1
            ElemID=aSide%connection%Elem%ID
            IF((ElemID.GT.0).AND.(ElemID.LE.nElems))THEN
              IF(aSide2%connection%Elem%ID.EQ.iElem)THEN
                ElemToElem(1,iSide2,NBElemID)=iElem
                ElemToElem(2,iSide2,NBElemID)=iSide
              END IF
            END IF
            aSide2 => aSide2%nextelemside
          END DO
        END IF !
        iSide=iSide+1
        aSide=>aSide%nextelemSide
      END DO
    END DO
  END IF
  ImplicitInitIsDone=.TRUE.
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE InitLinearSolver


SUBROUTINE GMRES_M(t,dt,Alpha,Beta,B,Norm_B,AbortCrit,DeltaX)
!===================================================================================================================================
! Uses matrix free to solve the linear system
! Attention: We use DeltaX=0 as our initial guess
!            X0 is allready stored in U
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY: nElems
USE MOD_LinearSolver_Vars,      ONLY: EpsGMRES,Xk,nKDim,nGMRESIterGlobal,nInnerGMRES,Precond
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t,dt,Alpha,Beta,Norm_B
REAL,INTENT(IN)   :: B(1:NVAR,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT):: AbortCrit
REAL,INTENT(OUT)  :: DeltaX(1:NVAR,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: V (1:NVAR,1:nElems,1:nKDim)
REAL              :: W (1:NVAR,1:nElems)
REAL              :: Z (1:NVAR,1:nElems,1:nKDim)
REAL              :: R0(1:NVAR,1:nElems)
REAL              :: Gam(1:nKDim+1),C(1:nKDim),S(1:nKDim),H(1:nKDim+1,1:nKDim+1),Alp(1:nKDim)
REAL              :: Norm_R0,Resu,Temp,Bet
INTEGER           :: m,nn,o
!===================================================================================================================================

!print*,'Norm_R',Norm_B

AbortCrit=EpsGMRES*Norm_B
!Restart=0
R0=B
Norm_R0=Norm_B
DeltaX=0.
nInnerGMRES=0

! GMRES(m)  
V(:,:,1)=R0/Norm_R0
Gam(1)=Norm_R0
DO m=1,nKDim
  nInnerGMRES=nInnerGMRES+1
  ! Preconditioner
  IF(Precond)THEN
    CALL LUSGS_FD(t,dt,V(:,:,m),Z(:,:,m))
  ELSE
    Z(:,:,m)=V(:,:,m)
  END IF
  ! matrix vector
  CALL MatrixVector(t,dt,Alpha,Beta,Z(:,:,m),W)
  ! Gram-Schmidt
  DO nn=1,m
    CALL VectorDotProduct(V(:,:,nn),W,H(nn,m))
    W=W-H(nn,m)*V(:,:,nn)
  END DO !nn
  CALL VectorDotProduct(W,W,Resu)
  H(m+1,m)=SQRT(Resu)
  ! Givens Rotation
  DO nn=1,m-1
    Temp     =   C(nn)*H(nn,m) + S(nn)*H(nn+1,m)
    H(nn+1,m) = - S(nn)*H(nn,m) + C(nn)*H(nn+1,m)
    H(nn,m)   =   Temp
  END DO !nn
  Bet=SQRT(H(m,m)*H(m,m)+H(m+1,m)*H(m+1,m))
  S(m)=H(m+1,m)/Bet
  C(m)=H(m,m)/Bet 
  H(m,m)=Bet
  Gam(m+1)=-S(m)*Gam(m)
  Gam(m)=C(m)*Gam(m)
  IF ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) THEN !converge or max Krylov reached
    DO nn=m,1,-1
       Alp(nn)=Gam(nn) 
       DO o=nn+1,m
         Alp(nn)=Alp(nn) - H(nn,o)*Alp(o)
       END DO !o
       Alp(nn)=Alp(nn)/H(nn,nn)
    END DO !nn
    DO nn=1,m
      DeltaX=DeltaX+Alp(nn)*Z(:,:,nn)
    END DO !nn
    ! Preconditioner back
    !Z=DeltaX
    !CALL ApplyPrecond(Z,DeltaX)
    nGMRESIterGlobal=nGMRESIterGlobal+nInnerGMRES 
    RETURN
  ELSE ! no convergence, next iteration   ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim)) 
    V(:,:,m+1)=W/H(m+1,m)
  END IF ! ((ABS(Gam(m+1)).LE.AbortCrit) .OR. (m.EQ.nKDim))
END DO ! m 

WRITE(*,*) "GMRES not converged with GMRES iterations: ", nInnerGMRES
WRITE(*,'(A22,E16.8)')   ' Norm_R0            : ',ABS(Gam(1))
WRITE(*,'(A22,E16.8)')   ' Norm_R             : ',ABS(Gam(m))
STOP
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE GMRES_M


SUBROUTINE BuildMatrix(t,dt)
!===================================================================================================================================
! compute the global jacobian matrix by use of finite differences
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY:nElems,nElems,Elems,tElem,tSide,FirstElem,nBCSides,BCSides,Sides,nSides
USE MOD_LinearSolver_Vars,      ONLY:EpsGMRES,Xk,rEps0,srEps0,Dinv,LowerUpper,ElemToElem,R_XK
USE MOD_FluxCalculation,        ONLY:FluxJacobianFD
USE MOD_EOS
USE MOD_ExactFunc
USE MOD_Boundary,               ONLY:SetBCatSides,Boundary
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t,dt!,Alpha,Beta,Norm_B
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: D(1:NVAR,1:NVAR,1:nElems)      ! Diagonal Matrix
LOGICAL             :: OK_FLAG
TYPE(tElem), POINTER:: aElem
TYPE(tSide), POINTER:: aSide
INTEGER             :: iSide,iElem,iVar,iVar2,NBElemID,NBSideID
!===================================================================================================================================

! zero matrices
D=0.
LowerUpper=0.

! fd for diagonal element
!$omp parallel do private(aElem,aSide,iSide,NBElemID,NBSideID,iVar,iVar2,OK_FLAG)
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  DO iVar=1,NVAR
    CALL FluxJacobianFD(t,iElem,iVar)
    D(:,iVar,iElem) = aElem%u_t(:) 
    aSide => aElem%FirstSide
    iSide=1
    DO WHILE(ASSOCIATED(aSide))
      NBElemID=aSide%connection%Elem%ID
      IF (NBElemID .GT. 0 .AND. NBElemID .LE. nElems) THEN
        NBSideID=ElemToElem(2,iSide,iElem)
        LowerUpper(:,iVar,NBSideID,NBElemID) =  aSide%connection%Elem%u_t(:)
      END IF
       aSide=>aSide%nextelemSide
       iSide=iSide+1
    END DO
  END DO ! ivar
  DO iVar=1,NVAR
    DO iVar2=1,NVAR
      D(iVar,iVar2,iElem)=-dt*D(iVar,iVar2,iElem)
    END DO ! iVar2
    D(iVar,iVar,iElem)=D(iVar,iVar,iElem)+1.0
  END DO
  ! backup
  CALL M44INV(D(:,:,iElem), DINV(:,:,iElem), OK_FLAG) 
  IF (.not. OK_FLAG) THEN
    WRITE(*,*) " LUSGS D-Matrix is singular."
    WRITE(*,*) " ElemID", iElem
    WRITE(*,*) " STOP!"
    STOP
  END IF        
END DO
!$omp end parallel do
LowerUpper=-dt*LowerUpper

END SUBROUTINE BuildMatrix


SUBROUTINE LUSGS_FD(t,dt,B,DeltaX)
!===================================================================================================================================
! LUSGS preconditioner, uses precomputed Block-LUSGS
! B is old vector and DeltaX is preconditioned vector
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY:nElems,nElems,Elems,tElem,tSide,FirstElem,nBCSides,BCSides,Sides,nSides
USE MOD_LinearSolver_Vars,      ONLY:EpsGMRES,DINV,LowerUpper,ElemToElem
USE MOD_Equation_Vars,          ONLY:gamma
USE MOD_Boundary,               ONLY:SetBCatSides,Boundary
USE MOD_FV,                     ONLY:FV_TimeDerivative
USE MOD_Flux_lax_friedrichs
USE MOD_EOS
USE MOD_ExactFunc
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t,dt!,Alpha,Beta,Norm_B
REAL,INTENT(IN)   :: B(1:NVAR,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL,INTENT(INOUT):: AbortCrit
REAL,INTENT(OUT)  :: DeltaX(1:NVAR,1:nElems)        ! DeltaX
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: DeltaXStar(1:NVAR,1:nElems)    ! DeltaXStar
TYPE(tElem), POINTER:: aElem
TYPE(tSide), POINTER:: aSide
INTEGER             :: iVar,iElem,iSide,NBElemID,NBSideID
!===================================================================================================================================

!CALL BuildMatrix(t,dt)
  
! compute LU-SGS with FDs
DeltaX=0.
DeltaXStar=0.
! forward sweep
!$omp parallel do private(aElem,aSide,iSide,NBElemID)
DO iElem=1,nElems
  aElem => Elems(iElem)%Elem
  aSide => aElem%FirstSide
  iSide=1
  DO WHILE(ASSOCIATED(aSide))
    NBElemID=aSide%connection%Elem%ID
    !NBSideID=ElemToElem(2,iSide,iElem)
    IF (NBElemID .GT. 0 .AND. NBElemID .LT. iElem) THEN
      ! rotate state into normal direction
      DeltaXStar(:,iElem) = DeltaXStar(:,iElem)+MATMUL(LowerUpper(:,:,iSide,iElem),DeltaXStar(:,NBElemID))
    END IF
    aSide => aSide%nextElemSide
    iSide=iSide+1
  END DO
  ! Calculate DeltaXStar  
  DeltaXStar(:,iElem) = MATMUL(DINV(:,:,iElem), B(:,iElem) - DeltaXStar(:,iElem)) 
END DO ! iElem
!$omp end parallel do
 
! backwards sweep
!$omp parallel do private(aElem,aSide,iSide,NBElemID)
DO iElem=nElems,1,-1
  aElem => Elems(iElem)%Elem
  ! backward sweep
  aSide => aElem%FirstSide
  iSide=1
  DO WHILE(ASSOCIATED(aSide))
    NBElemID=aSide%connection%Elem%ID
    IF((NBElemID.GT. iElem) .AND.(NBElemID .LE. nElems)) THEN
      DeltaX(:,iElem) = DeltaX(:,iElem) + MATMUL(LowerUpper(:,:,iSide,iElem),DeltaX(:,NBElemID))
    END IF
    aSide => aSide%nextElemSide
    iSide = iSide+1
  END DO
  ! Calculate DeltaX
  DeltaX(:,iElem) = DeltaXStar(:,iElem)- MATMUL(DINV(:,:,iElem),DeltaX(:,iElem))
END DO ! iElem
!$omp end parallel do
 
END SUBROUTINE LUSGS_FD


SUBROUTINE LUSGS_SLOW(t,dt,B,DeltaX)
!===================================================================================================================================
! compute LU-SGS, contains build in matrix 
! whole system matrix is computed via finite difference, used for debugging
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,              ONLY:nElems,nElems,Elems,tElem,tSide,FirstElem,nBCSides,BCSides,Sides,nSides
USE MOD_LinearSolver_Vars,      ONLY:EpsGMRES,Xk,R_XK
USE MOD_Equation_Vars,          ONLY:gamma
USE MOD_Boundary,               ONLY:SetBCatSides,Boundary
USE MOD_FV,                     ONLY: FV_TimeDerivative
USE MOD_Flux_lax_friedrichs
USE MOD_EOS
USE MOD_ExactFunc
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t,dt
REAL,INTENT(IN)   :: B(1:NVAR,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: DeltaX(1:NVAR,1:nElems)        ! DeltaX
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: D(1:NVAR,1:NVAR,1:nElems)      ! Diagonal Matrix
REAL                :: DINV_SLOW(1:NVAR,1:NVAR,1:nElems)   ! Diagonal Matrix Inverse
REAL                :: DeltaXStar(1:NVAR,1:nElems)    ! DeltaXStar
REAL                :: DeltaXStar_tmp(1:NVAR,1:nElems)! DeltaXStar
REAL                :: DeltaX_tmp(1:NVAR,1:nElems)    ! DeltaXStar
LOGICAL             :: OK_FLAG
TYPE(tElem), POINTER:: aElem,aElem2
TYPE(tSide), POINTER:: aSide
REAL                 :: dRdU(1:4*nElems,1:4*nElems)
REAL                 :: rEps0,srEps0
INTEGER              :: iVar,iElem,iSide,r,s,iVar2,iElem2,NBElemID
REAL                 :: R0  (1:NVAR,1:nElems)
!===================================================================================================================================


rEps0=SQRT(EPSILON(0.0))
srEps0=1./rEps0
  

! whole system matrix - required to get old values for preconditoner
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  !U = Xk+EpsFD*V
  aElem%cvar(:) = XK(:,iElem)
  CALL ConsPrim(aElem%pvar, aElem%cvar)
END DO ! iElem

dRdU=0.
s=0
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  DO iVar=1,NVAR
    s=s+1
    ! update of conservative variables
    aElem%cvar(iVar) = aElem%cvar(iVar)+reps0
    ! obtain thermodynamic consistent state
    CALL ConsPrim(aElem%pvar, aElem%cvar)
    CALL FV_TimeDerivative(t, 0) 
    ! reset conservative variables
    aElem%cvar(iVar) = aElem%cvar(iVar)-reps0
    ! obtain thermodynamic consistent state
    CALL ConsPrim(aElem%pvar, aElem%cvar)
    ! diagonal entry
    r=(iElem-1)*NVAR
    DO iVar2=1,NVAR
      r=r+1
      dRdU(r,s) = dRdU(r,s)+(aElem%u_t(iVar2)-R_XK(iVar2,iElem))*sreps0
    END DO !iiVar
    ! neighbor diagonal
    aSide => aElem%FirstSide
    DO WHILE(ASSOCIATED(aSide))
      IF((aSide%connection%Elem%ID .GT. 0) .AND.(aSide%connection%Elem%ID .LE. nElems)) THEN
        iElem2=aSide%connection%Elem%ID
        aElem2 => Elems(iElem2)%Elem
        r=(iElem2-1)*NVAR
        DO iVar2=1,NVAR
          r=r+1
          dRdU(r,s) = dRdU(r,s)+(aElem2%u_t(iVar2)-R_XK(iVar2,iElem2))*sreps0
        END DO !iiVar
      END IF
      aSide => aSide%nextElemSide
    END DO
  END DO ! iVar
END DO


DO r=1,4*nElems
  DO s=1,4*nElems
    dRdU(r,s)=-dt*dRdU(r,s)
  END DO ! s
  dRdU(r,r)=dRdU(r,r)+1.0
END DO ! r

! build dinv
DO iElem=1,nElems
  r=(iElem-1)*NVAR
  D(1:NVAR,1:NVAR,iElem) = dRdU(r+1:r+NVAR,r+1:r+NVAR)
  ! Invert D_ii Matrix
  CALL M44INV(D(:,:,iElem), DINV_SLOW(:,:,iElem), OK_FLAG) 
  IF (.not. OK_FLAG) THEN
    WRITE(*,*) " LUSGS D-Matrix is singular."
    WRITE(*,*) " ElemID", iElem
    WRITE(*,*) " STOP!"
    STOP
  END IF        
END DO 

!nInnerGMRES=0
DeltaX=0.
DeltaXStar=0.
DeltaX_tmp = 0.
DeltaXStar_tmp = 0.  
! forward sweep
! calculation of deltaXstar
DO iElem=1,nElems
  aElem => Elems(iElem)%Elem

  aSide => aElem%FirstSide
  DO WHILE(ASSOCIATED(aSide))
    NBElemID=aSide%connection%Elem%ID
    IF (NBElemID .GT. 0 .AND. NBElemID .LT. iElem) THEN
        r=(iElem-1)*NVAR
        s=(NBElemID-1)*NVAR
        DeltaXStar_tmp(:,iElem) = DeltaXStar_tmp(:,iElem) &
                                +MATMUL(dRdU(r+1:r+NVAR,s+1:s+NVAR),DeltaXStar(:,NBElemID))
    END IF
    aSide => aSide%nextElemSide
  END DO
  ! Calculate DeltaXStar  
  DeltaXStar(:,iElem) = MATMUL(DINV_SLOW(:,:,iElem), B(:,iElem) - DeltaXStar_tmp(:,iElem)) 
END DO ! iElem

! backwards sweep
! modified to caputure BC elems
DO iElem=nElems,1,-1
  aElem => Elems(iElem)%Elem
  ! backward sweep
  aSide => aElem%FirstSide
  DO WHILE(ASSOCIATED(aSide))
    NBElemID=aSide%connection%Elem%ID
    IF((NBElemID.GT. iElem) .AND.(NBElemID .LE. nElems)) THEN
      r=(iElem-1)*NVAR
      s=(NBElemID-1)*NVAR
      DeltaX_tmp(:,iElem) = DeltaX_tmp(:,iElem) &
                          +MATMUL(dRdU(r+1:r+NVAR,s+1:s+NVAR),DeltaX(:,NBElemID))
    END IF
    aSide => aSide%nextElemSide
  END DO
  ! Calculate DeltaX
  DeltaX(:,iElem) = DeltaXStar(:,iElem)- MATMUL(DINV_SLOW(:,:,iElem),DeltaX_tmp(:,iElem))
END DO ! iElem
 
END SUBROUTINE LUSGS_SLOW


FUNCTION FLUXJACOBIAN(PVAR,NormVec)
!===================================================================================================================================
! function to compute the flux jacobian
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Equation_Vars, ONLY : gamma,gamma1,gamma2,gamma1q
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: PVAR(NVAR)
REAL,INTENT(IN)                :: NormVec(1:2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(4,4)            :: FLUXJACOBIAN
! CHARACTER(LEN=255)             ::  ! the complete timestamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                           ::phi,grosV,a_3,a_2,a_1,E_loc
!===================================================================================================================================

! Parameters

!e_r  = gamma1q * p_r + 0.5 * rho_r * &
!       (v1_r * v1_r + v2_r * v2_r)

!E_loc = gamma * PVAR(P) + 0.5*PVAR(RHO)*(PVAR(V1)**2+PVAR(V2)**2)
E_loc = gamma1q * PVAR(P) + 0.5*PVAR(RHO)*(PVAR(V1)**2+PVAR(V2)**2)
E_loc = E_loc/pvar(RHO)
phi = 0.5*(gamma1)*(PVAR(V1)**2+PVAR(V2)**2)
grosV = NormVec(X_DIR)*PVAR(V1)+NormVec(Y_DIR)*PVAR(V2)
a_3 = gamma2
a_2 = gamma1
a_1 = gamma*E_loc-phi

! compute the flux jacobian

FLUXJACOBIAN(1,1) = 0
FLUXJACOBIAN(1,2) = NormVec(X_DIR)
FLUXJACOBIAN(1,3) = NormVec(Y_DIR)
FLUXJACOBIAN(1,4) = 0

FLUXJACOBIAN(2,1) = NormVec(X_DIR)*phi-PVAR(V1)*grosV
FLUXJACOBIAN(2,2) = grosV-a_3*NormVec(X_DIR)*PVAR(V1)
FLUXJACOBIAN(2,3) = NormVec(Y_DIR)*PVAR(V1)-a_2*NormVec(X_DIR)*PVAR(V2)
FLUXJACOBIAN(2,4) = a_2*NormVec(X_DIR)

FLUXJACOBIAN(3,1) = NormVec(Y_DIR)*phi-PVAR(V1)*grosV
FLUXJACOBIAN(3,2) = NormVec(X_DIR)*PVAR(V2)-a_2*NormVec(Y_DIR)*PVAR(V1)
FLUXJACOBIAN(3,3) = grosV-a_3*NormVec(X_DIR)*PVAR(V2)
FLUXJACOBIAN(3,4) = a_2*NormVec(Y_DIR)

FLUXJACOBIAN(4,1) = grosV*(phi-a_1)
FLUXJACOBIAN(4,2) = NormVec(X_DIR)*a_1-a_2*PVAR(V1)*grosV
FLUXJACOBIAN(4,3) = NormVec(Y_DIR)*a_1-a_2*PVAR(V2)*grosV
FLUXJACOBIAN(4,4) = gamma*grosV
!-----------------------------------------------------------------------------------------------------------------------------------
END FUNCTION FLUXJACOBIAN

SUBROUTINE  M44INV(A, AINV, OK_FLAG)
!===================================================================================================================================
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  M44INV  -  Compute the inverse of a 4x4 matrix.
!
!  A       = input 4x4 matrix to be inverted
!  AINV    = output 4x4 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!             function to compute the flux jacobian
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, DIMENSION(4,4), INTENT(IN)  :: A
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, DIMENSION(4,4), INTENT(OUT) :: AINV
LOGICAL, INTENT(OUT) :: OK_FLAG
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL, PARAMETER :: EPS = 1.0D-10
REAL :: DET
REAL, DIMENSION(4,4) :: COFACTOR

!===================================================================================================================================

DET =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)   &
      -A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))         &
      +A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)  &
      -A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))         &
      +A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

IF (ABS(DET) .LE. EPS) THEN
  AINV = 0.0D0
  OK_FLAG = .FALSE.
  RETURN
END IF

COFACTOR(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
COFACTOR(1,2) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
COFACTOR(1,3) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
COFACTOR(1,4) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
COFACTOR(2,1) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
COFACTOR(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
COFACTOR(2,3) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
COFACTOR(2,4) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
COFACTOR(3,1) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
COFACTOR(3,2) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
COFACTOR(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
COFACTOR(3,4) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
COFACTOR(4,1) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
COFACTOR(4,2) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
COFACTOR(4,3) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
COFACTOR(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

AINV = TRANSPOSE(COFACTOR) / DET

OK_FLAG = .TRUE.

END SUBROUTINE M44INV

SUBROUTINE MatrixVector(t,dt,Alpha,Beta,V,Resu)
!===================================================================================================================================
! Computes Matrix Vector Product using the spatiall operator and finite difference approach, see Imperator.pdf
! Computes resu=A*v 
! A is operator at linearization state xk (Newton iteration)
! Important: needs definition of epsFD before calling subroutine
!            needs definition of xk before calling subroutine
!            needs computation of R_xk before calling subroutine
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_LinearSolver_Vars,  ONLY: Xk,R_Xk,rEps0,iterGlobal
USE MOD_FV,                 ONLY: FV_TimeDerivative
USE MOD_Mesh_Vars,          ONLY: nElems,Elems
USE MOD_Mesh_Vars,          ONLY: tElem
USE MOD_EoS,                ONLY: ConsPrim
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: t,Alpha,Beta,dt
REAL,INTENT(IN)   :: V   (1:NVAR,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu(1:NVAR,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: EpsFD
TYPE(tElem), POINTER    :: aElem
INTEGER           :: iElem
!===================================================================================================================================

! needed for FD matrix vector approximation
CALL VectorDotProduct(V,V,EpsFD)
EpsFD= rEps0/SQRT(EpsFD)

!$omp parallel do private(aElem)
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  !U = Xk+EpsFD*V
  aElem%cvar(:) = XK(:,iElem)+EpsFD*V(:,iElem)
  CALL ConsPrim(aElem%pvar, aElem%cvar)
END DO
!$omp end parallel do
CALL FV_TimeDerivative(t, iterGlobal ) 
!$omp parallel do private(aElem)
DO iElem = 1, nElems
  aElem => Elems(iElem)%Elem
  Resu(:,iElem) = V(:,iElem) - Alpha*dt*(aElem%u_t(:) - R_Xk(:,iElem))/EpsFD
END DO
!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE MatrixVector

SUBROUTINE VectorDotProduct(A,B,Resu)
!===================================================================================================================================
! Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,          ONLY: nElems
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: A(1:NVAR,1:nElems)
REAL,INTENT(IN)   :: B(1:NVAR,1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iVar,iElem
!===================================================================================================================================

Resu=0.
!$omp parallel do reduction(+:Resu)
DO iElem=1,nElems
   DO iVar=1,NVAR
     Resu=Resu + A(iVar,iElem)*B(iVar,iElem)
  END DO
END DO
!$omp end parallel do
!-----------------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE VectorDotProduct

SUBROUTINE FinalizeLinearSolver()
!===================================================================================================================================
! Deallocate global variable U (solution) and Ut (dg time derivative).
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
!===================================================================================================================================
!ImplicitInitIsDone = .FALSE.
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE FinalizeLinearSolver

FUNCTION IDENT(A)
!===================================================================================================================================
! function to compute identity matrix
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: A
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(A,A)            :: IDENT
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                           :: phi,grosV,a_3,a_2,a_1
INTEGER                        :: i,j
!===================================================================================================================================
! compute identity matrix
   DO i = 1,A
       DO j = 1,A
           IDENT(i,j) = 0.
           IF (i .eq. j) IDENT(i,j ) = 1.
       END DO
   END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END FUNCTION IDENT

END MODULE MOD_LinearSolver
