MODULE MOD_Mesh
!===================================================================================================================================
! Mesh Readin
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide 
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitMesh
   MODULE PROCEDURE InitMesh
END INTERFACE
INTERFACE CreateMesh
   MODULE PROCEDURE CreateMesh
END INTERFACE
INTERFACE CreateCartMesh
   MODULE PROCEDURE CreateCartMesh
END INTERFACE
INTERFACE CreateElementInfo
  MODULE PROCEDURE CreateElementInfo
END INTERFACE
INTERFACE CreateSideInfo
  MODULE PROCEDURE CreateSideInfo
END INTERFACE
INTERFACE CreateReconstructionInfo
  MODULE PROCEDURE CreateReconstructionInfo
END INTERFACE
INTERFACE ConnectPeriodicBC
  MODULE PROCEDURE ConnectPeriodicBC
END INTERFACE
INTERFACE SortSideList
   MODULE PROCEDURE SortSideList
END INTERFACE
INTERFACE PartitionSideList
   MODULE PROCEDURE PartitionSideList
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE tSideList
  INTEGER              :: Node(2)
  INTEGER              :: BC = 0
  TYPE(tSide), POINTER :: Side
  LOGICAL              :: isRotated = .FALSE.
END TYPE tSideList
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: InitMesh, CreateElementInfo, CreateSideInfo, CreateReconstructionInfo
!-----------------------------------------------------------------------------------------------------------------------------------
PRIVATE :: CreateMesh,CreateCartMesh, SortSideList, PartitionSideList, ConnectPeriodicBC
!===================================================================================================================================

CONTAINS

SUBROUTINE InitMesh()
!===================================================================================================================================
! Init Mesh, calls read meash
!===================================================================================================================================
! MODULES
USE MOD_Globals 
USE MOD_Mesh_Vars    ,ONLY:dxRef,totalArea_q,nElems,GridFile 
USE MOD_Output_Vars  ,ONLY:iVisuProg,strOutFile 
USE MOD_Timedisc_Vars,ONLY: Restart
USE MOD_Output_cgns  ,ONLY:CGNS_WriteMesh
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

! Initialize Mesh
WRITE(*,*) "-[Mesh]-----------------------------------------------------"
CALL ReadMesh()
GridFile =  TRIM(strOutFile) // '_mesh.cgns'
CALL CreateMesh()
IF (iVisuProg == CGNS) THEN
  IF (.NOT.Restart) THEN
    ! Write CGNS Mesh file
    CALL CGNS_WriteMesh()
  END IF
END IF
WRITE(*,*) ' ... done.'
! Average Element Size
dxRef = SQRT(1./(totalArea_q*nElems))
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE InitMesh


SUBROUTINE readMesh()
!===================================================================================================================================
! read in mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Readin      ,ONLY:getFreeIOUnit
USE MOD_Readintools ,ONLY:GETINT,GETSTR,GETREALARRAY,GETINTARRAY
USE MOD_Mesh_Vars   ,ONLY:xMax,xMin,yMax,yMin
USE MOD_Mesh_Vars   ,ONLY:CartMesh,Unst
USE MOD_Mesh_Vars   ,ONLY:MeshType,strMeshFile,strMeshFormat
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: Unit
INTEGER             :: i, iSide
REAL                :: X0(2), DX(2)
!===================================================================================================================================

CALL getFreeIOUnit(100,Unit)
MeshType = GETINT('Meshtype','1')
SELECT CASE (MeshType)
CASE(UNST)
! Unstructured grid: Read Mesh file
  WRITE(*,'(a)') '|   Mesh type is UNSTRUCTURED'
  strMeshFormat = GETSTR('MeshFormat')
  WRITE(*,'(a,a)') '   Mesh file format is ', TRIM(strMeshFormat)
  strMeshFile = GETSTR('Meshfile')
  WRITE(*,'(a,a)') '   Mesh file: ', TRIM(strMeshFile) // TRIM(strMeshFormat)
  ! Set CGNS grid file name to be equal to the .mesh file name
  strMeshFile = TRIM(strMeshFile)//strMeshFormat
CASE(CART)
! Cartesian Grid generator: Read required data
  WRITE(*,*) '|   Mesh type is CARTESIAN'
  CartMesh%imax = GETINT('nElemsX')
  CartMesh%jmax = GETINT('nElemsY')
  X0 = GETREALARRAY('X0',2)
  xmin = X0(1)
  ymin = X0(2)
  DX = GETREALARRAY('Xmax',2)
  xmax = DX(1)
  ymax = DX(2)
  !WRITE(*,'(a,F6.2,a,F6.2,a,F6.2,a,F6.2,a)') '   Grid extents: (', xmin, ',', ymin, ') x (' , &
  !                                         xmax, ',', ymax, ')'
! Boundary Conditions
  CartMesh%nBC(:)=GETINTARRAY('nBCsegments',4)
  DO iSide = 1, 4
    IF(CartMesh%nBC(iSide) == 1) THEN
      CartMesh%BCType(iSide, 1) = GETINT('MeshBCtype')
      CartMesh%BCRange(iSide, 1, 1) = 1
      IF (MOD(iSide,2) == 0) THEN
        CartMesh%BCRange(iSide, 1, 2) = CartMesh%jmax
      ELSE
        CartMesh%BCRange(iSide, 1, 2) = CartMesh%imax
      END IF
    ELSE
      DO i = 1, CartMesh%nBC(iSide)
        CartMesh%BCType(iSide, i) = GETINT('MeshBCtype')
        CartMesh%BCRange(iSide, i, 1) = GETINT('BCstart')
        CartMesh%BCRange(iSide, i, 2) = GETINT('BCend')
       END DO
    END IF
  END DO
CASE DEFAULT
  PRINT*, 'Error: MeshType is only unstructured(=0) or cartesian(=1)'
  PRINT*, 'stop in readin/readin.f90'
  STOP
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE readMesh


SUBROUTINE CreateMesh()
!===================================================================================================================================
! creat cartesian, structured mesh
! read in of all supported mesh types:
! *msh, *emc2, *cgns
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Output
USE MOD_Mesh_Readin            ,ONLY:ReadEMC2,ReadGmsh,ReadGridgenCGNS
USE MOD_Boundary_Vars          ,ONLY:tBoundary
USE MOD_Mesh_Vars              ,ONLY:tElem,tSide,tPureSidePtr,tNode,tNodePtr
USE MOD_Mesh_Vars              ,ONLY:nTrias,nQuads,nBCSides,nElems,nNodes,nSides,nInnerSides
USE MOD_Mesh_Vars              ,ONLY:Sides,Elems,BCSides,FirstBCSide,firstElem,FirstNode,FirstSide
USE MOD_Mesh_Vars              ,ONLY:totalArea_q,xMax,xMin,yMax,yMin
USE MOD_Mesh_Vars              ,ONLY:MeshType,strMeshFormat,strMeshFile
USE MOD_Boundary_Vars          ,ONLY:FirstBC
USE MOD_Initialcondition_Vars  ,ONLY:nDomains
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
! Variables for MESH I/O
INTEGER                  :: nVertices, nBCEdges,nZones
INTEGER, POINTER         :: Tria(:,:), Quad(:,:), BCEdge(:,:),ZoneConnect(:)
REAL, POINTER            :: Vertex(:,:)
! Variables for Nodelist
INTEGER                   :: iNode
TYPE(tNode), POINTER      :: aNode => NULL()
TYPE(tNodePtr), POINTER   :: VertexPtr(:)
! Variables for Elementlist
INTEGER                   :: iElem, iTria, iQuad, iSide, iSidePtr
INTEGER                   :: iNode1, iNode2
TYPE(tElem),POINTER       :: aElem, prevElem
TYPE(tSide),POINTER       :: aSide, bSide, cSide
TYPE(tPureSidePtr),POINTER:: aBCSide
! Variables for Sidelist
TYPE(tSideList), POINTER  :: SideList(:)
REAL                      :: temp
INTEGER                   :: BCID, BCType, iBC, BC(4)
! Variables for extended MESH data
TYPE(tBoundary), POINTER  :: aBC
!===================================================================================================================================

nTrias = 0
nQuads = 0
! Create cartesian Mesh or read unstrFuctured Mesh from file
IF (MeshType .eq. CART) THEN
! Extract Boundary Conditions
  iBC = 4
  aBC => FirstBC
  DO WHILE(ASSOCIATED(aBC))
    BC(iBC) = aBC%BCType * 100 + aBC%BCID
    iBC = iBC - 1
    aBC => aBC%nextBC
  END DO
! Create cartesian Mesh
  CALL CreateCartMesh(Vertex, nVertices,                      &
                      BCEdge, nBCEdges,                       &
                      Quad, nQuads                            )
  ALLOCATE(ZoneConnect(nVertices))
  ZoneConnect(:)=0  ! because of Gridgen Zone Connectivity
ELSE
  ! read emc2,gmsh or cgns file
  SELECT CASE (strMeshFormat)
  CASE('.mesh')
    CALL ReadEMC2(strMeshFile,        &
                  Vertex, nVertices,  &
                  BCEdge, nBCEdges,   &
                  Tria, nTrias,       &
                  Quad, nQuads        )
    ALLOCATE(ZoneConnect(nVertices))
    ZoneConnect(:)=0 !because of Gridgen Zone Connectivity

    CASE('.msh')
     WRITE(*,*) '    reading Gmsh File...'
     CALL ReadGmsh(strMeshFile,        &
                   Vertex, nVertices,  &
                   BCEdge, nBCEdges,   &
                   Tria, nTrias,       &
                   Quad, nQuads        )
     ALLOCATE(ZoneConnect(nVertices))
     ZoneConnect(:)=0 !because of Gridgen Zone Connectivity

  CASE('.cgns')
    WRITE(*,*)'    read Gridgen CGNS file...'
    CALL ReadGridgenCGNS(strMeshFile,     &
                  Vertex, nVertices,  &
                  BCEdge, nBCEdges,   &
                  Tria, nTrias,       &
                  Quad, nQuads,       &
                  ZoneConnect,nZones  )
    IF ((nZones .NE. nDomains).AND.(nDomains > 1)) THEN
      WRITE(*,'(a,I3)')'number of Domains in Gridgen-File :', nZones
      WRITE(*,'(a,I3)')'and number of Domains in ParameterFile :', nDomains
      WRITE(*,*)'are contradictory!' 
      STOP
    END IF 
    WRITE(*,*)'    ... done. '
  END SELECT
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Generate Mesh information
nElems = nTrias + nQuads
nInnerSides = (nTrias*3 + nQuads*4 - nBCEdges) / 2
nBCSides    = nBCEdges
nSides      = nInnerSides + nBCEdges
nNodes      = 0 !nVertices
!-----------------------------------------------------------------------------------------------------------------------------------
! Nodes
!-----------------------------------------------------------------------------------------------------------------------------------
ALLOCATE(VertexPtr(nVertices))
DO iNode=1,nVertices
   NULLIFY(VertexPtr(iNode)%Node)
END DO
xMin = Vertex(1, X_DIR)
xMax = Vertex(1, X_DIR)
yMax = Vertex(1, Y_DIR)
yMin = Vertex(1, Y_DIR)
! Generate first node in Mesh
FirstNode => NULL()
! loop backwards over all vertices in order to create nodes
DO iNode = nVertices, 1, -1
 ! Zone Connection in Gridgen with Nodes
 IF (ZoneConnect(iNode) .EQ. 0) THEN! inner node of zone 
 ! allocate memory for the new node and assign data to it
   ALLOCATE(aNode)
   aNode%ID   = iNode
   aNode%x(1) = Vertex(iNode, X_DIR)
   aNode%x(2) = Vertex(iNode, Y_DIR)
   xmin  = MIN(aNode%x(1), xMin)
   xmax  = MAX(aNode%x(1), xMax)
   ymin  = MIN(aNode%x(2), yMin)
   ymax  = MAX(aNode%x(2), yMax)
 ! Link nodes to each other
   aNode%nextNode => FirstNode
   FirstNode => aNode
 ! save pointer to new node in Pointer array
   VertexPtr(iNode)%Node => aNode
 ! count nodes
   nNodes=nNodes+1
 ELSE
   IF (ASSOCIATED(VertexPtr(ZoneConnect(iNode))%Node)) THEN
     VertexPtr(iNode)%Node => VertexPtr(ZoneConnect(iNode))%Node ! already built
   ELSE
   ! allocate memory for the new node and assign data to it
     ALLOCATE(aNode)
     aNode%ID   = iNode
     aNode%x(1) = Vertex(iNode, X_DIR)
     aNode%x(2) = Vertex(iNode, Y_DIR)
     xmin  = MIN(aNode%x(1), xMin)
     xmax  = MAX(aNode%x(1), xMax)
     ymin  = MIN(aNode%x(2), yMin)
     ymax  = MAX(aNode%x(2), yMax)
   ! Link nodes to each other
     aNode%nextNode => FirstNode
     FirstNode => aNode
   ! save pointer to new node in Pointer array
     VertexPtr(iNode)%Node => aNode
     VertexPtr(ZoneConnect(iNode))%Node => VertexPtr(iNode)%Node
   ! count nodes
     nNodes=nNodes+1
   END IF
 END IF
END DO
DEALLOCATE(Vertex,ZoneConnect)
!-----------------------------------------------------------------------------------------------------------------------------------
! Elements and side pointers
!-----------------------------------------------------------------------------------------------------------------------------------
! Allocate SideList for later use
ALLOCATE(SideList(nSides*2))
! Initialize Element List and Element Counter
FirstElem => NULL()
iElem = 0
iSidePtr = 0
!-----------------------------------------------------------------------------------------------------------------------------------
! Loop over all triangles
DO iTria = 1, nTrias
! Element Counter
  iElem = iElem + 1
! allocate memory and fill it with data
  ALLOCATE(aElem)
! save Element ID and Domain
  aElem%ID     = iElem
  aElem%domain = Tria(iTria, 4)
  aElem%ElemType   = 3
! connect elemlist
  IF (.NOT.ASSOCIATED(FirstElem)) THEN
    FirstElem => aElem
    prevElem => aElem
  ELSE
    prevElem%nextElem => aElem
    prevElem => aElem
  END IF
  aElem%nextElem => NULL()
! Save Element Nodes (these nodes will alway remain identical and are only used for I/O)
  ALLOCATE(aElem%NodeArray(3))
  DO iNode = 1, 3
    aElem%NodeArray(iNode)%Node => VertexPtr(Tria(iTria, iNode))%Node
  END DO
! Compute Barycenter of Element
  aElem%Bary(:) = 0.
  DO iNode = 1, 3
    aElem%Bary(:) = aElem%Bary(:) + aElem%NodeArray(iNode)%Node%x(:) !old Vertex(Quad(iQuad, iNode), :)
  END DO
  aElem%Bary(:) = aElem%Bary(:) / 3.
! Build sides for the element
  aElem%firstSide => NULL()
  DO iSide = 1, 3
  ! initialize Side
    ALLOCATE(aSide)
    aSide%connection   => NULL()
    aSide%nextElemSide => NULL()
    aSide%nextSide     => NULL()
    aSide%BC           => NULL()
    aSide%Node(1)%node => aElem%NodeArray(iSide)%Node
    aSide%Node(2)%node => aElem%NodeArray(MOD(iSide,3)+1)%Node
    aSide%Elem         => aElem
    iSidePtr = iSidePtr + 1
  ! link Side to local element list
    aSide%nextElemSide => aElem%firstSide
    aElem%firstSide => aSide
  ! write Side into global array which will be used for generating the connectivity
    iNode1=aElem%NodeArray(iSide)%Node%ID
    iNode2=aElem%NodeArray(MOD(iSide,3)+1)%Node%ID
    SideList(iSidePtr)%Node(1) = MIN(iNode1, iNode2)
    SideList(iSidePtr)%Node(2) = MAX(iNode1, iNode2)
    SideList(iSidePtr)%BC = 0
    IF (MIN(iNode1, iNode2) == iNode2) THEN
      SideList(iSidePtr)%isRotated = .TRUE.
    END IF
    SideList(iSidePtr)%Side => aSide
    ! next side ...
  END DO
END DO
IF (nTrias .GT. 0 ) DEALLOCATE(Tria)
!-----------------------------------------------------------------------------------------------------------------------------------
! Loop over all quadrangles
DO iQuad = 1, nQuads
! Element Counter
  iElem = iElem + 1
! allocate memory and fill it with data
  ALLOCATE(aElem)
! save Element ID and Domain
  aElem%ID     = iElem
  aElem%domain = Quad(iQuad, 5)
  aElem%ElemType   = 4
! connect elemlist
  IF (.NOT.ASSOCIATED(FirstElem)) THEN
    FirstElem => aElem
    prevElem => aElem
  ELSE
    prevElem%nextElem => aElem
    prevElem => aElem
  END IF
  aElem%nextElem => NULL()
! Save Element Nodes (these nodes will alway remain identical and are only used for I/O)
  ALLOCATE(aElem%NodeArray(4))
  DO iNode = 1, 4
    aElem%NodeArray(iNode)%Node => VertexPtr(Quad(iQuad, iNode))%Node
  END DO
! Compute Barycenter of Element
  aElem%Bary(:) = 0.
  DO iNode = 1, 4
    aElem%Bary(:) = aElem%Bary(:) + aElem%NodeArray(iNode)%Node%x !Vertex(Quad(iQuad, iNode), :)
  END DO
  aElem%Bary(:) = aElem%Bary(:) * 0.25
! Build sidepointers for the element
  aElem%firstSide => NULL()
  DO iSide = 1, 4
  ! initialize SidePtr
    ALLOCATE(aSide)
    aSide%connection   => NULL()
    aSide%nextElemSide => NULL()
    aSide%nextSide     => NULL()
    aSide%BC           => NULL()
    aSide%Node(1)%node => aElem%NodeArray(iSide)%Node
    aSide%Node(2)%node => aElem%NodeArray(MOD(iSide,4)+1)%Node
    aSide%Elem         => aElem
    iSidePtr = iSidePtr + 1
  ! link Side to local element list
    aSide%nextElemSide => aElem%firstSide
    aElem%firstSide => aSide
  ! write SidePtr into global array which will be used for generating the connectivity
    iNode1=aElem%NodeArray(iSide)%Node%ID
    iNode2=aElem%NodeArray(MOD(iSide,4)+1)%Node%ID
    SideList(iSidePtr)%Node(1) = MIN(iNode1, iNode2)
    SideList(iSidePtr)%Node(2) = MAX(iNode1, iNode2)
    SideList(iSidePtr)%BC = 0
    IF (MIN(iNode1, iNode2) == iNode2) THEN
      SideList(iSidePtr)%isRotated = .TRUE.
    END IF
    SideList(iSidePtr)%Side => aSide
  ! next side ...
  END DO
END DO
IF (nQuads .GT. 0 ) DEALLOCATE(Quad)
!-----------------------------------------------------------------------------------------------------------------------------------
! Sides and connectivity
!-----------------------------------------------------------------------------------------------------------------------------------
! Save all BCEdges into the SideList Array
DO iSide = 1, nBCSides
  iSidePtr = iSidePtr + 1
! initialize Ghost Side
  ALLOCATE(aSide)
  aSide%connection   => NULL()
  aSide%nextElemSide => NULL()
  aSide%nextSide     => NULL()
  aSide%Node(1)%node => VertexPtr(BCEdge(iSide,1))%Node
  aSide%Node(2)%node => VertexPtr(BCEdge(iSide,2))%Node
! Initialize Ghost Elem
  ALLOCATE(aElem)
  aSide%Elem => aElem
! Save nodes of BC Sides
  iNode1 = VertexPtr(BCEdge(iSide,1))%Node%ID
  iNode2 = VertexPtr(BCEdge(iSide,2))%Node%ID
! Rotate nodes if necessary
  SideList(iSidePtr)%Node(1) = MIN(iNode1, iNode2)
  SideList(iSidePtr)%Node(2) = MAX(iNode1, iNode2)
  IF (MIN(iNode1, iNode2) == iNode2) THEN
    SideList(iSidePtr)%isRotated = .TRUE.
  END IF
! Save Boundary Condition
  BCType = BCEdge(iSide,3)/100
  BCID   = MOD(BCEdge(iSide,3), 100)
  aBC => FirstBC
  DO WHILE (ASSOCIATED(aBC))
    IF ((aBC%BCType == BCType).AND.(aBC%BCID == BCID)) THEN
      aSide%BC => aBC
      EXIT
    END IF
    aBC => aBC%nextBC
  END DO
  IF ( (BCType .NE. 0 ) .AND. (.not.ASSOCIATED(aBC))) THEN
     WRITE(*,*)'Boundary condition',BCTYPE*100 + BCID
     WRITE(*,'(a)')' from Gridgen Meshfile not defined in Parameterfile'
     STOP
  END  IF
! SidePtr does not exist and is set to NULL
  SideList(iSidePtr)%Side => aSide
! next side ...
END DO
DEALLOCATE(BCEdge)
!-----------------------------------------------------------------------------------------------------------------------------------
! Sort SideList Array in order to retrieve the connectivity information efficiently
CALL SortSideList(SideList)
! Initialize Sidelists in MESH
FirstSide => NULL()
FirstBCSide => NULL()
! Go through SideList Array and create connectivity
DO iSide = 1, 2 * nSides, 2
! extract side pair
  aSide => SideList(iSide)%Side
  bSide => SideList(iSide+1)%Side
! grid file read check
  IF (SUM(SideList(iSide)%Node(:)- SideList(iSide+1)%Node(:)) .NE. 0 ) THEN
    WRITE(*,*)'Problem in connect'
    STOP
  END IF
! Save sides into the global element array
  IF (ASSOCIATED(aSide%BC).OR.ASSOCIATED(bSide%BC)) THEN
  ! Boundary Side
    ! make sure aSide is the physical side and bSide the ghost side
    ! save only aSide into the list
    IF (ASSOCIATED(aSide%BC)) THEN
      cSide => aSide
      aSide => bSide
      bSide => cSide
    END IF
    ! Save boundary side (bSide) into Boundary Side List
    ALLOCATE(aBCSide)
    aBCSide%Side => bSide
    aBCSide%nextSide => FirstBCSide
    FirstBCSide => aBCSide
  END IF
  aSide%nextSide => FirstSide
  FirstSide => aSide
! Generate Connectivity
  aSide%connection => bSide
  bSide%connection => aSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Extended Element info: Area and projection of cell onto axes
!-----------------------------------------------------------------------------------------------------------------------------------
TotalArea_q = 0.
aElem => firstElem
DO WHILE (ASSOCIATED(aElem))
  CALL CreateElementInfo(aElem)
! next Elem
  aElem => aElem%NextElem
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Generate virtual barycenters of ghostcells
!-----------------------------------------------------------------------------------------------------------------------------------
aBCSide => FirstBCSide
DO WHILE (ASSOCIATED(aBCSide))
  aSide => aBCSide%Side
  temp = 2.*ABS(DOT_PRODUCT(aSide%connection%GP, aSide%connection%n))
  aSide%Elem%bary = aSide%connection%Elem%Bary + temp * aSide%connection%n
  aBCSide => aBCSide%nextSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Side info
!-----------------------------------------------------------------------------------------------------------------------------------
! Vectors from barycenter to barycenter, Gaussian integration points and normal vectors
aSide => FirstSide
DO WHILE (ASSOCIATED(aSide))
  CALL CreateSideInfo(aSide)
! Next Side
  aSide => aSide%NextSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Periodic BCs
!-----------------------------------------------------------------------------------------------------------------------------------
CALL ConnectPeriodicBC()
!-----------------------------------------------------------------------------------------------------------------------------------
! variables for reconstruction
!-----------------------------------------------------------------------------------------------------------------------------------
aElem => firstElem
DO WHILE(ASSOCIATED(aElem))
  CALL CreateReconstructionInfo(aElem)
! Next Elem
  aElem => aElem%nextElem
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Element and side lists
!-----------------------------------------------------------------------------------------------------------------------------------
! Allocate element and side lists
ALLOCATE(Sides(1:nSides))
ALLOCATE(Elems(1:nElems))
ALLOCATE(BCSides(1:nBCSides))
! Elems
aElem => firstElem
iElem = 1
DO WHILE(ASSOCIATED(aElem))
  Elems(iElem)%Elem => aElem
  aElem => aElem%nextElem
  iElem = iElem + 1
END DO
! Sides
aSide => firstSide
iSide = 1
DO WHILE(ASSOCIATED(aSide))
  Sides(iSide)%Side => aSide
  !aSide%ID=iSide
  aSide => aSide%nextSide
  iSide = iSide + 1
END DO
! Boundaries
aBCSide => FirstBCSide
iSide = 0
DO WHILE(ASSOCIATED(aBCSide))
  IF(aBCSide%Side%BC%BCType .NE. PERIODIC) THEN
    iSide = iSide + 1
    BCSides(iSide)%Side => aBCSide%Side
    !BCSides(iSide)%Side%ID=iSide
  END IF
  aBCSide => aBCSide%nextSide
END DO
nBCSides = iSide
!-----------------------------------------------------------------------------------------------------------------------------------
DEALLOCATE(VertexPtr)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CreateMesh


SUBROUTINE CreateElementInfo( aElem)
!===================================================================================================================================
! compute the cell specific values: barycenter, area, projection length of cell
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem
USE MOD_Mesh_Vars,ONLY:TotalArea_q 
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem), POINTER     :: aElem              ! current element            !
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                  :: i
REAL                     :: vec(2), n(4,2), length(4)
!===================================================================================================================================

! Compute Barycenter of Element
aElem%Bary(:) = 0.
DO i = 1, aElem%ElemType
  aElem%Bary(:) = aElem%Bary(:) + aElem%nodearray(i)%Node%x(:)
END DO
aElem%Bary(:) = aElem%Bary(:) / aElem%ElemType
! Element area and its inverse
! side lengths and normals are computed for the element's main sides since 
! a main side can be made up of a number of subsides (hanging nodes!)
DO i = 1, aElem%ElemType-1
  vec(:)    = aElem%nodearray(i+1)%Node%x(:) - aElem%nodearray(i)%Node%x(:)
  length(i) = SQRT(vec(1)*vec(1) + vec(2) * vec(2))
  n(i,1)    = aElem%nodearray(i+1)%Node%x(2) - aElem%nodearray(i)%Node%x(2)
  n(i,2)    = aElem%nodearray(i)%Node%x(1) - aElem%nodearray(i+1)%Node%x(1)
  n(i,:)    = n(i,:) / SQRT(n(i,1)*n(i,1) + n(i,2)*n(i,2))
END DO
vec(:)                 = aElem%nodearray(1)%Node%x(:) - aElem%nodearray(aElem%ElemType)%Node%x(:)
length(aElem%ElemType) = SQRT(vec(1)*vec(1) + vec(2) * vec(2))
n(aElem%ElemType,1)    = aElem%nodearray(1)%Node%x(2) - aElem%nodearray(aElem%ElemType)%Node%x(2)
n(aElem%ElemType,2)    = aElem%nodearray(aElem%ElemType)%Node%x(1) - aElem%nodearray(1)%Node%x(1)
n(aElem%ElemType,:)    = n(aElem%ElemType,:) / SQRT(n(i,1)*n(i,1) + n(i,2)*n(i,2))
aElem%Area = length(1) * length(2) * ABS(n(1,1)*n(2,2) - n(2,1)*n(1,2))
IF(aElem%ElemType == 3) THEN
  aElem%Area = 0.5 * aElem%Area
END IF
aElem%Area_q = 1. / aElem%Area
IF (TotalArea_q /= 0.) THEN
  TotalArea_q = 1./(1./TotalArea_q + aElem%Area)
ELSE
  TotalArea_q = aElem%Area_q
END IF
! projection of cell onto x- and y-axis (needed for timestep computation)
aElem%Sx = 0.5 * DOT_PRODUCT(ABS(n(1:aElem%ElemType,1)),length(1:aElem%ElemType))
aElem%Sy = 0.5 * DOT_PRODUCT(ABS(n(1:aElem%ElemType,2)),length(1:aElem%ElemType))
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CreateElementInfo


SUBROUTINE CreateSideInfo(aSide)
!===================================================================================================================================
! create side info, e.g. the normal vector, side length, ghost cells,...
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide), POINTER     :: aSide              ! current element            
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                     :: GP(2), temp,dotProduct
!===================================================================================================================================

! Compute normal vector
aSide%n(1) = aSide%Node(2)%Node%x(2) - aSide%Node(1)%Node%x(2)
aSide%n(2) = aSide%Node(1)%Node%x(1) - aSide%Node(2)%Node%x(1)
! Compute Side length
aSide%length = SQRT(aSide%n(1)*aSide%n(1) + aSide%n(2)*aSide%n(2))
! normalize normal vector
aSide%n = aSide%n / aSide%length
! Vectors from barycenter to side midpoint
GP(:) = 0.5 * (aSide%Node(1)%Node%x(:) + aSide%Node(2)%Node%x(:))
aSide%GP(:) = GP(:) - aSide%Elem%Bary(:)

! Compute dot-product of vector from barycenter to side midpoint with normal
! vector

dotProduct = aSide%n(1)*aSide%GP(1)+aSide%n(2)*aSide%GP(2)

! If the dot-product is lower zero, the normal vector points in the wrong
! direction and must be switcht

IF (dotProduct<0.0) THEN
  aSide%n(:)=aSide%n(:)*(-1.0)
END IF

! Compute barycenter of ghost cell in case this is a boundary side
IF (ASSOCIATED(aSide%connection%BC)) THEN
  temp = 2.*ABS(DOT_PRODUCT(aSide%GP(:), aSide%n(:)))
  aSide%connection%Elem%bary(:) = aSide%Elem%bary(:) + temp * aSide%n(:)
END IF
! Connection Side
aSide%connection%GP(:) = GP(:) - aSide%connection%Elem%Bary(:)
aSide%connection%n = - aSide%n
aSide%connection%length = aSide%length
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CreateSideInfo


SUBROUTINE CreateReconstructionInfo( aElem)
!===================================================================================================================================
! compute required vectors for reconstruction
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem), POINTER     :: aElem              ! current element
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                     :: r11, r12, r22, alpha1, alpha2
TYPE(tSide), POINTER     :: aSide
INTEGER                  :: iGP
!===================================================================================================================================

r11 = 0.
r12 = 0.
r22 = 0.
aSide => aElem%firstside
DO WHILE (ASSOCIATED(aSide))
! Vectors connecting cell barycenters
  aSide%BaryBaryVec = aSide%GP - aSide%connection%GP
  aSide%BaryBaryDist = SQRT(SUM(aSide%BaryBaryVec(:) * aSide%BaryBaryVec(:)))
! Reconstruction info
  r11 = r11 + aSide%BaryBaryVec(1)*aSide%BaryBaryVec(1)
  r12 = r12 + aSide%BaryBaryVec(1)*aSide%BaryBaryVec(2)
  r22 = r22 + aSide%BaryBaryVec(2)*aSide%BaryBaryVec(2)
  aSide => aSide%nextElemSide
END DO
r11 = SQRT(r11)
r12 = r12 / r11
r22 = SQRT(r22-r12*r12)
aSide => aElem%firstside
DO WHILE (ASSOCIATED(aSide))
!-----------------------------------------------------------------------------------------------------------------------------------
! since alpha1 and 2 contain directional information (delta x), we need to determine whether
! the element is to the left or right of the side (switch in sign!)
  alpha1 = aSide%BaryBaryVec(1) / (r11*r11)
  alpha2 = (aSide%BaryBaryVec(2) - r12/r11*aSide%BaryBaryVec(1)) / (r22*r22)
  aSide%w(1) = alpha1-r12/r11*alpha2
  aSide%w(2) = alpha2
  aSide => aSide%nextElemSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Gaussian integration points and weights for volume integrals
!  (the weights include the metric)
SELECT CASE (aElem%ElemType)
CASE(3)
  aElem%nGP = 3
  ALLOCATE(aElem%wGP(aElem%nGP), aElem%xGP(aElem%nGP, NDIM))
  aSide => aElem%FirstSide
  DO iGP = 1, aElem%nGP
    aElem%xGP(iGP,:) = aElem%Bary + aSide%GP
    aElem%wGP(iGP)   = aElem%Area / 3.
    aSide => aSide%nextElemSide
  END DO
CASE(4)
  aElem%nGP = 5
  ALLOCATE(aElem%wGP(aElem%nGP), aElem%xGP(aElem%nGP, NDIM))
  aSide => aElem%FirstSide
  DO iGP = 1, aElem%nGP - 1
    aElem%xGP(iGP,:) = aElem%Bary + aSide%GP
    aElem%wGP(iGP)   = aElem%Area / 6.
    aSide => aSide%nextElemSide
  END DO
   aElem%xGP(5,:) = 0.5 * (aElem%nodearray(1)%node%x(:) + aElem%nodearray(3)%node%x(:))
   aElem%wGP(5)   = aElem%Area / 3.
CASE DEFAULT
  WRITE(*,*) '  ERROR in CreateReconstructionInfo() in CreateMesh.f90:'
  WRITE(*,*) '  Currently there is no quadrature rule for elements other than triangles or quadrangles'
  STOP
END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CreateReconstructionInfo



SUBROUTINE ConnectPeriodicBC()
!===================================================================================================================================
! create connection for periodic BCs
!===================================================================================================================================
! MODULES
USE MOD_Globals 
USE MOD_Mesh_Vars    ,ONLY:tPureSidePtr,tSide
USE MOD_Mesh_Vars    ,ONLY:FirstBCSide,FirstSide,nSides
USE MOD_Boundary_Vars,ONLY:isPeriodic
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
! Local Variable Declaration
TYPE(tPureSidePtr), POINTER :: aBCSide, sBCSide
TYPE(tSide), POINTER        :: aSide, sSide, urSide, targetSide
REAL                        :: aGPPos(2), sGPPos(2)
LOGICAL                     :: connected
INTEGER                     :: nPeriodicBC
!===================================================================================================================================

isPeriodic = .FALSE.
aBCSide => FirstBCSide
DO WHILE(ASSOCIATED(aBCSide))
  aSide => aBCSide%Side
  IF ((aSide%BC%BCType == PERIODIC).AND.(MOD(aSide%BC%BCID,10) == 1)) THEN
    isPeriodic = .TRUE.
    connected = .FALSE.
    aGPPos = aSide%GP + aSide%Elem%Bary
    nPeriodicBC = aSide%BC%BCID / 10
    sBCSide => FirstBCSide
    DO WHILE (ASSOCIATED(sBCSide))
    sSide => sBCSide%Side
    ! Check if we have a potential connection
      IF ((sSide%BC%BCType == PERIODIC).AND.(MOD(sSide%BC%BCID,10) == 2).AND.(sSide%BC%BCID / 10 == nPeriodicBC)) THEN
      ! Check if the connection works out
        sGPPos = sSide%GP + sSide%Elem%Bary
        IF (SUM(ABS(aGPPos + aSide%BC%connection - sGPPos)) <= REALTOL) THEN
        ! We have a match
          connected = .TRUE.
        ! Reconnect sides
          urSide     => aSide%Connection
          targetSide => sSide%Connection
          urSide%connection => targetSide
          targetSide%connection => urSide
        ! Reorganise lists: targetSide has to be removed from the list
          nSides = nSides - 1
          IF (ASSOCIATED(FirstSide, targetSide)) THEN
            FirstSide => FirstSide%nextSide
          ELSE
            sSide => FirstSide
            DO WHILE (ASSOCIATED(sSide))
              IF (ASSOCIATED(sSide%nextSide, targetSide)) THEN
                sSide%nextSide => targetSide%nextSide
                EXIT
              END IF
              sSide => sSide%nextSide
            END DO
          END IF
        ! Exit search loop
          EXIT
        END IF
      END IF
      sBCSide => sBCSide%nextSide
    END DO
    IF (.NOT.connected) THEN
      WRITE(*,*) 'ERROR in ConnectPeriodicBC: No Connection found.'
      WRITE(*,'(a,2(F10.5))') 'Side GP was: ', aGPPos
      STOP
    END IF
  END IF
  aBCSide => aBCSide%nextSide
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ConnectPeriodicBC


SUBROUTINE CreateCartMesh(Vertex, nVertices, BCEdge, nBCEdges, Quad, nQuads)
!===================================================================================================================================
! create cartesian mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals 
USE MOD_Mesh_Vars,ONLY: CartMesh
USE MOD_Mesh_Vars,ONLY: xMax,xMin,yMax,yMin
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER          :: nVertices, nBCEdges, nQuads
INTEGER, POINTER :: BCEdge(:,:), Quad(:,:)
REAL, POINTER    :: Vertex(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL             :: dx, dy, x, y
INTEGER          :: i, j, k, iBC
INTEGER          :: iMax, jMax
!===================================================================================================================================

! Set up necessary variables
iMax = CartMesh%iMax
jMax = CartMesh%jMax

nVertices = (iMax+1) * (jMax+1)
nQuads    = iMax * jMax
nBCEdges  = 2 * (iMax + jMax)

dx = (xMax - xMin) / iMax
dy = (yMax - yMin) / jMax

ALLOCATE(Quad(nQuads,5), BCEdge(nBCEdges,3), Vertex(nVertices, 2))
!-----------------------------------------------------------------------------------------------------------------------------------

! Build Vertices
y = yMin
k = 1
DO j = 1, jMax + 1
  x = xMin
  DO i = 1, iMax + 1 
    Vertex(k, 1) = x
    Vertex(k, 2) = y
    k = k + 1
    x = x + dx
  END DO
  y = y + dy
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Build Elements
k = 1
DO j = 1, jMax
  DO i = 1, iMax
    Quad(k,1) = k + j - 1
    Quad(k,2) = Quad(k,1) + 1
    Quad(k,3) = Quad(k,2) + iMax + 1
    Quad(k,4) = Quad(k,3) - 1
    Quad(k,5) = 1
    k = k + 1
  END DO
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Build BCEdges
k = 1
! Bottom Side
WRITE(*,*) '|Bottom Side:'
WRITE(*,*) '|  # of BCs: ', CartMesh%nBC(BOTTOM)
DO iBC = 1, CartMesh%nBC(BOTTOM)
  WRITE(*,*) '|  BC Type: ', CartMesh%BCType(BOTTOM,iBC)
  WRITE(*,'(a, i5, a, i5)') ' |  BC Range: ', CartMesh%BCRange(BOTTOM,iBC,1), ' - ', CartMesh%BCRange(BOTTOM,iBC,2)
  DO i = CartMesh%BCRange(BOTTOM,iBC,1), CartMesh%BCRange(BOTTOM,iBC,2)
    BCEdge(k,1) = i
    BCEdge(k,2) = i + 1
    BCEdge(k,3) = CartMesh%BCType(BOTTOM,iBC)
    k = k + 1
  END DO
END DO
! Right Side
WRITE(*,*) '|Right Side:'
WRITE(*,*) '|  # of BCs: ', CartMesh%nBC(RIGHT)
DO iBC = 1, CartMesh%nBC(RIGHT)
  WRITE(*,*) '|  BC Type: ', CartMesh%BCType(RIGHT,iBC)
  WRITE(*,'(a, i5, a, i5)') ' |  BC Range: ', CartMesh%BCRange(RIGHT,iBC,1), ' - ', CartMesh%BCRange(RIGHT,iBC,2)
  DO j = CartMesh%BCRange(RIGHT,iBC,1), CartMesh%BCRange(RIGHT,iBC,2)
    BCEdge(k,1) = (iMax + 1) * j
    BCEdge(k,2) = (iMax + 1) * (j + 1)
    BCEdge(k,3) = CartMesh%BCType(RIGHT,iBC)
    k = k + 1
  END DO
END DO
! Top Side
WRITE(*,*) '|Top Side:'
WRITE(*,*) '|  # of BCs: ', CartMesh%nBC(TOP)
DO iBC = 1, CartMesh%nBC(TOP)
  WRITE(*,*) '|  BC Type: ', CartMesh%BCType(TOP,iBC)
  WRITE(*,'(a, i5, a, i5)') ' |  BC Range: ', CartMesh%BCRange(TOP,iBC,1), ' - ', CartMesh%BCRange(TOP,iBC,2)
  DO i = CartMesh%BCRange(TOP,iBC,1), CartMesh%BCRange(TOP,iBC,2)
    BCEdge(k,1) = (iMax+1)*jMax + i
    BCEdge(k,2) = BCEdge(k,1) + 1
    BCEdge(k,3) = CartMesh%BCType(TOP,iBC)
    k = k + 1
  END DO
END DO
! Left Side
WRITE(*,*) '|Left Side:'
WRITE(*,*) '|  # of BCs: ', CartMesh%nBC(LEFT)
DO iBC = 1, CartMesh%nBC(LEFT)
  WRITE(*,*) '|  BC Type: ', CartMesh%BCType(LEFT,iBC)
  WRITE(*,'(a, i5, a, i5)') ' |  BC Range: ', CartMesh%BCRange(LEFT,iBC,1), ' - ', CartMesh%BCRange(LEFT,iBC,2)
  DO j = CartMesh%BCRange(LEFT,iBC,1), CartMesh%BCRange(LEFT,iBC,2)
    BCEdge(k,1) = (j-1) * (iMax + 1) + 1
    BCEdge(k,2) = BCEdge(k,1) + iMax + 1
    BCEdge(k,3) = CartMesh%BCType(LEFT,iBC)
    k = k + 1
  END DO
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CreateCartMesh


!--------------------------------------------------------------------------------------------------!
RECURSIVE SUBROUTINE SortSideList(SideList)
!===================================================================================================================================
! SortSideList:  Sorts the side list with QuickSort
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSideList), DIMENSION(:) :: SideList
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                       :: iq
!===================================================================================================================================

IF(size(SideList) > 1) THEN
  CALL PartitionSideList(SideList, iq)
  CALL SortSideList (SideList(:iq-1))
  CALL SortSideList (SideList(iq:))
ENDIF
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE SortSideList

SUBROUTINE PartitionSideList(A, marker)
!===================================================================================================================================
! PartitionSideList:
!  Sorting routine used by SortSideList above.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSideList), DIMENSION(:) :: A
INTEGER                       :: marker
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER         :: i, j
TYPE(tSideList) :: temp
INTEGER         :: x(3)              ! pivot point
!===================================================================================================================================

x(1:2) = A(1)%Node(:) ! Node1, Node2
x(3)   = A(1)%BC      ! Boundary Condition, either 0 or BC
i = 0
j = size(A) + 1

DO
  j = j-1
  DO
    IF ((A(j)%Node(1) < x(1)).OR.(((A(j)%Node(1) == x(1)).AND.(A(j)%Node(2) <= x(2))))) EXIT
     IF (((A(j)%Node(1) < x(1)).OR.(((A(j)%Node(1) == x(1)).AND.(A(j)%Node(2) < x(2))))).OR. &
                       ((A(j)%Node(1) == x(1)).AND.(A(j)%Node(2) == x(2)).AND.(A(j)%BC <= x(3)))) EXIT
    j = j-1
  END DO
  i = i+1
  DO
    IF ((A(i)%Node(1) > x(1)).OR.(((A(i)%Node(1) == x(1)).AND.(A(i)%Node(2) >= x(2))))) EXIT
     IF (((A(i)%Node(1) > x(1)).OR.(((A(i)%Node(1) == x(1)).AND.(A(i)%Node(2) > x(2))))).OR. &
                       ((A(i)%Node(1) == x(1)).AND.(A(i)%Node(2) == x(2)).AND.(A(i)%BC >= x(3)))) EXIT
    i = i+1
  END DO
  IF (i < j) THEN
    ! exchange A(i) and A(j)
    temp = A(i)
    A(i) = A(j)
    A(j) = temp
  ELSEIF (i == j) THEN
    marker = i+1
    RETURN
  ELSE
    marker = i
    RETURN
  ENDIF
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE PartitionSideList

END MODULE MOD_Mesh
