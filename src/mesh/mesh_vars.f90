MODULE MOD_Mesh_Vars
!===================================================================================================================================
! Mesh Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Boundary_Vars,ONLY:tBoundary
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!===================================================================================================================================

  ! Container for boundary conditions
  CHARACTER(LEN=256)   :: ParameterFile           ! Name of Parameter File
  CHARACTER(LEN=256)   :: strMeshFormat           ! 1= *.mesh format
  CHARACTER(LEN=256)   :: strMeshFile             ! Name of Mesh File
  CHARACTER(LEN=256)   :: strIniCondFile          ! CGNS file containing the initial condition
!-----------------------------------------------------------------------------------------------------------------------------------
! This array is part of tElem-connectivity section. Each array entry contains a pointer to a node that
! is associated with the respective element
TYPE tNodePtr
  TYPE(tNode), POINTER :: Node                         ! Pointer to a Node in the NodeList
END TYPE tNodePtr
!-----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------!
TYPE tElemPtr
  TYPE(tElem), POINTER :: Elem                         ! Pointer to an Element in the Element list
END TYPE tElemPtr
!-----------------------------------------------------------------------------------------------------------------------------------
! This type constitutes the Node list. The list head can be addressed by =>MESH%FirstNode
TYPE tNode
  INTEGER   :: ID                                     ! Unique Node ID (from PREPMESH)
  TYPE(tNode), POINTER :: nextNode                    ! Pointer to the next Node in the NodeList
  REAL                 :: x(1:2)                      ! x- and y-coordinates of the Node in real space
END TYPE tNode
!-----------------------------------------------------------------------------------------------------------------------------------
! A pointer structure used for generating an array-bases list of sides
TYPE tSidePtr
  TYPE(tSide), POINTER        :: Side             ! Pointer to the side
END TYPE tSidePtr
!-----------------------------------------------------------------------------------------------------------------------------------
! A pure side pointer used for generating lists of sides (i.e. for the profile)
TYPE tPureSidePtr
  TYPE(tPureSidePtr), POINTER :: nextSide             ! Pointer to the next entry in the the sidepointer list
  TYPE(tSide), POINTER        :: Side                 ! Pointer to the actual Side in the Side List
END TYPE tPureSidePtr
!-----------------------------------------------------------------------------------------------------------------------------------
! This type constitutes the Side list. The list head can be addressed by =>MESH%FirstSide
TYPE tSide
  INTEGER                  :: ID              ! Unique Side ID (from PREPMESH)
  INTEGER                  :: BCType          ! Boundary Condition Type
  INTEGER                  :: BCID            ! Sub-ID of Boundary Condition
  TYPE(tBoundary), POINTER :: BC              ! Pointer to Boundary Condition
  REAL                     :: pvar(NVAR)      ! Side State
  REAL                     :: n(2)            ! normal vector
  REAL                     :: length          ! length of the edge
  REAL                     :: BaryBaryVec(2)  ! Vector (x,y) from Element Barycenter to Neighbor Barycenter
  REAL                     :: BaryBaryDist    ! Length of BaryBaryVec
  REAL                     :: GP(1:2)         ! x- and y-component of the vector from the cell barycenter to
                                              ! the gaussian point at the side midpoint
  REAL                     :: w(1:2)          ! omega1 and 2 entries for 2nd order gradient reconstruction
                                              ! see Blazek Chapter 5.3.4
  REAL                     :: flux(NVAR)      ! flux over side (with respect to element)
  TYPE(tSide), POINTER     :: connection      ! neighbor side
  TYPE(tSide), POINTER     :: nextElemSide    ! Pointer to the next Side in the local SideList
  TYPE(tSide), POINTER     :: nextSide        ! Pointer to the next Side in the global SideList
  TYPE(tNodePtr)           :: Node(2)         ! Pointer to the 2 nodes that belong to the Side
  TYPE(tElem), POINTER     :: Elem            ! Pointer to the left and right element adjacent to the side
END TYPE tSide
!-----------------------------------------------------------------------------------------------------------------------------------
! This type constitutes the Element list.
TYPE tElem
  INTEGER   :: ElemType             ! triangle (3) or quadrangle (4)
  INTEGER   :: ID                   ! Unique Element ID
  INTEGER   :: domain               ! flow domain (part of a mesh)
  REAL      :: Bary(2)              ! Barycenter of cell (xy)  (2nd order)
  REAL      :: Sx                   ! cell extension in x-dir for time step
  REAL      :: Sy                   ! cell extension in y-dir for time step
  REAL      :: Area                 ! area of the cell
  REAL      :: Area_q               ! 1/area of the cell
  REAL      :: pvar(NVAR)           ! Primitive Variables
  REAL      :: cvar(NVAR)           ! Conservative Variables
  REAL      :: cvar_stage(NVAR)     ! Conservative Variables at initial RK stage (needed for 2nd order RK time stepping)
  REAL      :: u_x(NVAR)            ! x-gradient of pvar
  REAL      :: u_y(NVAR)            ! y-gradient of pvar
  REAL      :: u_t(NVAR)            ! t-gradient of pvar
  REAL      :: source(NVAR)         ! source term
  REAL      :: dt                   ! Cell timestep
  REAL      :: dt_loc               ! Cell timestep
  REAL      :: venk_epsilon_sq      ! constant for venkatakrishnan's limiter
  INTEGER   :: innerSides           ! Number of non-BC Sides

  ! REAL      :: K(NVAR)              ! Runge-Kutta stage update

  INTEGER       :: nGP              ! Number of gaussian integration points
  REAL, POINTER :: xGP(:,:)         ! Gaussian points for volume integrals
  REAL, POINTER :: wGP(:)           ! Gaussian weights. The metric is already included.

  TYPE(tSide), POINTER        :: firstSide         ! Pointer to the first side of the element
  TYPE(tElem), POINTER        :: nextElem          ! Pointer to the next Element in the Element List
  TYPE(tNodePtr), POINTER     :: nodearray(:)      ! Array of pointers to the nodes of the Element
                                                   ! (size 3 or 4 for visualization output
END TYPE tElem
!-----------------------------------------------------------------------------------------------------------------------------------
! This Type contains all information for the cartesian mesh generator
TYPE tCartMesh
  INTEGER      :: imax, jmax
  INTEGER      :: nBC(4)
  INTEGER      :: BCType(4, 20)      ! Number of different BC per side are currently limited to 20
  INTEGER      :: BCRange(4, 20, 2)
END TYPE tCartMesh
!-----------------------------------------------------------------------------------------------------------------------------------
! The type constitutes the pointer-based mesh structure. It contains the global mesh data and the list heads for the
! element, side and nodelists
TYPE(tCartMesh)        :: CartMesh
CHARACTER(LEN=256)     :: GridFile                 ! Name of CGNS gridfile

INTEGER                :: MeshType                 ! 0= unstruct 1= cartesian mesh!
INTEGER                :: MeshFormat               ! 0= EMC2 .mesh, 1=CGNS

INTEGER                :: nNodes                   ! number of nodes

INTEGER                :: nElems                   ! number of cells
INTEGER                :: nTrias                   ! number of triangles
INTEGER                :: nQuads                   ! number of quadrangles

INTEGER                :: nSides                   ! number of Sides
INTEGER                :: nBCSides                 ! number of Boundary Sides
INTEGER                :: nInnerSides              ! number of Inner Sides

REAL                   :: totalArea_q              ! inverse of total area of all cells
REAL                   :: xMin, xMax
REAL                   :: yMin, yMax
REAL                   :: dxRef

! Lists
TYPE(tElemPtr), POINTER :: Elems(:)                ! Element List
TYPE(tSidePtr), POINTER :: Sides(:)                ! Side List
TYPE(tSidePtr), POINTER :: BCSides(:)              ! Ghost Side List


TYPE(tElem), POINTER   :: FirstElem                ! Pointer to the head of the element
                                                   ! List
TYPE(tNode), POINTER   :: FirstNode                ! Pointer to the head of the node list

TYPE(tSide), POINTER   :: FirstSide                ! head of the inner side list
TYPE(tPureSidePtr), POINTER :: FirstBCSide         ! head of the boundary side list
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE MOD_Mesh_Vars
