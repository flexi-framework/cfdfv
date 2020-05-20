MODULE MOD_Mesh_readin
!===================================================================================================================================
! Mesh Readin
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ReadEMC2
   MODULE PROCEDURE ReadEMC2
END INTERFACE
INTERFACE ReadGmsh
   MODULE PROCEDURE ReadGmsh
END INTERFACE
INTERFACE ReadCGNS
   MODULE PROCEDURE ReadCGNS
END INTERFACE
INTERFACE ReadGridgenCGNS
   MODULE PROCEDURE ReadGridgenCGNS
END INTERFACE
INTERFACE errorCGNS
   MODULE PROCEDURE errorCGNS
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: ReadCGNS, ReadEMC2,ReadGmsh, ReadGridgenCGNS,errorCGNS
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadEMC2(FileName, Vertex, nVertices, BCEdge, nBCEdges, Tria, nTrias, Quad, nQuads        )
!===================================================================================================================================
! Read EMC2 mesh
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:
USE MOD_Readin, ONLY: getFreeIOUnit,getCmdLine
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(256)   :: FileName
INTEGER          :: nVertices, nBCEdges, nTrias, nQuads
INTEGER, POINTER :: BCEdge(:,:), Tria(:,:), Quad(:,:), Temp(:,:)
REAL, POINTER    :: Vertex(:,:)
! Local variable declaration
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: i, MeshUnit, stat
INTEGER          :: nEdges, localEdge(3)
CHARACTER(256)   :: actLine
!===================================================================================================================================

!read from emc2 file in .mesh format
CALL getFreeIOUnit(100,MeshUnit )                                         ! Get a free IO-Unit greater then 100
OPEN(UNIT = MeshUnit, FILE = TRIM(FileName), &
     STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT=stat)
IF (stat .NE. 0) THEN                                                      ! Error Handler
   WRITE(*,*)'ERROR: cannot open meshfile, stopped in CreateMesh.f90!'
   STOP                                                                    ! Error Handler
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
!read in VERTICES
REWIND(UNIT=MeshUnit)
actLine = ''
DO WHILE((TRIM(actLine).ne.'Vertices').and.(TRIM(actLine).ne.'End'))
   CALL getCmdLine(MeshUnit, actLine)
END DO
IF (TRIM(actLine).eq.'End') THEN
   PRINT*, 'No Vertices found in meshfile , stop in CreatePrepMesh'
   STOP
END IF
! Read number of vertices
CALL getCmdLine(MeshUnit, actLine) ! nVertices
read(actLine, "(i6)") nVertices
! allocate memory for vertices
ALLOCATE(Vertex(1:nVertices, 2))
DO i = 1, nVertices
  CALL getCmdLine(MeshUnit, actLine)
  read(actLine, *) Vertex(i, 1:2)
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! read in BCEDGES (only Edges of boundary and domain boundaries)
REWIND(UNIT=MeshUnit)
actLine = ''
DO WHILE((TRIM(actLine).ne.'Edges').and.(TRIM(actLine).ne.'End'))
   CALL getCmdLine(MeshUnit, actLine)
END DO
IF (TRIM(actLine).eq.'End') THEN
   PRINT*, 'No Boundary Edges found in meshfile, stop in CreatePrepMesh'
   STOP
END IF
! Read number of edges
CALL getCmdLine(MeshUnit, actLine)
read(actLine, "(i6)") nEdges
! allocate memory for edge I/O
ALLOCATE(Temp(1:nEdges,1:3))
! Read edge information
! All Edges with ID 0 (emc2 garbage) are sorted out
nBCEdges = 0
DO i = 1, nEdges
  CALL getCmdLine(MeshUnit, actLine)
  read(actLine, *) localEdge(1:3)
  IF (localEdge(3) /= 0) THEN
    nBCEdges = nBCEdges + 1
    Temp(nBCEdges,:) = localEdge(:)
  END IF
END DO
! Allocate memory for edges
ALLOCATE(BCEdge(nBCEdges, 3))
! Save edge data
BCEdge(:,:) = Temp(1:nBCEdges,:)
! Deallocate temporary array
DEALLOCATE(Temp)
! Output
WRITE(*,'(a, I5, a)') '     ', nBCEdges, ' Boundary Edges read from .mesh File.'
!-----------------------------------------------------------------------------------------------------------------------------------
! read in Trias
REWIND(UNIT = MeshUnit)
actLine = ''
DO WHILE((TRIM(actLine).ne.'Triangles').and.(TRIM(actLine).ne.'End'))
   CALL getCmdLine(MeshUnit, actLine)
END DO
IF (TRIM(actLine).eq.'End') THEN
  nTrias = 0
ELSE
! read number of triangles
  CALL getCmdLine(MeshUnit, actLine)
  read(actLine, "(i6)") nTrias
! allocate memory for triangles
  ALLOCATE(Tria(nTrias, 4))
! read triangles
  DO i = 1, nTrias
    CALL getCmdLine(MeshUnit, actLine)
    read(actLine, *) Tria(i,1:4)
  END DO
END IF
! Output
WRITE(*,'(a, I8, a)') '     ', nTrias, ' Triangles read from .mesh File.'
!-----------------------------------------------------------------------------------------------------------------------------------
! read in Quads
REWIND(UNIT = MeshUnit)
actLine = ''
DO WHILE((TRIM(actLine).ne.'Quadrangles').and.(TRIM(actLine).ne.'End'))
   CALL getCmdLine(MeshUnit, actLine)
END DO
IF (TRIM(actLine).eq.'End') THEN
  nQuads = 0
ELSE
! read number of quadrangles
  CALL getCmdLine(MeshUnit, actLine)
  read(actLine, "(i6)") nQuads
! allocate memory for quadrangles
  ALLOCATE(Quad(nQuads, 5))
! read quadrangles
  DO i=1,nQuads
     CALL getCmdLine(MeshUnit, actLine)
     read(actLine, *) Quad(i,1:5)   !*** listengesteuert einlesen
  END DO
END IF
! Output
WRITE(*,'(a, I5, a)') '     ', nQuads, ' Quadrangles read from .mesh File.'
CLOSE(MeshUnit)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ReadEMC2

SUBROUTINE ReadGmsh(FileName, Vertex, nVertices, BCEdge, nBCEdges, Tria, nTrias, Quad, nQuads )
!===================================================================================================================================
! Read in GMSH meshes
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:
USE MOD_Readin, ONLY: getFreeIOUnit
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(256)   :: FileName
INTEGER          :: nVertices, nBCEdges, nTrias, nQuads
INTEGER, POINTER :: BCEdge(:,:), Tria(:,:), Quad(:,:)
INTEGER, POINTER :: BCEdgeTemp(:,:), TriaTemp(:,:), QuadTemp(:,:)
REAL, POINTER    :: Vertex(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: MeshUnit, stat, i, j, discard, nElem, elemType, nTag, &
                    pGroup, id
CHARACTER(256)   :: actLine
!===================================================================================================================================

!read from emc2 file in .mesh format
CALL getFreeIOUnit(100,MeshUnit)                                              ! Get a free IO-Unit greater then 100
OPEN(UNIT = MeshUnit, FILE = TRIM(FileName), &
     STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT=stat)
IF (stat .NE. 0) THEN                                                      ! Error Handler
   WRITE(*,*)'ERROR: cannot open meshfile, stopped in CreateMesh.f90!'
   STOP                                                                    ! Error Handler
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
!read in Nodes
REWIND(UNIT=MeshUnit)
actLine = ''
DO WHILE ( TRIM(ADJUSTL(actLine)) .NE. '$Nodes' .AND. stat .EQ. 0 )
  READ(MeshUnit, '(a)', IOSTAT = stat) actLine
END DO
IF ( stat .NE. 0 ) THEN
  WRITE(*,*) 'ERROR: no nodes in Gmsh file'
  STOP
END IF

! number of nodes
READ(MeshUnit,*) nVertices
! TODO error handler

ALLOCATE(Vertex(nVertices,2))

! nodes
! NOTE: This is an oversimplification. The Gmsh file format does not
! require the Nodes to appear in a particular order as each of them
! is assigned a unique id. This routine uses the fact that the Gmsh
! editor at the time of this writing conveniently does so.
DO i = 1, nVertices
  READ(MeshUnit,*, IOSTAT = stat) id, Vertex(i, 1:2) ! node-id, x, y
  IF ( stat .NE. 0 ) THEN
    WRITE(*,*) 'ERROR: I/O error while reading Gmsh file: ', stat
  END IF
  IF ( id .NE. i ) THEN
    WRITE(*,*) '    node-id ', id, ' does not match node position.'
    WRITE(*,*) '    This WILL cause trouble. Aborting.'
    STOP
  END IF
END DO

WRITE(*,'(a, I8, a)') '     ', nVertices, ' Nodes read from .msh File.'

nBCEdges = 0
nTrias = 0
nQuads = 0

! look for elements section
REWIND(UNIT = MeshUnit)
actLine = ''
DO WHILE ( TRIM(ADJUSTL(actLine)) .NE. '$Elements' .AND. stat .EQ. 0 )
  READ(MeshUnit, '(a)', IOSTAT=stat) actLine
END DO
IF ( stat .NE. 0 ) THEN
  WRITE(*,*) 'ERROR: no elements in Gmsh file'
  STOP
END IF

! number of elements
READ(MeshUnit,*, IOSTAT = stat) nElem
IF ( stat .NE. 0 ) THEN
  WRITE(*,*) 'ERROR: I/O error while reading Gmsh file: ', stat
  STOP
END IF

! allocate space for the elements
! NOTE: More memory than is acutally needed is allocated here as the
! actual number of the different element types is not known a priori.
! This memory is freed later on though.
i = 0
ALLOCATE(BCEdgeTemp(nElem, 3), STAT = stat)
i = i + ABS(stat)
ALLOCATE(TriaTemp(nElem, 4), STAT = stat)
i = i + ABS(stat)
ALLOCATE(QuadTemp(nElem, 5), STAT = stat)
i = i + ABS(stat)
IF ( i .NE. 0 ) THEN
  WRITE(*,*) 'ERROR: Could not allocate memory for mesh elements'
  STOP
END IF

! read in elements
DO i = 1, nElem
  ! buffer the active line for consumption
  READ(MeshUnit,'(a)', IOSTAT = stat) actLine
  IF ( stat .NE. 0 ) THEN
    WRITE(*,*) 'ERROR: I/O error while reading Gmsh file: ', stat
    STOP
  END IF
  ! read the element-id
  actLine = TRIM(ADJUSTL(actLine))
  CALL ConsumeInteger(actLine, id, stat)
  IF ( stat .LE. 0 ) THEN
    WRITE(*,*) 'ERROR: The ', i, 'th element in the Gmsh file is broken'
    STOP
  END IF

  ! element-type
  actLine = TRIM(ADJUSTL(actLine))
  CALL ConsumeInteger(actLine, elemType, stat)
  IF ( stat .LE. 0 ) THEN
    WRITE(*,*) 'ERROR: The element with id ', id, ' in the Gmsh file is broken'
    STOP
  END IF

  ! number of integer tags to read
  actLine = TRIM(ADJUSTL(actLine))
  CALL ConsumeInteger(actLine, nTag, stat)
  IF ( stat .LE. 0 ) THEN
    WRITE(*,*) 'ERROR: The element with id ', id, ' in the Gmsh file is broken'
    STOP
  END IF

  ! first integer tag denotes physical group of the element
  actLine = TRIM(ADJUSTL(actLine))
  CALL ConsumeInteger(actLine, pGroup, stat)
  IF ( stat .LE. 0 ) THEN
    WRITE(*,*) 'ERROR: The element with id ', id, ' in the Gmsh file is broken'
    STOP
  END IF

  ! all further integer tags are discarded
  DO j = 2, nTag
    actLine = TRIM(ADJUSTL(actLine))
    CALL ConsumeInteger(actLine, discard, stat)
    IF ( stat .LE. 0 ) THEN
      WRITE(*,*) 'ERROR: The element with id ', id, ' in the Gmsh file is broken'
      STOP
    END IF
  END DO

  ! handle different element types
  SELECT CASE(elemType)
  CASE(1) ! 2 node line
    IF ( pGroup .GT. 100 ) THEN ! physical group above 100 denotes BC
      nBCEdges = nBCEdges + 1
      BCEdgeTemp(nBCEdges, 3) = pGroup

      READ(actLine,*, IOSTAT = stat) BCEdgeTemp(nBCEdges, 1:2) ! read two end nodes
      IF ( stat .NE. 0 ) THEN
        WRITE(*,*) 'ERROR: The element with id ', id, ' in the Gmsh file is broken'
        STOP
      END IF
    END IF
  CASE(2) ! 3 node triangle
    nTrias = nTrias + 1
    TriaTemp(nTrias, 4) = pGroup

    READ(actLine,*, IOSTAT = stat) TriaTemp(nTrias, 1:3) ! read three corners
    IF ( stat .NE. 0 ) THEN
      WRITE(*,*) 'ERROR: The element with id ', id, ' in the Gmsh file is broken'
      STOP
    END IF
  CASE(3) ! 4 node quadangle
    nQuads = nQuads + 1
    QuadTemp(nQuads, 5) = pGroup

    READ(actLine,*, IOSTAT = stat) QuadTemp(nQuads, 1:4) ! read four corners
    IF ( stat .NE. 0 ) THEN
      WRITE(*,*) 'ERROR: The element with id ', id, ' in the Gmsh file is broken'
      STOP
    END IF
  END SELECT
END DO

! allocate right amount of space and free temporary space
i = 0
ALLOCATE(BCEdge(nBCEdges, 3), STAT = stat)
i = i + ABS(stat)
BCEdge(:,:) = BCEdgeTemp(1:nBCEdges,:)
DEALLOCATE(BCEdgeTemp)
WRITE(*,'(a, I8, a)') '     ', nBCEdges, ' Boundary Edges read from .msh File.'

ALLOCATE(Tria(nTrias, 4), STAT = stat)
i = i + ABS(stat)
Tria(:,:) = TriaTemp(1:nTrias,:)
DEALLOCATE(TriaTemp)
WRITE(*,'(a, I8, a)') '     ', nTrias, ' Triangles read from .msh File.'

ALLOCATE(Quad(nQuads, 5), STAT = stat)
i = i + ABS(stat)
Quad(:,:) = QuadTemp(1:nQuads,:)
DEALLOCATE(QuadTemp)
WRITE(*,'(a, I8, a)') '     ', nQuads, ' Quadrangles read from .msh File.'

IF ( i .NE. 0 ) THEN
  WRITE(*,*) 'ERROR: Could not allocate memory for mesh elements'
  STOP
END IF

CLOSE(MeshUnit)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ReadGmsh

SUBROUTINE ConsumeInteger(buf, n, stat)
!===================================================================================================================================
! Rad integer from string
! ConsumeInteger reads an integer from the beginning of a string and
! returns this integer and the rest of the string
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/Output  VARIABLES
CHARACTER(LEN=*)  :: buf    ! string to consume from
INTEGER           :: n      ! result is stored here
INTEGER           :: stat   ! number of characters read is returned,
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=*), PARAMETER :: digits = '0123456789' ! for matching
INTEGER                     :: i
CHARACTER(LEN=16)           :: t ! space for dynamic format string
!===================================================================================================================================

! determine length of integer string representation
i = 0
DO WHILE ( i .LE. len(buf) .AND. scan(digits, buf(i+1:i+1)) .NE. 0 )
  i = i + 1
END DO

! no integer at beginning of buffer, return with 0 characters read
IF ( i .EQ. 0 ) THEN
  stat = 0
  RETURN
END IF

! construct dynamic format string
WRITE(t,"(i16)") i

! read integer from file
READ(buf, "(i" // trim(adjustl(t)) // ")") n
! drop beginning of file
buf = buf(i+1:)
stat = i
END SUBROUTINE ConsumeInteger


SUBROUTINE ReadCGNS(FileName, Vertices, nVertices, BCEdge, nBCEdges, Tria, nTrias, Quad, nQuads )
!===================================================================================================================================
! Read in CGNS mesh for restart
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:
USE MOD_Readin, ONLY: getFreeIOUnit, getCmdLine
!-----------------------------------------------------------------------------------------------------------------------------------
! Include CGNS Library:
! (Please note that the CGNS library has to be installed in the computer's
! library and include path (see CGNS documentation for more information:
! www.cgns.org)
USE CGNS
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=256):: FileName
INTEGER           :: nVertices, nBCEdges, nTrias, nQuads
INTEGER, POINTER  :: BCEdge(:,:), Tria(:,:), Quad(:,:)
REAL, POINTER     :: Vertices(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ierr, i, iTria, iQuad, iBCEdge
INTEGER           :: CGNSUnit
INTEGER           :: BaseIndex, ZoneIndex, SectionIndex
INTEGER           :: iSize(1,3), nElems, SectionType, el_start, el_end
INTEGER           :: nbndry, parent_flag, DataSize, normallist, iBC
INTEGER           :: BCCode
CHARACTER(LEN=32) :: ZoneName, SectionName,BocoName
INTEGER, POINTER  :: Elems(:)
INTEGER, POINTER  :: bcPnts(:)

INTEGER :: BCType, ptset_type, nBCPnts, NormalIndex(2), NormalListFlag
INTEGER :: NormalDataType, nDataSets
!===================================================================================================================================

! open cgns file
CALL cg_open_f(TRIM(FileName), MODE_READ, CGNSUnit, ierr)
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_open_f')
! We know there is only one base and only one zone:
! This has to be changed in the future
BaseIndex    = 1
ZoneIndex    = 1
! Get Nodes
CALL cg_zone_read_f(CGNSUnit, BaseIndex, ZoneIndex, ZoneName, iSize, ierr)
nVertices = iSize(1,1)
ALLOCATE(Vertices(nVertices, 2))
CALL cg_coord_read_f(CGNSUnit, BaseIndex, ZoneIndex, &
                     'CoordinateX', RealDouble, 1, nVertices, Vertices(:,1), ierr)
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_coord_read_f')
CALL cg_coord_read_f(CGNSUnit, BaseIndex, ZoneIndex, &
                     'CoordinateY', RealDouble, 1, nVertices, Vertices(:,2), ierr)
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_coord_read_f')
! We know there are only two sections:
! Get number of entries in Element connectivity
SectionIndex = 1
CALL cg_section_read_f(CGNSUnit,     &
                       BaseIndex,    &
                       ZoneIndex,    &
                       SectionIndex, &
                       SectionName,  &
                       SectionType,  &
                       el_start,     &
                       el_end,       &
                       nbndry,       &
                       parent_flag,  &
                       ierr          )
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_section_read_f')
! Get number of entries in Element connectivity
CALL cg_ElementDataSize_f(CGNSUnit,     &
                          BaseIndex,    &
                          ZoneIndex,    &
                          SectionIndex, &
                          DataSize,     &
                          ierr          )
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_ElementDataSize_f')
! Allocate Data arrays for Mesh information
ALLOCATE(Elems(1:DataSize))
! Read Element connectivity
CALL cg_elements_read_f(CGNSUnit,     &
                        BaseIndex,    &
                        ZoneIndex,    &
                        SectionIndex, &
                        Elems,        &
                        NULL,         &
                        ierr          )
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_elements_read_f')
! Extract triangles and quadrangles from mixed Elems array
! 1st step: count elements
nTrias = 0
nQuads = 0
i = 1
DO WHILE (i < DataSize)
  SELECT CASE (Elems(i))
  CASE(TRI_3)
    nTrias = nTrias + 1
    i = i + 4
  CASE(QUAD_4)
    nQuads = nQuads + 1
    i = i + 5
  CASE DEFAULT
    WRITE(*,*) 'Error in ReadCGNS: Only Tri_3 and Quad_4 Elements allowed.'
    STOP
  END SELECT
END DO
nElems = nTrias + nQuads
! 2nd step: Allocate arrays for element connectivity
ALLOCATE(Tria(nTrias, 1:4), Quad(nQuads, 1:5))
i = 1
iTria = 1
iQuad = 1
DO WHILE (i < DataSize)
  SELECT CASE (Elems(i))
  CASE(TRI_3)
    Tria(iTria,1:3) = Elems(i+1:i+3)
    Tria(iTria,4) = ZoneIndex
    iTria = iTria + 1
    i = i + 4
  CASE(QUAD_4)
    Quad(iQuad,1:4) = Elems(i+1:i+4)
    Quad(iQuad,5) = ZoneIndex
    iQuad = iQuad + 1
    i = i + 5
  END SELECT
END DO
! Deallocate Elems array
DEALLOCATE(Elems)
! Read boundaries
SectionIndex = 2
CALL cg_section_read_f(CGNSUnit,     &
                       BaseIndex,    &
                       ZoneIndex,    &
                       SectionIndex, &
                       SectionName,  &
                       SectionType,  &
                       el_start,     &
                       el_end,       &
                       nbndry,       &
                       parent_flag,  &
                       ierr          )
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_section_read_f')
CALL cg_ElementDataSize_f(CGNSUnit,     &
                          BaseIndex,    &
                          ZoneIndex,    &
                          SectionIndex, &
                          DataSize,     &
                          ierr          )
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_ElementDataSize_f')
! Allocate Array for boundary elements
ALLOCATE(Elems(1:DataSize))
! Read Boundary Elements
CALL cg_elements_read_f(CGNSUnit,     &
                        BaseIndex,    &
                        ZoneIndex,    &
                        SectionIndex, &
                        Elems,        &
                        NULL,         &
                        ierr          )
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_elements_read_f')
! Read number of boundaries
CALL cg_nbocos_f(CGNSUnit,     &
                 BaseIndex,    &
                 ZoneIndex,    &
                 nbndry,       &
                 ierr          )
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_nbocos_f')
! Get number of boundary sides
nBCEdges = 0
DO iBC = 1, nbndry
! Get boundary information
  CALL cg_boco_info_f(CGNSUnit,       &
                      BaseIndex,      &
                      ZoneIndex,      &
                      iBC,            &
                      BocoName,       &
                      BCType,         &
                      ptset_type,     &
                      nBCPnts,        &
                      NormalIndex,    &
                      NormalListFlag, &
                      NormalDataType, &
                      nDataSets,      &
                      ierr            )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_boco_info_f')
  nBCEdges = nBCEdges + nBCPnts
END DO
! Allocate Array for boundary information
ALLOCATE(BCEdge(nBCEdges, 3))
! Read boundary sides
iBCEdge = 1
DO iBC = 1, nbndry
! Get boundary information
  CALL cg_boco_info_f(CGNSUnit,       &
                      BaseIndex,      &
                      ZoneIndex,      &
                      iBC,            &
                      BocoName,       &
                      BCType,         &
                      ptset_type,     &
                      nBCPnts,        &
                      NormalIndex,    &
                      NormalListFlag, &
                      NormalDataType, &
                      nDataSets,      &
                      ierr            )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_boco_info_f')
! Allocate storage for current boundary
  ALLOCATE(bcpnts(nBCPnts))
  CALL cg_boco_read_f(CGNSUnit,    &
                      BaseIndex,   &
                      ZoneIndex,   &
                      iBC,         &
                      bcpnts,      &
                      normallist,  &
                      ierr         )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_boco_read_f')
  READ(BocoName,'(I4)') BCCode
  DO i = 1, nBCPnts
    BCEdge(iBCEdge, 1) = Elems((bcpnts(i)-1-nElems)*3+2)
    BCEdge(iBCEdge, 2) = Elems((bcpnts(i)-1-nElems)*3+3)
    BCEdge(iBCEdge, 3) = BCCode
    iBCEdge = iBCEdge + 1
  END DO
  DEALLOCATE(bcpnts)
END DO
CALL cg_close_f(CGNSUnit, ierr)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ReadCGNS

SUBROUTINE ReadGridgenCGNS(FileName, Vertices, nVertices, BCEdge, nBCEdges, Tria, nTrias, Quad, nQuads, ZoneConnect,nZones )
!===================================================================================================================================
! read CGNS mesh from gridgen
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:
USE MOD_Readin, ONLY: getFreeIOUnit, getCmdLine
!-----------------------------------------------------------------------------------------------------------------------------------
! Include CGNS Library:
! (Please note that the CGNS library has to be installed in the computer's
! library and include path (see CGNS documentation for more information:
! www.cgns.org)
USE CGNS
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=256):: FileName
INTEGER           :: nVertices, nBCEdges, nTrias, nQuads
INTEGER           :: vS,vE,tS,tE,qS,qE,eS,eE
INTEGER, POINTER  :: BCEdge(:,:), Tria(:,:), Quad(:,:),ZoneConnect(:)
REAL, POINTER     :: Vertices(:,:)
INTEGER,POINTER   :: nVerticesZone(:),nBCEdgesZone(:), nTriasZone(:), nQuadsZone(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: ierr, i
INTEGER           :: iConn, nConns, nConnPnts
INTEGER           :: iBC,nbndry, BCCode,vOffset,zoneNumber
INTEGER           :: BaseIndex, ZoneIndex, SectionIndex , nZones
INTEGER           :: CGNSUnit
INTEGER           :: iSize(1,3),  SectionType, el_start, el_end
INTEGER           :: parent_flag, DataSize, normallist
CHARACTER(LEN=32) :: SectionName, BocoName,ConnName,DonorName
CHARACTER(LEN=32),POINTER :: ZoneName(:)
INTEGER, POINTER  :: Elems(:)
INTEGER, POINTER  :: bcPnts(:),ConnPnts(:,:)
REAL              :: orient

INTEGER           :: BCType, ptset_type, nBCPnts, NormalIndex(2)
INTEGER           ::    NormalListFlag, NormalDataType, nDataSets
INTEGER           :: location, connect_type, donor_zonetype, donor_ptset_type
INTEGER           ::    donor_datatype, nData_Donor
!===================================================================================================================================


! open cgns file
CALL cg_open_f(TRIM(FileName), MODE_READ, CGNSUnit, ierr)  ! output needed: CGNSunit
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_open_f')
! We know there is only one base:
BaseIndex    = 1
! We know there is  only one section:
 SectionIndex = 1
! number of zones
CALL cg_nzones_f(CGNSUnit,BaseIndex,nZones,ierr)           ! output needed: nZones
  IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_nzones_f')
! allocate size variables
ALLOCATE(nVerticesZone(nZones),nBCEdgesZone(nZones),nTriasZone(nZones))
ALLOCATE(nQuadsZone(nZones),ZoneName(nZones))
! to allocate vectors, loop over zones
nTriasZone(:)=0
nQuadsZone(:)=0
nBCEdgesZone(:)=0
nVertices=0
nTrias=0
nQuads=0
nBCEdges=0
DO ZoneIndex=1,nZones
  ! count Nodes
  CALL cg_zone_read_f(CGNSUnit,             &
                      BaseIndex,            &
                      ZoneIndex,            &
                      ZoneName(Zoneindex),  &     !output needed name of Zone
                      iSize,                &     !output needed: Sizes of Zones
                      ierr)
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_zoneread_f')
  nVerticesZone(ZoneIndex) = iSize(1,1)
  ! count Elements
  CALL cg_section_read_f(CGNSUnit,     &
                         BaseIndex,    &
                         ZoneIndex,    &
                         SectionIndex, &
                         SectionName,  &   ! Name of section
                         SectionType,  &   ! Type of Element distibution, in Gridgen only TRIA_3 or QUAD_4
                         el_start,     &
                         el_end,       &
                         nbndry,       &
                         parent_flag,  &
                         ierr          )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_section_read_f')
  CALL cg_ElementDataSize_f(CGNSUnit,     &
                            BaseIndex,    &
                            ZoneIndex,    &
                            SectionIndex, &
                            DataSize,     &  ! output needed: DataSize
                            ierr          )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_ElementDataSize_f')
  SELECT CASE (SectionType)
  CASE(TRI_3)
    nTriasZone(ZoneIndex)=DataSize/3
  CASE(QUAD_4)
    nQuadsZone(ZoneIndex)=DataSize/4
  CASE DEFAULT
    WRITE(*,*) 'Error in ReadGridgenCGNS: Only Tri_3 and Quad_4 Elements allowed.'
    STOP
  END SELECT
  ! Read number of boundaries
  CALL cg_nbocos_f(CGNSUnit,     &
                   BaseIndex,    &
                   ZoneIndex,    &
                   nbndry,       &
                   ierr          )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_nbocos_f')
  ! Get number of boundary sides
  nBCEdges=0
  DO iBC = 1, nbndry
  ! Get boundary information
    CALL cg_boco_info_f(CGNSUnit,       &
                        BaseIndex,      &
                        ZoneIndex,      &
                        iBC,            &
                        BocoName,       &
                        BCType,         &
                        ptset_type,     &
                        nBCPnts,        &
                        NormalIndex,    &
                        NormalListFlag, &
                        NormalDataType, &
                        nDataSets,      &
                        ierr            )
      IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_boco_info_f')
    nBCEdges = nBCEdges + (nBCPnts-1)
  END DO
  nBCEdgesZone(ZoneIndex)= nBCEdges
END DO !Zoneindex=1,nZones
!-----------------------------------------------------------------------------------------------------------------------------------
! Allocate vertex array
nVertices=SUM(nVerticesZone(:))
ALLOCATE(Vertices(nVertices, 2),ZoneConnect(nVertices))
ZoneConnect(:)=0
! Allocate arrays for element connectivity
nTrias=SUM(nTriasZone(:))
nQuads=SUM(nQuadsZone(:))
ALLOCATE(Tria(nTrias, 1:4), Quad(nQuads, 1:5))
! Allocate array for BCedge connectivity
nBCEdges=SUM(nBCEdgesZone(:))
ALLOCATE(BCEdge(nBCEdges, 3))
WRITE(*,*)'       total number of zones  : ',nZones
WRITE(*,*)'       total number of Trias  : ',nTrias
WRITE(*,*)'       total number of Quads  : ',nQuads
WRITE(*,*)'       total number of Nodes  : ',nVertices

!start and end indices for one zone
vS=1 ! Start index in Vertices(:,:)
vE=0 ! End index in Vertices(:,:)
tS=1 ! Start index Tria(:,:)
tE=0 ! End Index Tria(:,:)
qS=1 ! Start Index in Quad(:,:)
qE=0 ! End Index in Quad(:,:)
eS=1 ! Start index in BCedge(:,:)
eE=0 ! End index in BCEdge(:,:)
DO ZoneIndex=1,nZones

  WRITE(*,'(a,a)')'       reading zone : ',ZoneName(ZoneIndex)
  READ(ZoneName(Zoneindex), "(i3)")  zoneNumber
WRITE(*,*)'       number of Trias  : ',nTriasZone(ZoneIndex)
WRITE(*,*)'       number of Quads  : ',nQuadsZone(ZoneIndex)
WRITE(*,*)'       number of Nodes  : ',nVerticesZone(ZoneIndex)
  ! Get Nodes
  vS = vE + 1
  vE = vE + nVerticesZone(ZoneIndex)
  CALL cg_coord_read_f(CGNSUnit, BaseIndex, ZoneIndex, &
                       'CoordinateX', RealDouble, 1,   &
                       nVerticesZone(ZoneIndex),       &
                       Vertices(vS:vE,1),              & !Vertices x array
                       ierr)
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_coord_read_f')
  CALL cg_coord_read_f(CGNSUnit, BaseIndex, ZoneIndex, &
                       'CoordinateY', RealDouble, 1,   &
                       nVerticesZone(ZoneIndex),       &
                       Vertices(vS:vE,2),              & !Vertices y array
                       ierr)
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_coord_read_f')
  ! Get number of entries in Element connectivity
  CALL cg_section_read_f(CGNSUnit,     &
                         BaseIndex,    &
                         ZoneIndex,    &
                         SectionIndex, &
                         SectionName,  &   ! Name of section
                         SectionType,  &   ! Type of Element, in Gridgen only TRI_3 or QUAD_4
                         el_start,     &
                         el_end,       &
                         nbndry,       &
                         parent_flag,  &
                         ierr          )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_section_read_f')
  ! Allocate Data arrays for Mesh information
  SELECT CASE (SectionType)
  CASE(TRI_3) ! triangle zone
    tS = tE + 1
    tE = tE + nTriasZone(Zoneindex)
    DataSize= nTriasZone(Zoneindex)*3
  CASE(QUAD_4) ! quad zone
    qS = qE + 1
    qE = qE + nQuadsZone(Zoneindex)
    DataSize= nQuadsZone(Zoneindex)*4
  CASE DEFAULT
    WRITE(*,*) 'Error in ReadGridgenCGNS: Only Tri_3 and Quad_4 Elements allowed.'
    STOP
  END SELECT
  ALLOCATE(Elems(1:DataSize))
  ! Read Element connectivity
  CALL cg_elements_read_f(CGNSUnit,     &
                          BaseIndex,    &
                          ZoneIndex,    &
                          SectionIndex, &
                          Elems,        &  ! connectivity array
                          NULL,         &
                          ierr          )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_elements_read_f')
  ! write Element connectivity in Trias /Quads array
  vS=vS-1 ! index offset of vertices for multiple zones
  !check orientation (cross product) once for the whole zone!
  orient=  ( Vertices((Elems(2)+vS),X_DIR)-Vertices((Elems(1)+vS),X_DIR) ) &
         * ( Vertices((Elems(3)+vS),Y_DIR)-Vertices((Elems(1)+vS),Y_DIR) ) &
         - ( Vertices((Elems(2)+vS),Y_DIR)-Vertices((Elems(1)+vS),Y_DIR) ) &
         * ( Vertices((Elems(3)+vS),X_DIR)-Vertices((Elems(1)+vS),X_DIR) )
  IF (orient .LT. 0.) WRITE(*,'(A42,I3,A12)') ' !!!orientation of elements in zone ', zoneNumber, ' corrected!'
  SELECT CASE (SectionType)
  CASE(TRI_3) ! triangle zone
    Tria(tS:tE,1)=Elems(1:DataSize:3) + vS
    IF (orient .GT. 0.) THEN
      Tria(tS:tE,2)=Elems(2:DataSize:3) + vS
      Tria(tS:tE,3)=Elems(3:DataSize:3) + vS
    ELSE
      Tria(tS:tE,2)=Elems(3:DataSize:3) + vS
      Tria(tS:tE,3)=Elems(2:DataSize:3) + vS
    END IF
    Tria(tS:tE,4)=ZoneNumber
  CASE(QUAD_4)
    Quad(qS:qE,1)=Elems(1:DataSize:4) + vS
    Quad(qS:qE,3)=Elems(3:DataSize:4) + vS
    IF (orient .GT. 0.) THEN ! orientation correct
      Quad(qS:qE,2)=Elems(2:DataSize:4) + vS
      Quad(qS:qE,4)=Elems(4:DataSize:4) + vS
    ELSE !
      Quad(qS:qE,2)=Elems(4:DataSize:4) + vS
      Quad(qS:qE,4)=Elems(2:DataSize:4) + vS
    END IF
    Quad(qS:qE,5)=ZoneNumber
  CASE DEFAULT
    WRITE(*,*) 'Error in ReadGridgenCGNS: Only Tri_3 and Quad_4 Elements allowed.'
    STOP
  END SELECT
  ! Deallocate Elems array
  DEALLOCATE(Elems)
  ! Read number of boundaries
  CALL cg_nbocos_f(CGNSUnit,     &
                   BaseIndex,    &
                   ZoneIndex,    &
                   nbndry,       & ! number of boundaries in this zone
                   ierr          )
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_nbocos_f')
  ! Read boundary sides
  DO iBC = 1, nbndry
  ! Get boundary information
    CALL cg_boco_info_f(CGNSUnit,       &
                        BaseIndex,      &
                        ZoneIndex,      &
                        iBC,            &
                        BocoName,       & ! used to identify BC
                        BCType,         &
                        ptset_type,     &
                        nBCPnts,        & ! number of points
                        NormalIndex,    &
                        NormalListFlag, &
                        NormalDataType, &
                        nDataSets,      &
                        ierr            )
      IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_boco_info_f')
  ! Allocate storage for current boundary
    ALLOCATE(bcpnts(nBCPnts))
    eS = eE + 1
    eE = eE + (nBCPnts-1)
    CALL cg_boco_read_f(CGNSUnit,    &
                        BaseIndex,   &
                        ZoneIndex,   &
                        iBC,         &
                        bcpnts,      & ! edge vertices( sorted along BC edge)
                        normallist,  &
                        ierr         )
      IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_boco_read_f')
    READ(BocoName,'(I3)') BCCode
    WRITE(*,*)'       ...BC', BCCode,' read,'
    BCEdge(eS:eE, 1) = bcpnts(1:nBCPnts-1) + vS
    BCEdge(eS:eE, 2) = bcpnts(2:nBCPnts)   + vS
    BCEdge(eS:eE, 3) = BCCode
    DEALLOCATE(bcpnts)
  END DO ! iBC=1:nbdnry
  ! Zone connectivity
  ! number of interfaces in actual zone
  CALL cg_nconns_f(CGNSunit,        &
                   BaseIndex,       &
                   Zoneindex,       &
                   nConns,          & !number of connections
                   ierr)
    IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_nconns_f')
  DO iConn=1,nConns
    !
    CALL cg_conn_info_f(CGNSunit,        &
                        BaseIndex,       &
                        Zoneindex,       &
                        iConn,           &
                        ConnName,        &
                        location,        &
                        connect_type,    &
                        ptset_type,      &
                        nConnPnts,       &
                        donorname,       & ! name of donor zone
                        donor_zonetype,  &
                        donor_ptset_type,&
                        donor_datatype,  &
                        ndata_donor,     &
                        ierr)
      IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_conn_info_f')
    ! find zone number of donor zone for the vertex offset
    i=1
    vOffset=0
    DO WHILE (i .LE. nZones)
      IF( donorname .EQ. ZoneName(i)) EXIT
      vOffset=vOffset + nVerticesZone(i)
      i=i+1
    END DO
    IF (i .EQ. ZoneIndex) CYCLE ! zone Connectivity inside one zone
    ALLOCATE(ConnPnts(nConnPnts,2))
    CALL cg_conn_read_f(CGNSunit,         &
                        BaseIndex,        &
                        Zoneindex, iConn, &
                        ConnPnts(:,1),    & ! node Ids of actual zone
                        donor_datatype,   &
                        ConnPnts(:,2),    & ! node IDs of donor zone
                        ierr)
      IF (ierr .EQ. ERROR) CALL errorCGNS(ierr,'cg_conn_read_f')

    DO i=1,nConnPnts
      ZoneConnect(ConnPnts(i,1)+vS)=ConnPnts(i,2) + vOffset
    END DO
    DEALLocate(ConnPnts)
  END DO ! iConn=1,nConns
END DO ! ZoneIndex=1:nZones
CALL cg_close_f(CGNSUnit, ierr)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ReadGridgenCGNS

SUBROUTINE errorCGNS(ierr,routine)
!===================================================================================================================================
! errorCGNS: if an error occurs in ReadCGNS
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
  ! Include CGNS Library:
  ! (Please note that the CGNS library has to be installed in the computer's
  ! library and include path (see CGNS documentation for more information:
  ! www.cgns.org)
  USE CGNS
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER            :: ierr
CHARACTER(LEN=*)   :: routine
CHARACTER(LEN=255) :: errormsg
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

CALL cg_get_error_f(errormsg)
WRITE(*,*) 'error in ', TRIM(routine),' :'
WRITE(*,*) TRIM(errormsg)
STOP
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE errorCGNS

END MODULE MOD_Mesh_readin
