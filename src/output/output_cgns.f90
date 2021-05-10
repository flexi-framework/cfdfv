MODULE MOD_Output_cgns
!===================================================================================================================================
! Output moudels
!===================================================================================================================================
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Include CGNS Library:
! (Please note that the CGNS library has to be installed in the computer's
! library and include path (see CGNS documentation for more information:
! www.cgns.org)
INCLUDE 'cgnslib_f.h'
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CGNSOutput
   MODULE PROCEDURE CGNSOutput
END INTERFACE
INTERFACE FinalizeCGNSOutput
   MODULE PROCEDURE FinalizeCGNSOutput
END INTERFACE
INTERFACE CGNS_WriteMesh
   MODULE PROCEDURE CGNS_WriteMesh
END INTERFACE
!-----------------------------------------------------------------------------------------------------------------------------------
PUBLIC  :: CGNS_WriteMesh,CGNSOutput,FinalizeCGNSOutput 
!===================================================================================================================================

CONTAINS


SUBROUTINE CGNSOutput(FileName, time, iter, ExactSolution)
!===================================================================================================================================
! Write solution to CGNS file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ExactFunc    ,ONLY:ExactFunc
USE MOD_Equation_Vars,ONLY:intExactFunc
USE MOD_Mesh_Vars    ,ONLY:nNodes,nElems,FirstElem,GridFile
USE MOD_Mesh_Vars    ,ONLY:tElem
USE MOD_TimeDisc_Vars,ONLY:Time_Overall
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                   :: time
INTEGER,INTENT(IN)                :: iter
CHARACTER(LEN=200),INTENT(IN)     :: FileName
LOGICAL,INTENT(IN)                :: ExactSolution
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                 :: ierr, isize(1,3)
INTEGER                 :: CGNSUnit
INTEGER                 :: BaseIndex
INTEGER                 :: ZoneIndex
INTEGER                 :: SolutionIndex
INTEGER                 :: FieldIndex
CHARACTER(LEN=45)       :: text
CHARACTER(LEN=200)      :: SolName
TYPE(tElem), POINTER    :: aElem
REAL                    :: pvar(NVAR)
REAL, ALLOCATABLE       :: rho_array(:), v1_array(:), v2_array(:), v3_array(:), p_array(:) ! solution arrays
!===================================================================================================================================

! Open the solution file
CALL cg_set_file_type_f(CG_FILE_ADF2, ierr)
CALL cg_open_f(TRIM(FileName), MODE_WRITE, CGNSUnit, ierr)
!-----------------------------------------------------------------------------------------------------------------------------------
! Set up Data for CGNS
isize(1,1) = nNodes
isize(1,2) = nElems
isize(1,3) = 0
! Create Base
!-----------------------------------------------------------------------------------------------------------------------------------
CALL cg_base_write_f(CGNSUnit,       &
                     'Base',         &
                     2,              &
                     3,              &
                     BaseIndex,      &
                     ierr            )
! Create zone
!-----------------------------------------------------------------------------------------------------------------------------------
CALL cg_zone_write_f(CGNSUnit,       &
                     BaseIndex,      &
                     'Zone',         &
                     isize,          &
                     Unstructured,   &
                     ZoneIndex,      &
                     ierr            )
!-----------------------------------------------------------------------------------------------------------------------------------
! Link Vertices and connectivity from the CGNS grid file
CALL cg_goto_f(CGNSUnit,    &
               BaseIndex,   &
               ierr,        &
               'Zone_t',    &
               ZoneIndex,   &
               'end'        )
CALL cg_link_write_f('GridCoordinates',            &
                     TRIM(GridFile),               &
                     '/Base/Zone/GridCoordinates', &
                     ierr                          )
CALL cg_link_write_f('Elements',                   &
                     TRIM(GridFile),               &
                     '/Base/Zone/Elements',        &
                     ierr                          )
!-----------------------------------------------------------------------------------------------------------------------------------
! Create new solution zone
SolName = 'FlowSolution'
CALL cg_sol_write_f(CGNSUnit,      &
                    BaseIndex,     &
                    ZoneIndex,     &
                    TRIM(SolName), &
                    CellCenter,    &
                    SolutionIndex, &
                    ierr           )
!-----------------------------------------------------------------------------------------------------------------------------------
! Prepare Density Array, x-Velocity Array, y-Velocity Array & pressure Array
!-----------------------------------------------------------------------------------------------------------------------------------
! Allocate the solution arrays
ALLOCATE(rho_array(1:nElems))
ALLOCATE(v1_array(1:nElems))
ALLOCATE(v2_array(1:nElems))
ALLOCATE(v3_array(1:nElems))	! dummy_Array for compatibility with Paraview 5.6.1
ALLOCATE(p_array(1:nElems))
!-----------------------------------------------------------------------------------------------------------------------------------
! Save the solution in a CGNS-compatible format
IF (ExactSolution) THEN
  aElem => FirstElem
  DO WHILE(ASSOCIATED(aElem))
    CALL ExactFunc(intExactFunc, aElem%Bary, time, pvar)
    rho_array(aElem%id) = pvar(RHO)
    v1_array(aElem%id)  = pvar(V1)
    v2_array(aElem%id)  = pvar(V2)
    p_array(aElem%id)   = pvar(P)
    aElem=>aElem%nextElem
  END DO
ELSE
  aElem => FirstElem
  DO WHILE(ASSOCIATED(aElem))
    rho_array(aElem%id) = aElem%pvar(RHO)
    v1_array(aElem%id)  = aElem%pvar(V1)
    v2_array(aElem%id)  = aElem%pvar(V2)
    p_array(aElem%id)   = aElem%pvar(P)
    aElem=>aElem%nextElem
  END DO
END IF
! Write density (RHO)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL cg_field_write_f(CGNSUnit,      &
                      BaseIndex,     &
                      ZoneIndex,     &
                      SolutionIndex, &
                      RealDouble,    &
                      'Density',     &
                      rho_array,     &
                      FieldIndex,    &
                      ierr           )
! Write x-Velocity (V1)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL cg_field_write_f(CGNSUnit,      &
                      BaseIndex,     &
                      ZoneIndex,     &
                      SolutionIndex, &
                      RealDouble,    &
                      'VelocityX',   &
                      v1_array,      &
                      FieldIndex,    &
                      ierr           )
! Write y-Velocity (V2)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL cg_field_write_f(CGNSUnit,      &
                      BaseIndex,     &
                      ZoneIndex,     &
                      SolutionIndex, &
                      RealDouble,    &
                      'VelocityY',   &
                      v2_array,      &
                      FieldIndex,    &
                      ierr           )
! Write y-Velocity (V3)
!-----------------------------------------------------------------------------------------------------------------------------------
v3_array(:) = 0.0
CALL cg_field_write_f(CGNSUnit,      &
                      BaseIndex,     &
                      ZoneIndex,     &
                      SolutionIndex, &
                      RealDouble,    &
                      'VelocityZ',   &
                      v3_array,      &
                      FieldIndex,    &
                      ierr           )
! Write Pressure (P)
!-----------------------------------------------------------------------------------------------------------------------------------
CALL cg_field_write_f(CGNSUnit,      &
                      BaseIndex,     &
                      ZoneIndex,     &
                      SolutionIndex, &
                      RealDouble,    &
                      'Pressure',    &
                      p_array,       &
                      FieldIndex,    &
                      ierr           )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write Convergence information
CALL cg_goto_f(CGNSUnit, BaseIndex, ierr, 'end')
WRITE(text, '(F20.12,1X,F20.12)') time, Time_Overall
CALL cg_descriptor_write_f('ConvergenceInfo', TRIM(text), ierr)
!-----------------------------------------------------------------------------------------------------------------------------------
! Close file
CALL cg_close_f(CGNSUnit,ierr)
!-----------------------------------------------------------------------------------------------------------------------------------
! Deallocate Arrays
DEALLOCATE(rho_array)
DEALLOCATE(V1_array)
DEALLOCATE(V2_array)
DEALLOCATE(V3_array)
DEALLOCATE(p_array)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CGNSOutput



SUBROUTINE FinalizeCGNSOutput()
!===================================================================================================================================
! Finialize CGNSoutput
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars,ONLY:tOutputTime
USE MOD_Output_Vars,ONLY:strOutFile,OutputTimes
USE MOD_Mesh_Vars  ,ONLY:GridFile,nElems,nNodes
USE MOD_TimeDisc_Vars,ONLY:stationary
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tOutputTime), POINTER :: OutputTime
INTEGER :: nOutput = 0
INTEGER :: i
REAL, ALLOCATABLE :: times(:)
INTEGER, ALLOCATABLE :: iters(:)
INTEGER :: ierr, isize(1,3)
INTEGER :: CGNSUnit
INTEGER :: BaseIndex
INTEGER :: ZoneIndex
CHARACTER(LEN=15) :: cIter
CHARACTER(LEN=32), ALLOCATABLE :: solutionNames(:)
CHARACTER(LEN=255) :: FileName
!===================================================================================================================================

! Count number of data outputs
nOutput=0
OutputTime => OutputTimes
DO WHILE ( ASSOCIATED(OutputTime) )
  nOutput = nOutput + 1
  OutputTime => OutputTime%Next
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Allocate and fill arrays
ALLOCATE(times(nOutput))
ALLOCATE(iters(nOutput))
ALLOCATE(solutionNames(nOutput))

i = nOutput
OutputTime => OutputTimes
DO WHILE ( ASSOCIATED(OutputTime) )
  times(i) = OutputTime%Time
  iters(i) = OutputTime%Iter
  WRITE(solutionNames(i), '(A,I9.9)') 'FlowSolution', OutputTime%Iter
  i = i - 1
  OutputTime => OutputTime%Next
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Open CGNS file
CALL cg_set_file_type_f(CG_FILE_ADF2, ierr)
CALL cg_open_f(TRIM(strOutFile) // '_Master.cgns', MODE_WRITE, CGNSUnit, ierr)
!-----------------------------------------------------------------------------------------------------------------------------------
! Set up Data for CGNS
isize(1,1) = nNodes
isize(1,2) = nElems
isize(1,3) = 0
!-----------------------------------------------------------------------------------------------------------------------------------
! Create Base
CALL cg_base_write_f(CGNSUnit,       &
                     'Base',         &
                     2,              &
                     3,              &
                     BaseIndex,      &
                     ierr            )
!-----------------------------------------------------------------------------------------------------------------------------------
! Create zone
CALL cg_zone_write_f(CGNSUnit,       &
                     BaseIndex,      &
                     'Zone',         &
                     isize,          &
                     Unstructured,   &
                     ZoneIndex,      &
                     ierr            )
!-----------------------------------------------------------------------------------------------------------------------------------
! Link Vertices and connectivity from the CGNS grid file
CALL cg_goto_f(CGNSUnit,    &
               BaseIndex,   &
               ierr,        &
               'Zone_t',    &
               ZoneIndex,   &
               'end'        )
CALL cg_link_write_f('GridCoordinates',            &
                     TRIM(GridFile),               &
                     '/Base/Zone/GridCoordinates', &
                     ierr                          )
CALL cg_link_write_f('Elements',                   &
                     TRIM(GridFile),               &
                     '/Base/Zone/Elements',        &
                     ierr                          )
!-----------------------------------------------------------------------------------------------------------------------------------
! Link solutions
CALL cg_goto_f(CGNSUnit, BaseIndex, ierr, 'Zone_t', ZoneIndex, 'end')
DO i = 1, nOutput
  IF (stationary) THEN
    WRITE(cIter, '(I9.9)') iters(i)
    CALL cg_link_write_f(solutionNames(i), &
                        TRIM(strOutFile) // '_' // TRIM(cIter) // '.cgns', &
                        '/Base/Zone/FlowSolution', ierr)
   ELSE
     FileName = TIMESTAMP(strOutFile,times(i))
     CALL cg_link_write_f(solutionNames(i), &
                        TRIM(FileName) // '.cgns', &
                        '/Base/Zone/FlowSolution', ierr)
   END IF 
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Create BaseIter Node
CALL cg_biter_write_f(CGNSUnit,         &
                      BaseIndex,        &
                      'TimeIterValues', &
                      nOutput,          &
                      ierr              )
CALL cg_goto_f(CGNSUnit,                &
               BaseIndex,               &
               ierr,                    &
               'BaseIterativeData_t',   &
               1,                       &
               'end'                    )
CALL cg_array_write_f("TimeValues",   &
                      RealDouble,     &
                      1,              &
                      (/ nOutput /),  &
                      times,          &
                      ierr            )
CALL cg_array_write_f("IterationValues",  &
                      Integer,            &
                      1,                  &
                      (/ nOutput /),      &
                      iters,              &
                      ierr                )
!-----------------------------------------------------------------------------------------------------------------------------------
! Create ZoneIter Node
CALL cg_ziter_write_f(CGNSUnit,             &
                      BaseIndex,            &
                      ZoneIndex,            &
                      'ZoneIterativeData',  &
                      ierr                  )
CALL cg_goto_f(CGNSUnit,              &
               BaseIndex,             &
               ierr,                  &
               'Zone_t',              &
               ZoneIndex,             &
               'ZoneIterativeData_t', &
               1,                     &
               'end'                  )
CALL cg_array_write_f('FlowSolutionPointers', &
                      Character,              &
                      2,                      &
                      (/ 32, nOutput /),     &
                      solutionNames,              &
                      ierr                    )
!-----------------------------------------------------------------------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------------------------------------------------------------------
CALL cg_close_f(CGNSUnit,ierr)
!-----------------------------------------------------------------------------------------------------------------------------------
! Deallocate data
DEALLOCATE(times, iters, solutionNames)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE FinalizeCGNSOutput


SUBROUTINE CGNS_WriteMesh()
!===================================================================================================================================
! Write CGNS mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Boundary_Vars,  ONLY:nBC,isPeriodic,FirstBC
USE MOD_Boundary_Vars,  ONLY:tBoundary
USE MOD_Mesh_Vars,      ONLY:tNode,tElem,tPureSidePtr
USE MOD_Mesh_Vars,      ONLY:nBCSides,firstNode,firstElem,GridFile
USE MOD_Mesh_Vars,      ONLY:FirstBCSide,nNodes,nElems,nTrias,nQuads
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Include CGNS Library:
! (Please note that the CGNS library has to be installed in the computer's
! library and include path (see CGNS documentation for more information:
! www.cgns.org)
INCLUDE 'cgnslib_f.h'
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                     :: ierr, nCount, iNode
INTEGER                     :: isize(1,3)
INTEGER                     :: CGNSUnit
INTEGER                     :: BaseIndex
INTEGER                     :: ZoneIndex
INTEGER                     :: GridIndex
INTEGER                     :: BoundaryIndex
INTEGER                     :: CoordinateIndex
INTEGER                     :: SectionIndex
INTEGER, POINTER            :: Elems(:), bElems(:)

REAL, POINTER               :: Nodes(:,:)
INTEGER                     :: allocStat, ElemIntSize, n
INTEGER                     :: BCPartition(nBC+1), nSide, iBC, ipnts(nBCSides)

TYPE(tNode), POINTER        :: aNode
TYPE(tElem), POINTER        :: aElem

TYPE(tPureSidePtr), POINTER :: aBCSide
TYPE(tBoundary),    POINTER :: aBC
CHARACTER(LEN=15)           :: cBC
!===================================================================================================================================

! Set up Data
isize(1,1) = nNodes
isize(1,2) = nElems
isize(1,3) = 0
!-----------------------------------------------------------------------------------------------------------------------------------
! Determine the size of the element array and allocate it
ElemIntSize = nTrias * 4 + nQuads * 5
ALLOCATE(Elems(1:ElemIntSize), STAT = allocStat)
ALLOCATE(Nodes(1:nNodes, 3), STAT = allocStat)
ElemIntSize = nBCSides * 3
ALLOCATE(bElems(1:ElemIntSize), STAT = allocStat)
!-----------------------------------------------------------------------------------------------------------------------------------
! Save vertices in a CGNS compatible format
iNode=1
aNode => firstNode
DO WHILE(ASSOCIATED(aNode))
  Nodes(iNode,:) = (/ aNode%x(:), 0.0 /)	! 3D Point extension for compatibility with Paraview 5.6.1
  aNode%ID=iNode
  aNode => aNode%NextNode
  iNode=iNode+1
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Save element connectivity in a CGNS compatible format
nCount = 1
aElem => FirstElem
DO WHILE (ASSOCIATED(aElem))
  IF (aElem%ElemType == 3) THEN
    Elems(nCount) = TRI_3
  ELSEIF (aElem%ElemType == 4) THEN
    Elems(nCount) = QUAD_4
  END IF
  nCount = nCount + 1
  DO iNode = 1, aElem%ElemType
    Elems(nCount) = aElem%NodeArray(iNode)%Node%ID
    nCount = nCount + 1
  END DO
  aElem => aElem%nextElem
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! Save Boundary Elements
IF (.NOT.isPeriodic) THEN
  BCPartition(1) = nElems
  nSide = nElems
  nCount = 1
  iBC = 2
  aBC => FirstBC
  DO WHILE(ASSOCIATED(aBC))
    aBCSide => FirstBCSide
    DO WHILE (ASSOCIATED(aBCSide))
      IF (ASSOCIATED(aBCSide%Side%BC,aBC)) THEN
        bElems(nCount)   = BAR_2
        bElems(nCount+1) = aBCSide%Side%Node(1)%Node%ID
        bElems(nCount+2) = aBCSide%Side%Node(2)%Node%ID
        nCount = nCount + 3
        nSide  = nSide  + 1
      END IF
      aBCSide => aBCSide%nextSide
    END DO
    BCPartition(iBC) = nSide
    iBC = iBC + 1
    aBC => aBC%NextBC
  END DO
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Open the CGNS grid file for writing:
CALL cg_open_f(TRIM(GridFile),           &
               MODE_WRITE,               &
               CGNSUnit,                 &
               ierr                      )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write coordinate base to CGNS file
CALL cg_base_write_f(CGNSUnit,  &
                     'Base',    &
                     2,         &
                     3,         &
                     BaseIndex, &
                     ierr       )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write the computational zone to the CGNS file
CALL cg_zone_write_f(CGNSUnit,     &
                     BaseIndex,    &
                     'Zone',       &
                     isize,        &
                     Unstructured, &
                     ZoneIndex,    &
                     ierr          )
call cg_grid_write_f(CGNSUnit,          &
                     BaseIndex,         &
                     ZoneIndex,         &
                     'GridCoordinates', &
                     GridIndex,         &
                     ierr               )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write the vertices' x-coordinates to the file
CALL cg_coord_write_f(CGNSUnit,              &
                      BaseIndex,             &
                      ZoneIndex,             &
                      RealDouble,            &
                      'CoordinateX',         &
                      Nodes(: ,X_DIR),       &
                      CoordinateIndex,       &
                      ierr                   )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write the vertices' y-coordinates to the file
CALL cg_coord_write_f(CGNSUnit,              &
                      BaseIndex,             &
                      ZoneIndex,             &
                      RealDouble,            &
                      'CoordinateY',         &
                      Nodes(: ,Y_DIR),       &
                      CoordinateIndex,       &
                      ierr                   )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write the vertices' z-coordinates to the file -> dummy Points
CALL cg_coord_write_f(CGNSUnit,              &
                      BaseIndex,             &
                      ZoneIndex,             &
                      RealDouble,            &
                      'CoordinateZ',         &
                      Nodes(: ,3),			 &
                      CoordinateIndex,       &
                      ierr                   )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write the element connectivity to the CGNS file
CALL cg_section_write_f(CGNSUnit,     &
                        BaseIndex,    &
                        ZoneIndex,    &
                        'Elements',   &
                        MIXED,        &
                        1,            &
                        nElems,       &
                        0,            &
                        Elems,        &
                        SectionIndex, &
                        ierr          )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write the boundary connectivity to the CGNS file
IF (.NOT.isPeriodic) THEN
  CALL cg_section_write_f(CGNSUnit,                    &
                          BaseIndex,                   &
                          ZoneIndex,                   &
                          'Boundaries',                &
                          MIXED,                       &
                          nElems + 1,                  &
                          nElems + nBCSides,           &
                          0,                           &
                          bElems,                      &
                          SectionIndex,                &
                          ierr                         )
!-----------------------------------------------------------------------------------------------------------------------------------
! Write the boundary connectivity to the CGNS file
  iBC = 1
  aBC => FirstBC
  DO WHILE (ASSOCIATED(aBC))
    WRITE(cBC, '(I3)') aBC%BCType * 100 + aBC%BCID
    nCount = 0
    DO n = BCPartition(iBC)+1, BCPartition(iBC+1)
      nCount = nCount + 1
      ipnts(nCount) = n
    ENDDO
    CALL cg_boco_write_f(CGNSUnit,       &
                         BaseIndex,      &
                         ZoneIndex,      &
                         TRIM(cBC),      &
                         BCGeneral,      &
                         ElementList,    &
                         nCount,         &
                         ipnts,          &
                         BoundaryIndex,  &
                         ierr            )
    aBC => aBC%nextBC
    iBC = iBC+1
  END DO
END IF
! Close CGNS file
CALL cg_close_f(CGNSUnit, ierr)
! Deallocate temporary arrays
DEALLOCATE(Elems, Nodes, STAT = allocStat)
!-----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE CGNS_WriteMesh

FUNCTION TIMESTAMP(Filename,Time)
!===================================================================================================================================
! Creates a timestamp, consistent of a filename (project name + processor) and current time niveau
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: Filename  ! (file)name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=255)             :: TimeStamp ! the complete timestamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                        :: i         ! loop variable
REAL                           :: Time      ! time
!===================================================================================================================================

WRITE(TimeStamp,'(F15.7)')Time
! Replace spaces with 0's
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i).EQ.' ') TimeStamp(i:i)='0'
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)
!-----------------------------------------------------------------------------------------------------------------------------------
END FUNCTION TIMESTAMP

END MODULE MOD_Output_cgns
