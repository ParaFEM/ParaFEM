PROGRAM d2off

  !/****h* /d2off
  !*  NAME
  !*    PROGRAM: d2off
  !*  SYNOPSIS
  !*    Usage:   ./d2off <file>
  !*  FUNCTION
  !*    Converts ParaFEM (.d) files into Object File Format (.off) files
  !*    
  !*    The documentation below was taken from:
  !*
  !*    http://shape.cs.princeton.edu/benchmark/documentation/off_format.html
  !*
  !*    Object File Format (.off) files are used to represent the geometry of
  !*    a model by specifying the polygons of the model's surface. The polygons
  !*    can have any number of vertices.
  !*    
  !*    The .off files used in ParaFEM conform to the following standard. OFF
  !*    files are all ASCII files beginning with the keyword OFF. The next line
  !*    states the number of vertices, the number of faces and the number of 
  !*    edges. The number of edges can be safely ignored.
  !*
  !*    The vertices are listed with x, y, z coordinates, written one per line. 
  !*    After te list of vertice, the faces are listed, with one face per line.
  !*    For each face, the number of vertices in specified, followed by indices
  !*    into the list of vertices. See the examples below.
  !*
  !*    OFF
  !*    numVertices numFaces numEdges
  !*    x1 y1 z1
  !*    x2 y2 z2
  !*    ...
  !*    xn xn xn
  !*    nvertices v1 v2 v3 ... vn
  !*    mvertices v1 v2 v3 ... vm
  !*    ... 
  !*    
  !*    A simple example for a cube:
  !*    
  !*    OFF
  !*    8 6 0
  !*    -0.5  -0.5   0.5
  !*     0.5  -0.5   0.5
  !*    -0.5   0.5   0.5
  !*     0.5   0.5   0.5
  !*    -0.5   0.5  -0.5
  !*     0.5   0.5  -0.5
  !*    -0.5  -0.5  -0.5
  !*     0.5  -0.5  -0.5
  !*     4 0 1 3 2
  !*     4 2 3 5 4
  !*     4 4 5 7 6
  !*     4 6 7 1 0
  !*     4 1 7 5 3
  !*     4 6 0 2 4
  !*
  !*  NOTE: The current version only works with 4 node tetrahedra. It also
  !*        creates repeated faces. The output OFF file can be viewed, cleaned
  !*        up and converted to other graphics formats using "Meshlab". See:
  !*
  !*        http://meshlab.sourceforge.net (site accessed 26.05.2010)
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  USE precision
  USE global_variables
  USE timing

  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER     :: ndim  = 3
  INTEGER               :: nod
  INTEGER               :: nn 
  INTEGER               :: nels
  INTEGER               :: nr 
  INTEGER               :: nodof
  INTEGER               :: i,j,k,id,iel
  INTEGER,PARAMETER     :: one   = 1
  INTEGER,PARAMETER     :: three = 3 
  INTEGER               :: nedges 
  INTEGER               :: argc,iargc
  INTEGER               :: meshgen
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp
  CHARACTER(LEN=50)     :: job_name
  CHARACTER(LEN=50)     :: fname
  CHARACTER(LEN=15)     :: element

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  INTEGER,ALLOCATABLE   :: g_num(:,:)
  REAL(iwp),ALLOCATABLE :: g_coord(:,:)
  REAL(iwp),ALLOCATABLE :: timest(:)

!------------------------------------------------------------------------------
! 3. Read job_name from the command line.
!    Read control data and mesh data
!------------------------------------------------------------------------------

  ALLOCATE(timest(1))
  timest    = zero
  timest(1) = elap_time()

  argc = iargc()

  IF(argc /= 1) THEN
    WRITE(*,*)
    WRITE(*,'(A)') " Usage: d2off <job_name> "
    WRITE(*,*)
    WRITE(*,'(A)') "        d2off expects argument <job_name> and 2 input files"
    WRITE(*,'(A)') "        <job_name>.dat and <job_name>.d"
    WRITE(*,*)
    WRITE(*,'(A)') "********************* PROGRAM ABORTED *********************"
    WRITE(*,*)
    STOP 
  END IF

  CALL GETARG(1, job_name)

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
  READ(10,*) element,meshgen,nels,nn,k,k,nod
  CLOSE(10)

  ALLOCATE(g_num(nod,nels))
  ALLOCATE(g_coord(ndim,nn))

  g_num   = 0
  g_coord = zero

!------------------------------------------------------------------------------
! 4. Read ".d" geometry file
!------------------------------------------------------------------------------

  ALLOCATE(g_num(nod,nels))
  ALLOCATE(g_coord(ndim,nn))
  g_num   = 0
  g_coord = zero

  fname = job_name(1:INDEX(job_name, " ") -1) // ".d"
  OPEN(11,FILE=fname,STATUS='OLD',ACTION='READ')

  READ(11,*)
  READ(11,*)

  DO i = 1,nn
    READ(11,*) id, g_coord(:,i)    
  END DO

  READ(11,*)
  DO iel = 1,nels 
    READ(11,*) id, k, k, k, g_num(:,iel)
  END DO

!------------------------------------------------------------------------------
! 5. Write ".off" geometry file
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ") -1) // ".off"
  OPEN(20,FILE=fname,STATUS='REPLACE',ACTION='WRITE') 
  
  nedges = 0
 
  WRITE(20,'(A)') "OFF"
  WRITE(20,'(3I8)') nn, nels*4, nedges 

  DO i = 1, nn
    WRITE(20,'(3E20.12)') g_coord(:,i)
  END DO

  DO iel = 1, nels
     WRITE(20,'(I3,3I8)') three, g_num(1,iel)-1, g_num(2,iel)-1, g_num(3,iel)-1
     WRITE(20,'(I3,3I8)') three, g_num(1,iel)-1, g_num(4,iel)-1, g_num(2,iel)-1
     WRITE(20,'(I3,3I8)') three, g_num(2,iel)-1, g_num(4,iel)-1, g_num(3,iel)-1
     WRITE(20,'(I3,3I8)') three, g_num(3,iel)-1, g_num(4,iel)-1, g_num(1,iel)-1
  END DO
  
END PROGRAM d2off
