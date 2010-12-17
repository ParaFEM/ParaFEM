MODULE INPUT

  !/****h* /input
  !*  NAME
  !*    MODULE: input
  !*  SYNOPSIS
  !*    Usage:      USE input
  !*  FUNCTION
  !*    Contains subroutines that handle input data. These subroutines are 
  !*    parallel and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    READ_G_COORD_PP        Reads the global coordinates
  !*    READ_G_NUM_PP          Reads the element nodal steering array
  !*    READ_LOADS             Reads nodal forces
  !*    READ_LOADS_NS          Reads lid velocities for p126
  !*    READ_FIXED             Reads fixed freedoms for displacement control
  !*    READ_REST              Reads the restraints
  !*    READ_P121              Reads the control data for program p121
  !*    READ_P126              Reads the control data for program p126
  !*    READ_P128AR            Reads the control data for program p128ar
  !*    READ_P129              Reads the control data for program p129
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  USE precision
  USE mp_interface

  CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_G_COORD_PP(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)

  !/****f* input/read_g_coord_pp
  !*  NAME
  !*    SUBROUTINE: read_g_coord_pp
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,     &
  !*                                     g_coord_pp)
  !*  FUNCTION
  !*    
  !*  INPUTS
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    21 May 2008
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010
  !******
  !* The algorithm needs to be improved to make best use of the available
  !* memory resources. Perhaps there could be a switch for problem size.
  !* When the problem is smaller than a certain threshold, we could just
  !* read and distribute the global array. 
  !*
  !* As this subroutine repeatedly goes around an elements and nodes 
  !* loop, it is going to do a lot of unnecessary work. Expect this routine
  !* to be very time consuming. It will benefit from optimisation.
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50),INTENT(IN) :: job_name
  CHARACTER(LEN=50)            :: fname
  INTEGER, INTENT(IN)          :: nn, npes, numpe
  INTEGER, INTENT(IN)          :: g_num_pp(:,:)
  INTEGER                      :: ndim
  INTEGER                      :: nels_pp
  INTEGER                      :: nod
  INTEGER                      :: iel, i, j, k, l ! loop counters
  INTEGER                      :: bitBucket
  REAL(iwp), INTENT(INOUT)     :: g_coord_pp(:,:,:) 
  REAL(iwp), ALLOCATABLE       :: g_coord(:,:)    ! temporary array
  REAL(iwp)                    :: zero = 0.0_iwp

!------------------------------------------------------------------------------
! 1. Allocate temporary array
!------------------------------------------------------------------------------

  nod      = UBOUND(g_coord_pp,1)
  ndim     = UBOUND(g_coord_pp,2)
  nels_pp  = UBOUND(g_coord_pp,3)

  ALLOCATE(g_coord(ndim,nn))

!------------------------------------------------------------------------------
! 2. Open the data file
!------------------------------------------------------------------------------

  fname     = job_name(1:INDEX(job_name, " ")-1) // ".d"  
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')

!------------------------------------------------------------------------------
! 3. Read data and populate g_coord_pp
!------------------------------------------------------------------------------

  g_coord_pp= zero

  READ(10,*)   !headers
  READ(10,*)   !headers

  DO i = 1,nn
    READ(10,*) bitBucket,g_coord(:,j)
  END DO
  
  DO iel = 1, nels_pp
    DO k = 1, nod
      l                  = g_num_pp(k,iel)
      g_coord_pp(k,:,iel)= g_coord(:,l)
    END DO
  END DO
  
!------------------------------------------------------------------------------
! 4. Deallocate global arrays and close data file
!------------------------------------------------------------------------------

  DEALLOCATE(g_coord)    
  CLOSE(10)

  RETURN
  
  END SUBROUTINE READ_G_COORD_PP
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_G_NUM_PP(job_name,iel_start,nn,npes,numpe,g_num_pp)

  !/****f* input/read_g_num_pp
  !*  NAME
  !*    SUBROUTINE: read_g_num_pp
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_num_pp(job_name,iel_start,nn,npes,numpe,     &
  !*                                   g_num_pp)
  !*  FUNCTION
  !*    Reads the element nodal steering array. Serial version.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    job_name              : Character
  !*                          : Used to create file name to read
  !*
  !*    iel_start             : Integer
  !*                          : First element number in a process
  !*
  !*    nn                    : Integer
  !*                          : Total number of nodes 
  !*
  !*    npes                  : Integer
  !*                          : Total number of processors
  !*
  !*    numpe                 : Integer
  !*                          : Process number
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    g_num_pp(nod,nels_pp) : Integer
  !*                          : Elements connectivity
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname
  INTEGER, INTENT(IN)           :: iel_start, nn, npes, numpe
  INTEGER, INTENT(INOUT)        :: g_num_pp(:,:)
  INTEGER                       :: nod, nels_pp, iel, i, k
  INTEGER                       :: ndim=3

!------------------------------------------------------------------------------
! 1. Initiallize variables
!------------------------------------------------------------------------------

  nod       = UBOUND(g_num_pp,1)
  nels_pp   = UBOUND(g_num_pp,2)

!------------------------------------------------------------------------------
! 2. Open the data file and advance to the start of the element steering array
!------------------------------------------------------------------------------

  fname     = job_name(1:INDEX(job_name, " ")-1) // ".d"
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')

  READ(10,*)   !header
  READ(10,*)   !header

  DO i = 1,nn 
    READ(10,*) !read nodes until reaching the elements 
  END DO

  READ(10,*)   !keyword

!------------------------------------------------------------------------------
! 3. Read element steering array
!------------------------------------------------------------------------------

  DO iel = 1,nels_pp
    READ(10,*)k,k,k,k,g_num_pp(:,iel),k
  END DO
        
  CLOSE(10)

  RETURN

  END SUBROUTINE READ_G_NUM_PP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_LOADS(JOB_NAME,NUMPE,NODE,VALUE)

  !/****f* input/read_loads
  !*  NAME
  !*    SUBROUTINE: read_loads
  !*  SYNOPSIS
  !*    Usage:      CALL read_loads(job_name,numpe,node,value)
  !*  FUNCTION
  !*    Master processor reads the global array of nodal forces 
  !*    and broadcasts them to the slave processors.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    job_name                 : Character
  !*                             : Used to generate file name to read
  !*
  !*    numpe                    : Integer
  !*                             : Processor number used for I/O
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    node(loaded_nodes)       : Integer
  !*                             : Nodes that have an applied load 
  !*
  !*    value(ndim,loaded_nodes) : Real
  !*                             : Force applied on each node
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010
  !******
  !*  The method is based on reading a global array, although it's not a
  !*  worry because the nodes with loads are a small proportion of the 
  !*  whole mesh.
  !*
  !*/
  
  IMPLICIT NONE
  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname  
  INTEGER                       :: i
  INTEGER,INTENT(IN)            :: numpe
  INTEGER,INTENT(INOUT)         :: node(:)
  REAL(iwp),INTENT(INOUT)       :: value(:,:)

  fname = job_name(1:INDEX(job_name, " ")-1) // ".lds"
  OPEN(22, FILE=fname, STATUS='OLD', ACTION='READ')

  DO i = 1,ubound(node,1)   !ubound(node,1)=loaded_nodes 
    READ(22,*) node(i),value(:,i)
  END DO

  CLOSE(22)

  RETURN

  END SUBROUTINE READ_LOADS

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE READ_FIXED(job_name,numpe,node,sense,valf)

    !/****f* input/read_fixed
    !*  NAME
    !*    SUBROUTINE: read_fixed
    !*  SYNOPSIS
    !*    Usage:      CALL read_fixed(job_name,numpe,node,sense,valf)
    !*  FUNCTION
    !*    Master process reads the global array of nodes with fixed degrees
    !*    of freedom, the fixed degrees of freedom and the value.
    !*    Master process broadcasts to slave processes.
    !*
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    job_name               : Character
    !*                           : File name to read
    !*
    !*    numpe                  : Integer
    !*                           : Process number
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    sense(fixed_freedoms)  : Integer
    !*                           : Degree of freedom fixed
    !*
    !*    node(fixed_freedoms)   : Integer
    !*                           : Node fixed
    !*
    !*    valf(fixed_freedoms)   : Real
    !*                           : Value fixed
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    25.10.2010
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  The method is based on reading a global array, although it's not a
    !*  worry because the nodes with boundary conditions are a small
    !*  proportion of the whole mesh.
    !*
    !*  Variable names need to be changed to those used in Smith & Griffiths
    !*/

    IMPLICIT NONE

    CHARACTER(LEN=50), INTENT(IN) :: job_name
    CHARACTER(LEN=50)             :: fname  
    INTEGER,INTENT(IN)            :: numpe
    INTEGER,INTENT(OUT)           :: sense(:), node(:)
    REAL(iwp),INTENT(OUT)         :: valf(:)
    INTEGER                       :: i, fixed_freedoms, ier
 
    fixed_freedoms = UBOUND(sense,1)

    fname = job_name(1:INDEX(job_name, " ")-1) // ".fix"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')

    DO i = 1,fixed_freedoms
      READ(10,*)node(i),sense(i),valf(i)
    END DO

    CLOSE(10)

    RETURN

  END SUBROUTINE READ_FIXED
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE READ_REST(JOB_NAME,NUMPE,REST)

  !/****f* input/read_rest
  !*  NAME
  !*    SUBROUTINE: read_rest
  !*  SYNOPSIS
  !*    Usage:      CALL read_rest(job_name,numpe,rest)
  !*  FUNCTION
  !*    Master process reads the global array of nodes with restrained 
  !*    degrees of freedom.
  !*    Master process broadcasts to slave processes.
  !*
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    job_name               : Character
  !*                           : Used to create file name to read
  !*
  !*    numpe                  : Integer
  !*                           : Process number
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    rest(nr,nodof+1)       : Integer
  !*                           : Nodes (column 1) and degrees of freedom
  !*                             restrained (columns 2,3,4). 
  !*                             The criterion is:  0 fixed, ! free
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/   
  
  IMPLICIT NONE
  
  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname
  INTEGER                       :: i
  INTEGER,INTENT(IN)            :: numpe
  INTEGER,INTENT(INOUT)         :: rest(:,:)

  fname     = job_name(1:INDEX(job_name, " ")-1) // ".bnd"
  OPEN(23, FILE=fname, STATUS='OLD', ACTION='READ')

  DO i = 1, UBOUND(rest,1)
    READ(23,*) rest(i,:)
  END DO

  CLOSE(23)

  RETURN
  END SUBROUTINE READ_REST

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P121(job_name,numpe,e,element,fixed_freedoms,limit,        &
                       loaded_nodes,mesh,nels,nip,nn,nod,nr,tol,v)

  !/****f* input/read_p121
  !*  NAME
  !*    SUBROUTINE: read_p121
  !*  SYNOPSIS
  !*    Usage:      CALL read_p121(job_name,numpe,e,element,fixed_freedoms,
  !*                               limit,loaded_nodes,mesh,nels,nip,nn,nod,nr,
  !*                               tol,v)
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    job_name               : Character
  !*                           : File name that contains the data to be read
  !*
  !*    numpe                  : Integer
  !*                           : Processor number
  !*
  !*    The following arguments have the INTENT(INOUT) attribute:
  !*
  !*    e                      : Real
  !*                           : Young's modulus
  !*
  !*    element                : Character
  !*                           : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*    fixed_freedoms         : Integer
  !*                           : Number of fixed displacements
  !*
  !*    limit                  : Integer
  !*                           : Maximum number of PCG iterations allowed
  !*
  !*    loaded_nodes           : Integer
  !*                           : Number of nodes with applied forces
  !*
  !*    mesh                   : Integer
  !*                           : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*
  !*    nels                   : Integer
  !*                           : Total number of elements
  !*
  !*    nip                    : Integer
  !*                           : Number of Gauss integration points
  !*
  !*    nn                     : Integer
  !*                           : Total number of nodes in the mesh
  !*
  !*    nod                    : Integer
  !*                           : Number of nodes per element
  !*
  !*    nr                     : Integer
  !*                           : Number of nodes with restrained degrees of
  !*                             freedom 
  !*
  !*    tol                    : Real
  !*                           : Tolerance for PCG
  !*
  !*    v                      : Real
  !*                           : Poisson coefficient
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    03.03.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: limit,mesh,fixed_freedoms 
  REAL(iwp), INTENT(INOUT)         :: e,v,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
  READ(10,*) element,mesh,nels,nn,nr,nip,nod,loaded_nodes,fixed_freedoms,  &
             e,v,tol,limit
  CLOSE(10)
   
  IF(fixed_freedoms .AND. loaded_nodes > 0) THEN
    PRINT *
    PRINT *, "Error - model has", fixed_freedoms, " fixed freedoms and"
    PRINT *, loaded_nodes, " loaded nodes"
    PRINT *, "Mixed displacement and load control not supported"
    PRINT *
    CALL shutdown()
  END IF

  RETURN
  END SUBROUTINE READ_P121

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END MODULE INPUT
