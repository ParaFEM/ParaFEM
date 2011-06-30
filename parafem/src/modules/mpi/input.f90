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
  !*    READ_MATERIALID_PP     Reads the material ID for each element
  !*    READ_MATERIALVALUE     Reads property values for each material ID
  !*    READ_NELS_PP           Reads number of elements assigned to processor
  !*    READ_P121              Reads the control data for program p121
  !*    READ_P126              Reads the control data for program p126
  !*    READ_P129              Reads the control data for program p129
  !*    READ_P1212             Reads the control data for program p1212
  !*    READ_XX2               Reads the control data for program ed5/xx2
  !*    BCAST_INPUTDATA_XX4    Reads the control data for program ed5/xx4
  !*    CHECK_INPUTDATA_XX4    Checks the control data for program ed5/xx4
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2004-2011 University of Manchester
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
  INTEGER                      :: nnStart         ! first node ID in g_coord
  INTEGER                      :: nnEnd           ! last node ID in g_coord
  INTEGER                      :: bufsize         ! packet size for bcast data
  INTEGER                      :: ier             ! MPI error code
  INTEGER                      :: iel, i, j, k, l ! loop counters
  INTEGER                      :: bitBucket
  INTEGER                      :: readSteps
  INTEGER                      :: readCount
  INTEGER                      :: readRemainder    
  REAL(iwp), INTENT(INOUT)     :: g_coord_pp(:,:,:) 
  REAL(iwp), ALLOCATABLE       :: g_coord(:,:)    ! temporary array
  REAL(iwp)                    :: zero = 0.0_iwp

!------------------------------------------------------------------------------
! 1. Find READSTEPS, the number of steps in which the read will be carried
!    out, READCOUNT, the size of each read and READREMAINDER, the number
!    of entries to be read after the last READSTEP.
!------------------------------------------------------------------------------

  fname     = job_name(1:INDEX(job_name, " ")-1) // ".d"  

  IF(npes > nn) THEN
    readSteps     = 1
    readRemainder = 0
    readCount     = nn
  ELSE
    readSteps     = npes
    readRemainder = MOD(nn,readSteps)
    readCount     = (nn - readRemainder) / npes
  END IF
  
!------------------------------------------------------------------------------
! 2. Allocate temporary array
!------------------------------------------------------------------------------

  nod      = UBOUND(g_coord_pp,1)
  ndim     = UBOUND(g_coord_pp,2)
  nels_pp  = UBOUND(g_coord_pp,3)

  ALLOCATE(g_coord(ndim,readCount))

!------------------------------------------------------------------------------
! 3. Master process opens the data file
!------------------------------------------------------------------------------

  IF(numpe==1)THEN
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*)   !headers
    READ(10,*)   !headers
  END IF

!------------------------------------------------------------------------------
! 4. Go round READSTEPS loop, read data, broadcast and populate g_coord_pp
!------------------------------------------------------------------------------

  g_coord_pp= zero
  bufsize   = ndim * readCount
  
  DO i = 1, readSteps
    g_coord = zero
    nnStart = (i-1) * readCount + 1
    nnEnd   =  i    * readCount 
    IF(numpe == 1) THEN
      DO j = 1,readCount
        READ(10,*) bitBucket,g_coord(:,j)
      END DO
    END IF
    CALL MPI_BCAST(g_coord,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
    DO iel = 1, nels_pp
      DO k = 1, nod
        IF(g_num_pp(k,iel) < nnStart) CYCLE
        IF(g_num_pp(k,iel) > nnEnd)   CYCLE
        l                  = g_num_pp(k,iel) - nnStart + 1
        g_coord_pp(k,:,iel)= g_coord(:,l)
      END DO
    END DO
  END DO
  
!------------------------------------------------------------------------------
! 5. If READREMAINDER > 0, collect remaining entries
!------------------------------------------------------------------------------
  
  IF(readRemainder > 0) THEN
    DEALLOCATE(g_coord)
    ALLOCATE(g_coord(ndim,readRemainder))
    g_coord  = zero
    bufsize  = ndim * readRemainder
    nnStart  = (readSteps * readCount) + 1
    nnEnd    = nnStart + readRemainder - 1
    IF(nnEnd > nn) THEN
      PRINT *, "Too many nodes"
      CALL shutdown()
      STOP
    END IF
    IF(numpe == 1) THEN
      DO j = 1,readRemainder
        READ(10,*) bitBucket,g_coord(:,j)
      END DO
    END IF
    CALL MPI_BCAST(g_coord,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
    DO iel = 1, nels_pp
      DO k = 1, nod
        IF(g_num_pp(k,iel) < nnStart) CYCLE
        IF(g_num_pp(k,iel) > nnEnd)   CYCLE
        l                  = g_num_pp(k,iel) - nnStart + 1
        g_coord_pp(k,:,iel)= g_coord(:,l)
      END DO
    END DO
  END IF
  
!------------------------------------------------------------------------------
! 6. Deallocate global arrays and close data file
!------------------------------------------------------------------------------

  DEALLOCATE(g_coord)    
  IF(numpe==1) CLOSE(10)

  RETURN
  
  END SUBROUTINE READ_G_COORD_PP
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_G_NUM_PP2(job_name,iel_start,nn,npes,numpe,g_num_pp)

  !/****f* input/read_g_num_pp2
  !*  NAME
  !*    SUBROUTINE: read_g_num_pp2
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,     &
  !*                                    g_num_pp)
  !*  FUNCTION
  !*    Master process reads the global array of elements and broadcasts
  !*    to slave processes.
  !*    Processes record only its local part of elements.
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
  !*    (c) University of Manchester 2007-2011
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  This method avoids the allocation of a global array.
  !*  Delete READ_G_NUM_PP and replace it with this subroutine once its 
  !*  stability has been demonstrated.
  !*
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname
  INTEGER, INTENT(IN)           :: iel_start, nn, npes, numpe
  INTEGER, INTENT(INOUT)        :: g_num_pp(:,:)
  INTEGER                       :: nod, nels_pp, iel, i, k
  INTEGER                       :: bufsize, ielpe, ier, ndim=3
  INTEGER                       :: readSteps,max_nels_pp
  INTEGER                       :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: g_num(:,:),localCount(:),readCount(:)

!------------------------------------------------------------------------------
! 1. Initiallize variables
!------------------------------------------------------------------------------

  nod       = UBOUND(g_num_pp,1)
  nels_pp   = UBOUND(g_num_pp,2)

!------------------------------------------------------------------------------
! 2. Find READSTEPS, the number of steps in which the read will be carried
!    out and READCOUNT, the size of each read.
!------------------------------------------------------------------------------

  ALLOCATE(readCount(npes))
  ALLOCATE(localCount(npes))
  
  readCount         = 0
  localCount        = 0
  readSteps         = npes
  localCount(numpe) = nels_pp

  CALL MPI_ALLREDUCE(localCount,readCount,npes,MPI_INTEGER,MPI_SUM,           &
                     MPI_COMM_WORLD,ier) 
 
!------------------------------------------------------------------------------
! 3. Allocate the array g_num into which the steering array for each processor
!    is to be read
!------------------------------------------------------------------------------

  max_nels_pp = MAXVAL(readCount,1)
  
  ALLOCATE(g_num(nod,max_nels_pp))
  
  g_num = 0             ! different value for each processor

!------------------------------------------------------------------------------
! 4. Master processor opens the data file and advances to the start of the 
!    element steering array
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname     = job_name(1:INDEX(job_name, " ")-1) // ".d"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*)   !header
    READ(10,*)   !header
    DO i = 1,nn  
      READ(10,*) !skip nodes until reaching the elements
    END DO
    READ(10,*)   !header
  END IF

!------------------------------------------------------------------------------
! 5. Go around READSTEPS loop, read data, and send to appropriate processor
!------------------------------------------------------------------------------

  DO i=1,npes
    IF(i == 1) THEN  ! local data
      IF(numpe == 1) THEN
        DO iel = 1,readCount(i)
          READ(10,*)k,k,k,k,g_num_pp(:,iel),k
        END DO
      END IF
    ELSE
      bufsize = readCount(i)*nod
      IF(numpe == 1) THEN
        DO iel = 1,readCount(i)
          READ(10,*)k,k,k,k,g_num(:,iel),k
        END DO
        CALL MPI_SEND(g_num(:,1:readCount(i)),bufsize,MPI_INTEGER,i-1,i,     &
                      MPI_COMM_WORLD,status,ier)
      END IF
      IF(numpe == i) CALL MPI_RECV(g_num_pp,bufsize,MPI_INTEGER,0,i,         &
                                   MPI_COMM_WORLD,status,ier)
    END IF
  END DO
  
!------------------------------------------------------------------------------
! 6. Close file and deallocate global arrays
!------------------------------------------------------------------------------

  IF(numpe==1) CLOSE(10)

  DEALLOCATE(g_num,readCount,localCount)
 
  RETURN

  END SUBROUTINE READ_G_NUM_PP2
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_ELEMENTS(job_name,iel_start,nn,npes,numpe,etype_pp,g_num_pp)

  !/****f* input/read_elements
  !*  NAME
  !*    SUBROUTINE: read_elements
  !*  SYNOPSIS
  !*    Usage:      CALL read_elements(job_name,iel_start,nn,npes,numpe,      &
  !*                                   etype_pp,g_num_pp)
  !*  FUNCTION
  !*    Master process reads the global array of elements and broadcasts
  !*    to slave processes.
  !*    Processes record only its local part of elements.
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
  !*    The following scalar array arguments have the INTENT(OUT) attribute:
  !*
  !*    g_num_pp(nod,nels_pp) : Element connectivity
  !*    etype_pp(nels_pp)     : Element property type vector     
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2011
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  This method avoids the allocation of a global array.
  !*  Delete READ_G_NUM_PP and replace it with this subroutine once its 
  !*  stability has been demonstrated.
  !*
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname
  INTEGER, INTENT(IN)           :: iel_start, nn, npes, numpe
  INTEGER, INTENT(INOUT)        :: g_num_pp(:,:),etype_pp(:)
  INTEGER                       :: nod, nels_pp, iel, i, k
  INTEGER                       :: bufsize, ielpe, ier, ndim=3
  INTEGER                       :: readSteps,max_nels_pp
  INTEGER                       :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: g_num(:,:),localCount(:),readCount(:)

!------------------------------------------------------------------------------
! 1. Initiallize variables
!------------------------------------------------------------------------------

  nod       = UBOUND(g_num_pp,1)
  nels_pp   = UBOUND(g_num_pp,2)

!------------------------------------------------------------------------------
! 2. Find READSTEPS, the number of steps in which the read will be carried
!    out and READCOUNT, the size of each read.
!------------------------------------------------------------------------------

  ALLOCATE(readCount(npes))
  ALLOCATE(localCount(npes))
  
  readCount         = 0
  localCount        = 0
  readSteps         = npes
  localCount(numpe) = nels_pp

  CALL MPI_ALLREDUCE(localCount,readCount,npes,MPI_INTEGER,MPI_SUM,           &
                     MPI_COMM_WORLD,ier) 
 
!------------------------------------------------------------------------------
! 3. Allocate the array g_num into which the steering array for each processor
!    is to be read
!------------------------------------------------------------------------------

  max_nels_pp = MAXVAL(readCount,1)
  
  ALLOCATE(g_num(nod+1,max_nels_pp))
  
  g_num = 0             ! different value for each processor

!------------------------------------------------------------------------------
! 4. Master processor opens the data file and advances to the start of the 
!    element steering array
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname     = job_name(1:INDEX(job_name, " ")-1) // ".d"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*)   !header
    READ(10,*)   !header
    DO i = 1,nn  
      READ(10,*) !skip nodes until reaching the elements
    END DO
    READ(10,*)   !header
  END IF

!------------------------------------------------------------------------------
! 5. Go around READSTEPS loop, read data, and send to appropriate processor
!------------------------------------------------------------------------------

  DO i=1,npes
    IF(i == 1) THEN  ! local data
      IF(numpe == 1) THEN
        DO iel = 1,readCount(i)
          READ(10,*)k,k,k,k,g_num_pp(:,iel),etype_pp(iel)
        END DO
      END IF
    ELSE
      bufsize = readCount(i)*(nod+1)
      IF(numpe == 1) THEN
        DO iel = 1,readCount(i)
          READ(10,*)k,k,k,k,g_num(:,iel)
        END DO
        CALL MPI_SEND(g_num(:,1:readCount(i)),bufsize,MPI_INTEGER,i-1,i,     &
                      MPI_COMM_WORLD,status,ier)
      END IF
      IF(numpe == i) THEN
        g_num = 0
        CALL MPI_RECV(g_num,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,status,ier)
        g_num_pp = g_num(1:nod,:)
        etype_pp = g_num(nod+1,:)
      END IF
    END IF
  END DO
  
!------------------------------------------------------------------------------
! 6. Close file and deallocate global arrays
!------------------------------------------------------------------------------

  IF(numpe==1) CLOSE(10)

  DEALLOCATE(g_num,readCount,localCount)

  RETURN

  END SUBROUTINE READ_ELEMENTS
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_G_NUM_PP(job_name,iel_start,nels,nn,numpe,g_num_pp)

  !/****f* input/read_g_num_pp
  !*  NAME
  !*    SUBROUTINE: read_g_num_pp
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_num_pp(job_name,iel_start,nels,nn,numpe,
  !*                                   g_num_pp) 
  !*  FUNCTION
  !*    The master processor reads sections of the global steering array and 
  !*    sends these to the slave processors. A loop over number of processors
  !*    NPES means that each slave processor is dealt with sequentially and
  !*    a global array of size NELS is not required.  
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    job_name              : Character
  !*                          : Used to create file name to read
  !*
  !*    iel_start             : Integer
  !*                          : First element number in a process
  !*
  !*    nels                  : Integer
  !*                          : Total number of elements
  !*
  !*    nn                    : Integer
  !*                          : Total number of nodes 
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
  !*  The method based on building global arrays.
  !*  
  !*  An improvement would be to avoid allocating the global array
  !*
  !*  Can copy the method in subroutine READ_G_COORD_PP
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: header
  CHARACTER(LEN=50)             :: fname
  INTEGER, INTENT(IN)           :: iel_start, nels, nn, numpe
  INTEGER, INTENT(INOUT)        :: g_num_pp(:,:)
  INTEGER                       :: nod, nels_pp, iel, i, k
  INTEGER                       :: bufsize, ielpe, ier, ndim=3
  INTEGER, ALLOCATABLE          :: g_num(:,:)
  REAL(iwp), ALLOCATABLE        :: dummy(:)

!------------------------------------------------------------------------------
! 1. Allocate global arrays
!------------------------------------------------------------------------------
  
  nod = UBOUND(g_num_pp,1)

  ALLOCATE(g_num(nod,nels))
  ALLOCATE(dummy(ndim))
  g_num = 0
  dummy = 0.0_iwp

!------------------------------------------------------------------------------
! 2. Master processor reads the data
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname     = job_name(1:INDEX(job_name, " ")-1) // ".d"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*)   header  !header
    READ(10,*)   header  !header
    DO i = 1,nn   !read nodes until reaching the elements
      READ(10,*) k, dummy(:)
    END DO
    READ(10,*)   header !keyword element
    DO iel = 1,nels
      READ(10,*)k,k,k,k,g_num(:,iel),k
    END DO
    CLOSE(10)
  END IF

!------------------------------------------------------------------------------
! 3. Master process broadcasts the data to slave processes
!------------------------------------------------------------------------------
 
  bufsize = nod*nels
  CALL MPI_BCAST(g_num,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Each process records only its local elements
!------------------------------------------------------------------------------
 
  nels_pp = UBOUND(g_num_pp,2)
  ielpe   = iel_start

  DO iel = 1,nels_pp
    g_num_pp(:,iel) = g_num(:,ielpe) 
    ielpe           = ielpe + 1
  END DO
 
!------------------------------------------------------------------------------
! 5. Deallocate global arrays
!------------------------------------------------------------------------------

  DEALLOCATE(g_num)
 
  RETURN

  END SUBROUTINE READ_G_NUM_PP
                
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

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
  INTEGER                       :: i,bufsize1,bufsize2,ier,ielpe
  INTEGER,INTENT(IN)            :: numpe
  INTEGER,INTENT(INOUT)         :: node(:)
  REAL(iwp),INTENT(INOUT)       :: value(:,:)

  IF(numpe==1)THEN
    fname = job_name(1:INDEX(job_name, " ")-1) // ".lds"
    OPEN(22, FILE=fname, STATUS='OLD', ACTION='READ')
    DO i = 1,ubound(node,1)   !ubound(node,1)=loaded_nodes 
      READ(22,*) node(i),value(:,i)
    END DO
    CLOSE(22)
  END IF

  bufsize1 = ubound(node,1)
  bufsize2 = ubound(value,1)*ubound(value,2)
  CALL MPI_BCAST(node,bufsize1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  CALL MPI_BCAST(value,bufsize2,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  RETURN

  END SUBROUTINE READ_LOADS

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE READ_LOADS_NS(JOB_NAME,NUMPE,NO,VAL)

  !/****f* input/read_loads_ns
  !*  NAME
  !*    SUBROUTINE: read_loads_ns
  !*  SYNOPSIS
  !*    Usage:      CALL read_loads_ns(job_name,numpe,no,val)
  !*  FUNCTION
  !*    Master processor reads the lid velocities for program p126 of Smith 
  !*    and Griffiths and sends a copy to all the slave processors. 
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
  !*    no(fixed_equations)       : Integer
  !*                              : Equations that have a prescribed velocity 
  !*
  !*    val(ndim,fixed_equations) : Real
  !*                              : Velocity that is applied to each equation
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
  INTEGER                       :: i,bufsize,ier,ielpe
  INTEGER,INTENT(IN)            :: numpe
  INTEGER,INTENT(INOUT)         :: no(:)
  REAL(iwp),INTENT(INOUT)       :: val(:)

  IF(numpe==1)THEN
    fname = job_name(1:INDEX(job_name, " ")-1) // ".lid"
    OPEN(22, FILE=fname, STATUS='OLD', ACTION='READ')
    DO i = 1,ubound(no,1)   !ubound(node,1)=fixed_equations 
      READ(22,*) no(i),val(i)
    END DO
    CLOSE(22)
  END IF

  bufsize = ubound(no,1)
  CALL MPI_BCAST(no,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  CALL MPI_BCAST(val,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  RETURN

  END SUBROUTINE READ_LOADS_NS

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
 
    !----------------------------------------------------------------------
    ! 1. Master process reads the data
    !----------------------------------------------------------------------

    fixed_freedoms = UBOUND(sense,1)

    IF(numpe==1)THEN
      fname = job_name(1:INDEX(job_name, " ")-1) // ".fix"
      OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
      DO i = 1,fixed_freedoms
        READ(10,*)node(i),sense(i),valf(i)
      END DO
      CLOSE(10)
    END IF

    !----------------------------------------------------------------------
    ! 2. Master process broadcasts the data to slave processes
    !----------------------------------------------------------------------
   
    CALL MPI_BCAST(node,fixed_freedoms,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(sense,fixed_freedoms,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(valf,fixed_freedoms,MPI_REAL8,0,MPI_COMM_WORLD,ier)

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
  !*  The method is based on reading a global array, although it's not a
  !*  worry because the nodes with constraints are a small proportion of the 
  !*  whole mesh.
  !*
  !*/   
  
  IMPLICIT NONE
  
  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname
  INTEGER                       :: i,bufsize1,bufsize2,ier,ielpe
  INTEGER,INTENT(IN)            :: numpe
  INTEGER,INTENT(INOUT)         :: rest(:,:)

  IF(numpe==1)THEN
    fname     = job_name(1:INDEX(job_name, " ")-1) // ".bnd"
    OPEN(23, FILE=fname, STATUS='OLD', ACTION='READ')
    DO i = 1, UBOUND(rest,1)
      READ(23,*) rest(i,:)
    END DO
    CLOSE(23)
  END IF

  bufsize1 = ubound(rest,1)*ubound(rest,2)

  CALL MPI_BCAST(rest,bufsize1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  RETURN
  END SUBROUTINE READ_REST

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE READ_MATERIALID_PP(MATERIALID_PP,FNAME,NN,NDIM,NELS,  &
                                NOD,IEL_START,NUMPE,NPES)
                                   
    !Master process reads the global array from a single file and 
    !sends a copy to each slave. (I/O method 1: refer to notes)
    !Each processor records only each correspondent part of the elements
    !The strategy is to send all the data to all the processors, and then
    !each processor records only what it is interested in

    IMPLICIT NONE
    CHARACTER(LEN=50), INTENT(in) :: fname
    INTEGER                  :: i,k,iel,bufsize,ier,ielpe
    INTEGER,INTENT(IN)       :: numpe,npes,iel_start,ndim,nn,nels,nod
    INTEGER,INTENT(INOUT)    :: materialID_pp(:)  !- local material IDs
    INTEGER,ALLOCATABLE      :: num(:)            !- node numbers
    INTEGER,ALLOCATABLE      :: materialID(:)     !- global material IDs
    REAL(iwp),ALLOCATABLE    :: coord(:)          !- coord
    CHARACTER(LEN=40)        :: keyword

    ALLOCATE(materialID(nels))
    ALLOCATE(coord(ndim))
    ALLOCATE(num(nod))

    IF(numpe==1)THEN

     OPEN(21,FILE=fname, STATUS='OLD', ACTION='READ')

     READ(21,*) keyword  !the headers of the file  (rubbish)
     READ(21,*) keyword  !the headers of the file  (rubbish)

!    the nodes are read, but nothing is done, just we wait until arriving
!    to the elements

     DO i = 1,nn 
      READ(21,*) k,coord(:)
     END DO

     READ(21,*) keyword  !here we read the word "elements"

     DO iel = 1,nels
       READ(21,*)k,k,k,k,num(:),materialID(iel)
     END DO

     CLOSE(21)
    END IF
    
    bufsize       = nels
    CALL MPI_BCAST(materialID,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    materialID_pp = 0 
    ielpe         = iel_start

    DO iel = 1,ubound(materialID_pp,1)
      materialID_pp(iel) = materialID(ielpe) 
      ielpe              = ielpe + 1
    END DO
   
    DEALLOCATE(materialID)
    DEALLOCATE(num)
    DEALLOCATE(coord)
   
    RETURN

  END SUBROUTINE READ_MATERIALID_PP
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_MATERIALVALUE(MATERIALVALUES,FNAME,NUMPE,NPES)

    !Master process reads the global array from a single file and 
    !sends a copy to each slave. (I/O method 1: refer to notes)
    !k is the material number in the file, it is read and discarded

    IMPLICIT NONE
    CHARACTER(LEN=50), INTENT(in) :: fname
    INTEGER                  :: i,k,nmats,bufsize,ier,ielpe
    INTEGER,INTENT(IN)       :: numpe,npes
    REAL(iwp),INTENT(INOUT)  :: materialValues(:,:)

    IF(numpe==1)THEN
     OPEN(21,FILE=fname, STATUS='OLD', ACTION='READ')
     READ(21,*) nmats
     DO i = 1,nmats
       READ(21,*)k, materialValues(:,i)
     END DO
     CLOSE(21)
    END IF
    
    bufsize       = ubound(materialValues,1)*ubound(materialValues,2)
    CALL MPI_BCAST(materialValues,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
   
    RETURN

  END SUBROUTINE READ_MATERIALVALUE

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_NELS_PP(job_name,nels_pp,npes,numpe)

  !/****f* input/read_nels_pp
  !*  NAME
  !*    SUBROUTINE: read_nels_pp
  !*  SYNOPSIS
  !*    Usage:      CALL read_nels_pp(job_name,nels_pp,npes,numpe)
  !*  FUNCTION
  !*    Reads NELS_PP from a file called job_name.psize 
  !*  INPUTS
  !*    The following scalar character argument has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer arguments have the INTENT(IN) attribute:
  !*
  !*    numpe                  : ID of calling processor
  !*    npes                   : Number of processors being used
  !*
  !*    The following scalar integer argument has the INTENT(OUT) attribute:
  !*
  !*    nels_pp                : Number of elements assigned to a processor
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    18.05.2011
  !*  COPYRIGHT
  !*    (c) University of Manchester 2011
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=50)                :: fname
  INTEGER, INTENT(IN)              :: npes,numpe
  INTEGER, INTENT(OUT)             :: nels_pp
  INTEGER, ALLOCATABLE             :: psize(:)
  INTEGER                          :: i,p

!------------------------------------------------------------------------------
! 1. Master processor reads the data into a temporary array 
!------------------------------------------------------------------------------

  ALLOCATE(psize(npes))
  psize = 0

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".psize"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) p,psize(:)
    IF(p/=npes) THEN
      PRINT*, "Number of partitions is different to the number of processes."
      CALL shutdown
    END IF
    CLOSE(10)
  END IF

!------------------------------------------------------------------------------
! 2. Array copied to all processors 
!------------------------------------------------------------------------------

  CALL MPI_BCAST(psize,npes,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 3. Each processor takes the value it needs
!------------------------------------------------------------------------------

  nels_pp = psize(numpe)
 
  DEALLOCATE(psize)

  RETURN
  END SUBROUTINE READ_NELS_PP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P121(job_name,numpe,e,element,fixed_freedoms,limit,        &
                       loaded_nodes,mesh,nels,nip,nn,nod,nr,partition,tol,v)

  !/****f* input/read_p121
  !*  NAME
  !*    SUBROUTINE: read_p121
  !*  SYNOPSIS
  !*    Usage:      CALL read_p121(job_name,numpe,e,element,fixed_freedoms,
  !*                               limit,loaded_nodes,mesh,nels,nip,nn,nod,nr,
  !*                               partition,tol,v)
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
  !*    partition              : Integer
  !*                           : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
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
  INTEGER, INTENT(INOUT)           :: limit,mesh,fixed_freedoms,partition 
  REAL(iwp), INTENT(INOUT)         :: e,v,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(10)
  REAL(iwp)                        :: real_store(3)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,nels,nn,nr,nip,nod,loaded_nodes,        &
               fixed_freedoms,e,v,tol,limit
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nr 
    integer_store(5)   = nip
    integer_store(6)   = nod
    integer_store(7)   = loaded_nodes
    integer_store(8)   = fixed_freedoms
    integer_store(9)   = limit
    integer_store(10)  = partition

    real_store         = 0.0_iwp

    real_store(1)      = e  
    real_store(2)      = v  
    real_store(3)      = tol  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 10
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 3
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  bufsize = 15
  CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    mesh            = integer_store(1)
    nels            = integer_store(2)
    nn              = integer_store(3)
    nr              = integer_store(4)
    nip             = integer_store(5)
    nod             = integer_store(6)
    loaded_nodes    = integer_store(7)
    fixed_freedoms  = integer_store(8)
    limit           = integer_store(9)
    partition       = integer_store(10)

    e               = real_store(1)
    v               = real_store(2)
    tol             = real_store(3)

  END IF

  IF(fixed_freedoms > 0 .AND. loaded_nodes > 0) THEN
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

  SUBROUTINE READ_P126(job_name,numpe,cjits,cjtol,ell,fixed_equations,kappa,  &
                       limit,mesh,nels,nip,nn,nr,partitioner,penalty,rho,tol, &
                       x0,visc)

  !/****f* input/read_p126
  !*  NAME
  !*    SUBROUTINE: read_p126
  !*  SYNOPSIS
  !*    Usage:      CALL read_p126(job_name,numpe,cjits,cjtol,ell,
  !*                               fixed_equations,kappa,limit,mesh,nels,     &
  !*                               nip,nn,nr,partitioner,penalty,rho,tol,x0,visc
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
  !*    cjits                  : Integer
  !*                           : Iteration limit for BiCGSTAB(l)
  !*
  !*    cjtol                  : Real
  !*                           : Tolerance for BiCGSTAB(l)
  !*
  !*    ell                    : Integer
  !*                           : l in BiCGSTAB(l) process
  !*
  !*    fixed_equations        : Integer
  !*                           : Number of fixed equations in the lid 
  !*
  !*    limit                  : Integer
  !*                           : Maximum number of iterations allowed
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
  !*    nr                     : Integer
  !*                           : Number of nodes with restrained degrees of
  !*                             freedom 
  !*
  !*    partitioner            : Integer
  !*                           : Type of partitioning
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*
  !*    penalty                : Real
  !*                           : Penalty
  !*
  !*    rho                    : Real
  !*                           : Fluid viscosity
  !*
  !*    tol                    : Real
  !*                           : Tolerance for outer loop
  !*
  !*    x0                     : Real
  !*                           : Starting value
  !*
  !*    visc                   : Real
  !*                           : Fluid viscosity
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    03.03.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010-11
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nip,cjits,fixed_equations
  INTEGER, INTENT(INOUT)           :: limit,mesh,ell,partitioner 
  REAL(iwp), INTENT(INOUT)         :: cjtol,kappa,penalty,rho,tol,x0,visc

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(10)
  REAL(iwp)                        :: real_store(7)
  CHARACTER(LEN=50)                :: fname
  CHARACTER(LEN=50)                :: program_name
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) program_name
    READ(10,*) mesh,partitioner,nels,nn,nr,fixed_equations,nip,visc,rho,tol,  &
               limit,cjtol,cjits,penalty,x0,ell,kappa
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nr 
    integer_store(5)   = fixed_equations 
    integer_store(6)   = nip
    integer_store(7)   = limit
    integer_store(8)   = cjits
    integer_store(9)   = ell
    integer_store(10)  = partitioner

    real_store         = 0.0_iwp

    real_store(1)      = visc  
    real_store(2)      = rho  
    real_store(3)      = tol
    real_store(4)      = cjtol
    real_store(5)      = penalty
    real_store(6)      = x0
    real_store(7)      = kappa
    
  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 10 
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 7
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)


!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    mesh            = integer_store(1)
    nels            = integer_store(2)
    nn              = integer_store(3)
    nr              = integer_store(4)
    fixed_equations = integer_store(5)
    nip             = integer_store(6)
    limit           = integer_store(7)
    cjits           = integer_store(8)
    ell             = integer_store(9)
    partitioner     = integer_store(10)

    visc            = real_store(1)
    rho             = real_store(2)
    tol             = real_store(3)
    cjtol           = real_store(4)
    penalty         = real_store(5)
    x0              = real_store(6)
    kappa           = real_store(7)
 
  END IF

  RETURN
  END SUBROUTINE READ_P126

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P1212(job_name,numpe,bmat,e,element,maxitr,mesh,ncv,nels,   &
                         nev,nip,nn,nod,nr,rho,tol,v,which)

  !/****f* input/read_p1212
  !*  NAME
  !*    SUBROUTINE: read_p1212
  !*  SYNOPSIS
  !*    Usage:      CALL read_p1212(job_name,numpe,bmat,e,element,maxitr,     &
  !*                                 mesh,ncv,nels,nev,nip,nn,nod,nr,rho,tol, &
  !*                                 v,which) 
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*    The following scalar character argument has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer argument has the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor ID of calling processor
  !*
  !*    The following scalar real arguments have the INTENT(INOUT) attribute:
  !*
  !*    e                      : Young's modulus
  !*    rho                    : Density
  !*    tol                    : Convergence tolerance
  !*    v                      : Poisson's ratio
  !*
  !*    The following scalar integer arguments have an INTENT(INOUT) attribute:
  !*
  !*    maxitr                 : Maximum number iterations allowed
  !*    mesh                   : Mesh numbering scheme
  !*                           : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    ncv                    : Arpack parameter
  !*    nels                   : Total number of elements
  !*    nev                    : Number of eigenvalues
  !*    nip                    : Number of integration points
  !*    nn                     : Number of nodes in the mesh
  !*    nod                    : Number of nodes in the element
  !*    nr                     : Number of restrained nodes
  !*
  !*    The following scalar character argument has an INTENT(INOUT) attribute:
  !*
  !*    bmat                   : Arpack parameter
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*    which                  : Arpack parameter
  !*  
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    05.07.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010-2011
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  CHARACTER(LEN=1), INTENT(INOUT)  :: bmat
  CHARACTER(LEN=2), INTENT(INOUT)  :: which
  CHARACTER(LEN=50)                :: fname
  CHARACTER(LEN=50)                :: program_name
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nod,nr,nip,nev,ncv
  INTEGER, INTENT(INOUT)           :: maxitr,mesh 
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(9)
  REAL(iwp)                        :: real_store(4)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) program_name
    READ(10,*) element,mesh,nels,nn,nr,nod,nip,rho,e,v,nev,ncv,bmat,which,   &
               tol,maxitr
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nod
    integer_store(5)   = nr 
    integer_store(6)   = nip
    integer_store(7)   = nev
    integer_store(8)   = ncv
    integer_store(9)   = maxitr

    real_store         = 0.0_iwp

    real_store(1)      = rho  
    real_store(2)      = e  
    real_store(3)      = v  
    real_store(4)      = tol  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 9 
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 4
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  bufsize = 15
  CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

  bufsize = 1
  CALL MPI_BCAST(bmat,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

  bufsize = 2
  CALL MPI_BCAST(which,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    mesh         = integer_store(1)
    nels         = integer_store(2)
    nn           = integer_store(3)
    nod          = integer_store(4)
    nr           = integer_store(5)
    nip          = integer_store(6)
    nev          = integer_store(7)
    ncv          = integer_store(8)
    maxitr       = integer_store(9)

    rho          = real_store(1)
    e            = real_store(2)
    v            = real_store(3)
    tol          = real_store(4)

  END IF

  RETURN
  END SUBROUTINE READ_P1212

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_XX2(job_name,numpe,element,fixed_freedoms,limit,            &
                      loaded_nodes,mesh,nels,nip,nn,nod,np_types,nr,          &
                      partition,tol)

  !/****f* input/read_xx2
  !*  NAME
  !*    SUBROUTINE: read_xx2
  !*  SYNOPSIS
  !*    Usage:      CALL read_xx2(job_name,numpe,element,fixed_freedoms,
  !*                              limit,loaded_nodes,mesh,nels,nip,nn,nod,
  !*                              np_types,nr,partition,tol)
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
  !*    partition              : Integer
  !*                           : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*
  !*    The following arguments have the INTENT(INOUT) attribute:
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
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    30.06.2011
  !*  COPYRIGHT
  !*    (c) University of Manchester 2011
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: limit,mesh,fixed_freedoms,partition 
  INTEGER, INTENT(INOUT)           :: np_types 
  REAL(iwp), INTENT(INOUT)         :: tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(11)
  REAL(iwp)                        :: real_store(1)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,np_types,nels,nn,nr,nip,nod,            &
               loaded_nodes,fixed_freedoms,tol,limit
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nr 
    integer_store(5)   = nip
    integer_store(6)   = nod
    integer_store(7)   = loaded_nodes
    integer_store(8)   = fixed_freedoms
    integer_store(9)   = limit
    integer_store(10)  = partition
    integer_store(11)  = np_types

    real_store         = 0.0_iwp

    real_store(1)      = tol  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 11
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 1
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  bufsize = 15
  CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    mesh            = integer_store(1)
    nels            = integer_store(2)
    nn              = integer_store(3)
    nr              = integer_store(4)
    nip             = integer_store(5)
    nod             = integer_store(6)
    loaded_nodes    = integer_store(7)
    fixed_freedoms  = integer_store(8)
    limit           = integer_store(9)
    partition       = integer_store(10)
    np_types        = integer_store(11)

    tol             = real_store(1)

  END IF

  IF(fixed_freedoms > 0 .AND. loaded_nodes > 0) THEN
    PRINT *
    PRINT *, "Error - model has", fixed_freedoms, " fixed freedoms and"
    PRINT *, loaded_nodes, " loaded nodes"
    PRINT *, "Mixed displacement and load control not supported"
    PRINT *
    CALL shutdown()
  END IF

  RETURN
  END SUBROUTINE READ_XX2
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE READ_P129(job_name,numpe,alpha1,beta1,e,element,                 &
                       limit,loaded_nodes,mesh,nels,nip,nn,nod,               &
                       npri,nr,nstep,omega,partitioner,rho,theta,tol,v)

  !/****f* input/read_p129
  !*  NAME
  !*    SUBROUTINE: read_p129
  !*  SYNOPSIS
  !*    Usage:      CALL read_p129(job_name,numpe,alpha1,beta1,e,element,     &
  !*                               limit,loaded_nodes,mesh,nels,nip,nn,nod,   &
  !*                               npri,nr,nstep,omega,partitioner,rho,theta, &
  !*                               tol,v)
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*    The following scalar character argument has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer argument has the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor ID of calling processor
  !*
  !*    The following scalar real arguments have the INTENT(INOUT) attribute:
  !*
  !*    alpha1                 : Rayleigh damping parameter
  !*    beta1                  : Rayleigh damping parameter
  !*    e                      : Young's modulus
  !*    omega                  : Intermediate value
  !*    rho                    : Density
  !*    theta                  : Parameter in "theta" integrator
  !*    tol                    : Convergence tolerance for PCG
  !*    v                      : Poisson's ratio
  !*
  !*    The following scalar integer arguments have an INTENT(INOUT) attribute:
  !*
  !*    limit                  : Maximum number of PCG iterations allowed
  !*    loaded_nodes           : Number of nodes with applied forces
  !*    mesh                   : Mesh numbering scheme
  !*                           : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    nels                   : Total number of elements
  !*    nip                    : Number of integration points
  !*    nn                     : Number of nodes in the mesh
  !*    nod                    : Number of nodes per element
  !*    npri                   : Print interval
  !*    nr                     : Number of restrained nodes
  !*    nstep                  : Number of time steps in analysis
  !*    partitioner            : Type of partitioning
  !*                           : 1 = Smith and Griffiths internal partitioning
  !*                           : 2 = External partitioning with .psize file
  !*
  !*    The following scalar character argument has an INTENT(INOUT) attribute:
  !*
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*  
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    25.02.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010-11
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  CHARACTER(LEN=50)                :: fname
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: nstep,npri,limit,mesh,partitioner 
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,alpha1,beta1,theta,omega,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(11)
  REAL(iwp)                        :: real_store(8)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partitioner,nels,nn,nr,nip,nod,loaded_nodes,     &
               rho,e,v,alpha1,beta1,nstep,npri,theta,omega,tol,limit
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nr 
    integer_store(5)   = nip
    integer_store(6)   = nod
    integer_store(7)   = loaded_nodes
    integer_store(8)   = nstep
    integer_store(9)   = npri
    integer_store(10)  = limit
    integer_store(11)  = partitioner

    real_store         = 0.0_iwp

    real_store(1)      = rho  
    real_store(2)      = e  
    real_store(3)      = v  
    real_store(4)      = alpha1  
    real_store(5)      = beta1  
    real_store(6)      = theta  
    real_store(7)      = omega  
    real_store(8)      = tol  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 11
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 8
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  bufsize = 15
  CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    mesh         = integer_store(1)
    nels         = integer_store(2)
    nn           = integer_store(3)
    nr           = integer_store(4)
    nip          = integer_store(5)
    nod          = integer_store(6)
    loaded_nodes = integer_store(7)
    nstep        = integer_store(8)
    npri         = integer_store(9)
    limit        = integer_store(10)
    partitioner  = integer_store(11)

    rho          = real_store(1)
    e            = real_store(2)
    v            = real_store(3)
    alpha1       = real_store(4)
    beta1        = real_store(5)
    theta        = real_store(6)
    omega        = real_store(7)
    tol          = real_store(8)

  END IF

  RETURN
  END SUBROUTINE READ_P129
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P1210(job_name,numpe,dtim,e,element,loaded_nodes,meshgen,   &
                        nels,nip,nn,nod,npri,nr,nstep,partitioner,pload,rho,  &
                        sbary,v)

  !/****f* input/read_p1210
  !*  NAME
  !*    SUBROUTINE: read_p1210
  !*  SYNOPSIS
  !*    Usage:      CALL read_p1210(job_name,numpe,dtim,e,element,            &
  !*                                loaded_nodes,meshgen,nels,nip,nn,nod,npri,&
  !*                                nr,nstep,partitioner,pload,rho,sbary,v)
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*    The following scalar character argument has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer arguments have the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor ID of calling processor
  !*    partitioner            : Type of partitioning used
  !*
  !*    The following scalar real arguments have the INTENT(INOUT) attribute:
  !*
  !*    dtim                   : Timestep
  !*    e                      : Young's modulus
  !*    pload                  : Load multiple
  !*    sbary                  : Shear yield stress
  !*    rho                    : Density
  !*    v                      : Poisson's ratio
  !*
  !*    The following scalar integer arguments have an INTENT(INOUT) attribute:
  !*
  !*    loaded_nodes           : Number of nodes with applied forces
  !*    meshgen                : Mesh numbering scheme
  !*                           : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    nels                   : Total number of elements
  !*    nip                    : Number of integration points
  !*    nn                     : Number of nodes in the mesh
  !*    nod                    : Number of nodes per element
  !*    npri                   : Print interval
  !*    nr                     : Number of restrained nodes
  !*    nstep                  : Number of time steps in analysis
  !*
  !*    The following scalar character argument has an INTENT(INOUT) attribute:
  !*
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*  
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    26.03.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  CHARACTER(LEN=50)                :: fname
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nip,nn,nod,npri,nr,nstep
  INTEGER, INTENT(INOUT)           :: loaded_nodes,meshgen,partitioner 
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,sbary,dtim,pload

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(10)
  REAL(iwp)                        :: real_store(6)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,meshgen,partitioner,nels,nip,nn,nr,nod,loaded_nodes,   &
               rho,e,v,sbary,pload,dtim,nstep,npri
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = loaded_nodes
    integer_store(2)   = meshgen
    integer_store(3)   = nels
    integer_store(4)   = nip
    integer_store(5)   = nn
    integer_store(6)   = nod
    integer_store(7)   = npri
    integer_store(8)   = nr 
    integer_store(9)   = nstep
    integer_store(10)  = partitioner

    real_store         = 0.0_iwp

    real_store(1)      = dtim
    real_store(2)      = e  
    real_store(3)      = pload
    real_store(4)      = rho  
    real_store(5)      = sbary
    real_store(6)      = v  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 10
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 6
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  bufsize = 15
  CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN
    
    loaded_nodes  = integer_store(1)
    meshgen       = integer_store(2)
    nels          = integer_store(3)
    nip           = integer_store(4)
    nn            = integer_store(5)
    nod           = integer_store(6)
    npri          = integer_store(7)
    nr            = integer_store(8)
    nstep         = integer_store(9)
    partitioner   = integer_store(10)

    dtim          = real_store(1)
    e             = real_store(2)
    pload         = real_store(3)
    rho           = real_store(4)
    sbary         = real_store(5)
    v             = real_store(6)
    
  END IF

  RETURN
  END SUBROUTINE READ_P1210
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    SUBROUTINE BCAST_INPUTDATA_XX4(numpe,npes,element,meshgen,nels,nn,nr,     &
                                     nip,plasitersMax,plasitersMin,           &
                                     loadIncrementsMax,                       &
                                     cjits,plastol,cjtol,fftol,ltol,          &
                                     numMaterials,numSteps)

    !Transfers input data from DAT file to all slave processors

    IMPLICIT NONE
    INTEGER,INTENT(INOUT)           :: numpe,npes,nels,nn,nr,nip
    INTEGER,INTENT(INOUT)           :: plasitersMax,plasitersMin,             &
                                       loadIncrementsMax,cjits,numMaterials,  &
                                       numSteps,meshgen
    INTEGER                         :: bufsizer,position,bufsize
    INTEGER                         :: bufdecl,recbufsize,ier
!   INTEGER, PARAMETER              :: ilength=4, rlength=8
    INTEGER, PARAMETER              :: ilength=8, rlength=8
    INTEGER, ALLOCATABLE            :: tempbuf(:)
    REAL(IWP),INTENT(INOUT)         :: plastol,cjtol,fftol,ltol
    CHARACTER(LEN=15),INTENT(INOUT) :: element

    bufsizer=11*ilength + 4*rlength

!   CALL MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(bufsizer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    bufdecl=bufsizer/4
    allocate(tempbuf(bufdecl))

!   IF(numpe==npes)THEN
    IF(numpe==1)THEN
      position = 0
      CALL MPI_PACK(meshgen,1,MPI_INTEGER,tempbuf,bufsizer,position,         &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(nels,1,MPI_INTEGER,tempbuf,bufsizer,position,            &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(nn,1,MPI_INTEGER,tempbuf,bufsizer,position,              &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(nr,1,MPI_INTEGER,tempbuf,bufsizer,position,              &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(nip,1,MPI_INTEGER,tempbuf,bufsizer,position,             &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(numMaterials,1,MPI_INTEGER,tempbuf,bufsizer,position,    &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(numSteps,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(plasitersMax,1,MPI_INTEGER,tempbuf,bufsizer,position,    &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(plasitersMin,1,MPI_INTEGER,tempbuf,bufsizer,position,    &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(loadIncrementsMax,1,MPI_INTEGER,tempbuf,bufsizer,        &
                    position,MPI_COMM_WORLD,ier)
      CALL MPI_PACK(cjits,1,MPI_INTEGER,tempbuf,bufsizer,position,           &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(plastol,1,MPI_REAL8,tempbuf,bufsizer,position,           &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(cjtol,1,MPI_REAL8,tempbuf,bufsizer,position,             &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(fftol,1,MPI_REAL8,tempbuf,bufsizer,position,             &
                    MPI_COMM_WORLD,ier)
      CALL MPI_PACK(ltol,1,MPI_REAL8,tempbuf,bufsizer,position,              &
                    MPI_COMM_WORLD,ier)
    END IF

!   CALL MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)
    CALL MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,0,MPI_COMM_WORLD,ier)

!   IF(numpe/=npes)THEN
    IF(numpe/=1)THEN
      position = 0
      CALL MPI_UNPACK(tempbuf,bufsizer,position,meshgen,1,MPI_INTEGER,       &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,nels,1,MPI_INTEGER,          &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,nn,1,MPI_INTEGER,            &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,nr,1,MPI_INTEGER,            &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,nip,1,MPI_INTEGER,           &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,numMaterials,1,MPI_INTEGER,  &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,numSteps,1,MPI_INTEGER,      &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,plasitersMax,1,MPI_INTEGER,  & 
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,plasitersMin,1,MPI_INTEGER,  & 
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,loadIncrementsMax,1,         &
                      MPI_INTEGER,MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,cjits,1,MPI_INTEGER,         &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,plastol,1,MPI_REAL8,         &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,cjtol,1,MPI_REAL8,           &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,fftol,1,MPI_REAL8,           &
                      MPI_COMM_WORLD,ier) 
      CALL MPI_UNPACK(tempbuf,bufsizer,position,ltol,1,MPI_REAL8,            &
                      MPI_COMM_WORLD,ier) 
    END IF

    bufsizer = 15
    CALL MPI_BCAST(element,bufsizer,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

    RETURN
    END SUBROUTINE BCAST_INPUTDATA_XX4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    SUBROUTINE CHECK_INPUTDATA_XX4(numpe,npes,element,meshgen,nels,nn,nr,   &
                                     nip,plasitersMax,plasitersMin,           &
                                     loadIncrementMax,                        &
                                     cjits,plastol,cjtol,fftol,ltol,          &
                                     numMaterials,numSteps)

    ! Checks input data from DAT file on master processor only

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)           :: numpe,npes,meshgen,nels,nn,nr,nip
    INTEGER,INTENT(IN)           :: plasitersMax,plasitersMin,                &
                                    loadIncrementMax,cjits,numMaterials,      &
                                    numSteps
    REAL(IWP),INTENT(IN)         :: plastol,cjtol,fftol,ltol
    CHARACTER(LEN=15),INTENT(IN) :: element

!------------------------------------------------------------------------------
! 1. I/O trap. Not enough NELS for the available processors
!------------------------------------------------------------------------------

    IF(nels < npes) THEN
      IF(numpe == 1) THEN
        WRITE(11,'(/2A)')   "------------------------------------------",     &
                            "----------"
        WRITE(11,'(A)')     "Program has detected a fatal error"
        WRITE(11,'(/A)')    "  The number of elements NELS must be greater"
        WRITE(11,'(A/)')    "  than the number of processors NPES"
        WRITE(11,'(A)')     "Analysis aborted"
        WRITE(11,'(2A/)')   "------------------------------------------",     &
                            "----------"
      END IF
      CALL shutdown()
    END IF

!------------------------------------------------------------------------------
! 2. Master processor echos the input data.
!------------------------------------------------------------------------------

    IF(numpe==1) THEN
            
      WRITE(11,'(/2A)')   "------------------------------------------",         &
                          "----------"
      WRITE(11,'(2A)')    "                     MODEL DATA           ",         &
                          "          "
      WRITE(11,'(2A/)')   "------------------------------------------",         &
                          "----------"
      WRITE(11,'(A,I10)') "Number of processors:                     ", npes
      WRITE(11,'(A,I10)') "Number of elements:                       ", nels
      WRITE(11,'(A,I10)') "Number of nodes:                          ", nn
      WRITE(11,'(A,I10)') "Number of restrained nodes:               ", nr
      WRITE(11,'(A,I10)') "Number of materials:                      ",         &
                           numMaterials
      WRITE(11,'(A,I10)') "Number of load steps:                     ", numSteps
      WRITE(11,'(A,I10)') "Upper limit on number of load increments: ",         &
                           loadIncrementMax 
      WRITE(11,'(A,I10)') "Lower limit of plastic iters per load inc:",         &
                           plasitersMin
      WRITE(11,'(A,I10)') "Upper limit of plastic iters per load inc:",         &
                           plasitersMax
    END IF 

!------------------------------------------------------------------------------
! 3. I/O trap for plasiters limits
!------------------------------------------------------------------------------
    
    IF(plasitersMin == plasitersMax .OR. plasitersMin > plasitersMax) THEN
      IF(numpe == 1) THEN
        WRITE(11,'(/2A)')   "------------------------------------------",       &
                            "----------"
        WRITE(11,'(A)')     "Program has detected a fatal error"
        WRITE(11,'(/A)')    "  Input variable PLASITERSMAX"
        WRITE(11,'(A/)')    "  must be greater than PLASITERSMIN"
        WRITE(11,'(A)')     "Analysis aborted"
        WRITE(11,'(2A/)')   "------------------------------------------",       &
                            "----------"
      END IF
      CALL shutdown()
    END IF

    RETURN
    END SUBROUTINE CHECK_INPUTDATA_XX4
    
END MODULE INPUT
