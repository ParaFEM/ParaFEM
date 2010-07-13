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
  !*    READ_NUM               Reads the element nodal steering array
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
  !*        Available under commercial licence
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
    
  SUBROUTINE READ_G_NUM_PP(job_name,iel_start,nels,nn,numpe,g_num_pp)

  !/****f* input/read_g_num_pp
  !*  NAME
  !*    SUBROUTINE: read_g_num_pp
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_num_pp(job_name,iel_start,nels,nn,numpe,      &
  !*                                   g_num_pp)
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
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  The method based on building global arrays.
  !*  
  !*  An improvement would be to avoid allocating the global array
  !*
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
    DO i = 1,ubound(rest,1)  !ubound(rest,1) = nr
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

  SUBROUTINE READ_P121(job_name,numpe,e,element,limit,loaded_nodes,mesh,nels, &
                       nip,nn,nod,nr,tol,v)

  !/****f* input/read_p121
  !*  NAME
  !*    SUBROUTINE: read_p121
  !*  SYNOPSIS
  !*    Usage:      CALL read_p121(job_name,numpe,e,element,limit,             &
  !*                               loaded_nodesmesh,nels,nip,nn,nod,nr,tol,v
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
  INTEGER, INTENT(INOUT)           :: limit,mesh 
  REAL(iwp), INTENT(INOUT)         :: e,v,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(8)
  REAL(iwp)                        :: real_store(3)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,nels,nn,nr,nip,nod,loaded_nodes,e,v,              &
               tol,limit
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nr 
    integer_store(5)   = nip
    integer_store(6)   = nod
    integer_store(7)   = loaded_nodes
    integer_store(8)   = limit

    real_store         = 0.0_iwp

    real_store(1)      = e  
    real_store(2)      = v  
    real_store(3)      = tol  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 8
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 3
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
    limit        = integer_store(8)

    e            = real_store(1)
    v            = real_store(2)
    tol          = real_store(3)

  END IF

  RETURN
  END SUBROUTINE READ_P121

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P126(job_name,numpe,cjits,cjtol,ell,kappa,limit,meshgen,    &
                       nels,nip,nn,nr,nres,penalty,rho,tol,x0,visc)

  !/****f* input/read_p126
  !*  NAME
  !*    SUBROUTINE: read_p126
  !*  SYNOPSIS
  !*    Usage:      CALL read_p126(job_name,numpe,cjits,cjtol,ell,kappa,      &
  !*                               limit,meshgen,nels,nip,nn,nr,nres,penalty, &
  !*                               rho,tol,x0,visc
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
  !*                           : BiCGSTAB(l) parameter
  !*
  !*    cjtol                  : Real
  !*                           : BiCGSTAB(l) parameter
  !*
  !*    ell                    : Integer
  !*                           : BiCGSTAB(l) parameter
  !*
  !*    limit                  : Integer
  !*                           : Maximum number of iterations allowed
  !*
  !*    meshgem                : Integer
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
  !*    nres                   : Integer
  !*                           :
  !*
  !*    penalty                : Real
  !*                           :
  !*
  !*    rho                    : Real
  !*                           : Fluid parameter
  !*
  !*    tol                    : Real
  !*                           : Tolerance for PCG
  !*
  !*    x0                     : Real
  !*                           :
  !*
  !*    visc                   : Real
  !*                           : Fluid viscosity
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
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nres,nip,cjits
  INTEGER, INTENT(INOUT)           :: limit,meshgen,ell 
  REAL(iwp), INTENT(INOUT)         :: cjtol,kappa,penalty,rho,tol,x0,visc

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(9)
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
    READ(10,*) meshgen,nels,nn,nr,nres,nip,visc,rho,tol,limit,                &
               cjtol,cjits,penalty,x0,ell,kappa
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = meshgen
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nr 
    integer_store(5)   = nres 
    integer_store(6)   = nip
    integer_store(7)   = limit
    integer_store(8)   = cjits
    integer_store(9)   = ell

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

  bufsize = 9 
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 7
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)


!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    meshgen      = integer_store(1)
    nels         = integer_store(2)
    nn           = integer_store(3)
    nr           = integer_store(4)
    nres         = integer_store(5)
    nip          = integer_store(6)
    limit        = integer_store(7)
    cjits        = integer_store(8)
    ell          = integer_store(8)

    visc         = real_store(1)
    rho          = real_store(2)
    tol          = real_store(3)
    cjtol        = real_store(4)
    penalty      = real_store(5)
    x0           = real_store(6)
    kappa        = real_store(7)

  END IF

  RETURN
  END SUBROUTINE READ_P126

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P128ar(job_name,numpe,bmat,e,element,maxitr,mesh,ncv,nels,  &
                         nev,nip,nn,nod,nr,rho,tol,v,which)

  !/****f* input/read_p128ar
  !*  NAME
  !*    SUBROUTINE: read_p128ar
  !*  SYNOPSIS
  !*    Usage:      CALL read_p128ar(job_name,numpe,bmat,e,element,maxitr,    &
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
  !*    (c) University of Manchester 2010
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
  END SUBROUTINE READ_P128ar

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P129(job_name,numpe,alpha1,beta1,e,element,                 &
                       limit,loaded_nodes,mesh,nels,nip,nn,nod,               &
                       npri,nr,nstep,omega,rho,theta,tol,v)

  !/****f* input/read_p129
  !*  NAME
  !*    SUBROUTINE: read_p129
  !*  SYNOPSIS
  !*    Usage:      CALL read_p129(job_name,numpe,alpha1,beta1,e,element,     &
  !*                               limit,loaded_nodes,mesh,nels,nip,nn,nod,   &
  !*                               npri,nr,nstep,omega,rho,theta,tol,v)
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
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: nstep,npri,limit,mesh 
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,alpha1,beta1,theta,omega,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(10)
  REAL(iwp)                        :: real_store(8)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,nels,nn,nr,nip,nod,loaded_nodes,rho,e,v,          &
               alpha1,beta1,nstep,npri,theta,omega,tol,limit
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

  bufsize = 10
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
                        nels,nip,nn,nod,npri,nr,nstep,pload,rho,sbary,v)

  !/****f* input/read_p1210
  !*  NAME
  !*    SUBROUTINE: read_p1210
  !*  SYNOPSIS
  !*    Usage:      CALL read_p1210(job_name,numpe,dtim,e,element,            &
  !*                                loaded_nodes,meshgen,nels,nip,nn,nod,npri,&
  !*                                nr,nstep,pload,rho,sbary,v)
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
  INTEGER, INTENT(INOUT)           :: loaded_nodes,meshgen 
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,sbary,dtim,pload

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(9)
  REAL(iwp)                        :: real_store(6)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,meshgen,nels,nip,nn,nr,nod,loaded_nodes,rho,e,v,       &
               sbary,pload,dtim,nstep,npri
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

  bufsize = 9
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

END MODULE INPUT
