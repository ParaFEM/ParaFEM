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
  !*    GETNAME                Gets the base name of a ".dat" data file
  !*    GETNAME_MG             Gets the base name of an ".mg" data file
  !*    READ_G_COORD_PP        Reads the global coordinates
  !*    READ_G_COORD_PP_E   
  !*    READ_G_COORD_PP_BE     Reads the global coordinates (binary ensi)
  !*    READ_NODES
  !*    READ_G_NUM_PP          Reads the element nodal steering array
  !*    READ_G_NUM_PP_E 
  !*    READ_G_NUM_PP_BE       Reads element nodal steering array (binary ensi)
  !*    READ_ELEMENTS          Reads the element nodal steering array
  !*    READ_ELEMENTS_2        Reads the element nodal steering array
  !*    READ_ETYPE_PP          Reads the element material ID
  !*    READ_ETYPE_PP_BE       Reads the element material ID (binary ensi)
  !*    READ_G_NUM_PP2
  !*    READ_LOADS             Reads nodal forces
  !*    READ_LOADS_NS          Reads lid velocities for p126
  !*    READ_FIXED             Reads fixed freedoms for displacement control
  !*    READ_REST              Reads the restraints
  !*    READ_RESTRAINTS        Reads the restraints
  !*    READ_MATERIALID_PP     Reads the material ID for each element
  !*    READ_MATERIALVALUE     Reads property values for each material ID
  !*    READ_MATERIAL          Reads property values for each material ID
  !*    READ_NELS_PP           Reads number of elements assigned to processor
  !*    READ_LF
  !*    READ_P121              Reads the control data for program p121
  !*    READ_P122              Reads the control data for program p122
  !*    READ_P123              Reads the control data for program p123
  !*    READ_P124_4            Reads the control data for program p124, 4th ed
  !*    READ_P124              Reads the control data for program p124, 5th ed
  !*    READ_P125              Reads the control data for program p125
  !*    READ_P126              Reads the control data for program p126
  !*    READ_P127              Reads the control data for program p127
  !*    READ_P128
  !*    READ_P129              Reads the control data for program p129
  !*    READ_P1210
  !*    READ_XX1               Reads the control data for program xx1
  !*    READ_XX2               Reads the control data for program xx2
  !*    BCAST_INPUTDATA_XX5    Reads the control data for program xx5
  !*    CHECK_INPUTDATA_XX5    Checks the control data for program xx5
  !*    READ_XX6               Reads the control data for program xx6
  !*    READ_DATA_XX7          Reads the control data for program xx7
  !*    READ_XX12
  !*    READ_TIMESTEPS
  !*    READ_X_PP
  !*    READ_AMPLITUDE
  !*    BCAST_INPUTDATA_P123   Reads the control data for program p123
  !*    BCAST_INPUTDATA_P127   Reads the control data for program p127
  !*    READ_RFEMSOLVE
  !*    MESH_ENSI              Creates ASCII ensight gold files
  !*    MESH_ENSI_BIN          Creates BINARY ensight gold files *not tested*
  !*    MESH_ENSI_GEO_BIN      Creates BINARY ensight geo files *not tested*
  !*    MESH_ENSI_MATID_BIN    Creates BINARY ensight MATID files *not tested*
  !*    MESH_ENSI_NDBND_BIN    Creates BINARY ensight NDBND files *not tested*
  !*    MESH_ENSI_NDLDS_BIN    Creates BINARY ensight LDS files *not tested*
  !*    MESH_ENSI_CASE
  !*    MESH_ENSI_GEO
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2004-2014 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  If you contribute to this module, add your author name.
  !*
  !*/

  USE mpi_wrapper
  USE precision
  USE mp_interface

  CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

 SUBROUTINE getname(argv,nlen)
 
   !/****f* input/getname
   !*  NAME
   !*    SUBROUTINE: getname
   !*  SYNOPSIS
   !*    Usage:      CALL getname(argv,nlen)
   !*  FUNCTION
   !*    Returns the base name of the data file.
   !*  OUTPUTS
   !*    Scalar integer
   !*    nlen               : number of characters in data file base name 
   !*
   !*    Scalar character
   !*    argv               : holds data file base name
   !*  AUTHOR
   !*    I.M. Smith
   !*    D.V. Griffiths
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*  
   !*  Subroutine required by the 5th Edition of "Programming the Finite
   !*  Element Method". Take care when modifying
   !*/
  
    IMPLICIT NONE
    INTEGER                  :: narg
    INTEGER,INTENT(OUT)      :: nlen
    INTEGER                  :: lnblnk,iargc
    CHARACTER(*),INTENT(OUT) :: argv
    LOGICAL                  :: found=.false.
 
    narg=iargc()
    IF(narg.lt.1)THEN
      WRITE(*,*)'Please enter the base name of data file: '
      READ(*,*) argv
    ELSE
      CALL getarg(1,argv)
    ENDIF

   !nlen=lnblnk(argv)
    nlen=len_trim(argv)
    INQUIRE(file=argv(1:nlen)//'.dat',exist=found)
!   INQUIRE(file=argv(1:nlen)//'.mg',exist=found)   !hack

    IF(.not.found)THEN
      WRITE(*,*)'Data file not found'
      WRITE(*,*)'Please create or check spelling.'
      STOP
    ENDIF

  RETURN
  END SUBROUTINE getname

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE getname_mg(argv,nlen)
 
   !/****f* input/getname_mg
   !*  NAME
   !*    SUBROUTINE: getname_mg
   !*  SYNOPSIS
   !*    Usage:      CALL getname_mg(argv,nlen)
   !*  FUNCTION
   !*    Returns the base name of the data file.
   !*  OUTPUTS
   !*    Scalar integer
   !*    nlen               : number of characters in data file base name 
   !*
   !*    Scalar character
   !*    argv               : holds data file base name
   !*  AUTHOR
   !*    I.M. Smith
   !*    D.V. Griffiths
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*  
   !*  Duplicate of subroutine required by the 5th Edition of "Programming the Finite
   !*  Element Method". Needs to be merged with SUBROUTINE getname
   !*/
  
    IMPLICIT NONE
    INTEGER                  :: narg
    INTEGER,INTENT(OUT)      :: nlen
    INTEGER                  :: lnblnk,iargc
    CHARACTER(*),INTENT(OUT) :: argv
    LOGICAL                  :: found=.false.
 
    narg=iargc()
    IF(narg.lt.1)THEN
      WRITE(*,*)'Please enter the base name of data file: '
      READ(*,*) argv
    ELSE
      CALL getarg(1,argv)
    ENDIF

    nlen=len_trim(argv)
    INQUIRE(file=argv(1:nlen)//'.mg',exist=found) 

    IF(.not.found)THEN
      WRITE(*,*)'Data file not found'
      WRITE(*,*)'Please create or check spelling.'
      STOP
    ENDIF

  RETURN
  END SUBROUTINE getname_mg

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
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE READ_G_COORD_PP_E(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)

  !/****f* input/read_g_coord_pp_e
  !*  NAME
  !*    SUBROUTINE: read_g_coord_pp_e
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_coord_pp_e(job_name,g_num_pp,nn,npes,numpe,  &
  !*                                       g_coord_pp)
  !*  FUNCTION
  !*    
  !*  INPUTS
  !*
  !*  AUTHOR(s)
  !*    Lee Margetts
  !*    Llion Evans
  !*  CREATION DATE
  !*    22 September 2014
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2014
  !******
  !*  THIS SUBROUTINE IS WORK IN PROGRESS
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50),INTENT(IN) :: job_name
  CHARACTER(LEN=50)            :: fname
  CHARACTER(LEN=80)            :: cbuffer
  INTEGER, INTENT(IN)          :: nn, npes, numpe
  INTEGER, INTENT(IN)          :: g_num_pp(:,:)
  INTEGER                      :: ndim
  INTEGER                      :: nels_pp
  INTEGER                      :: nod
  INTEGER                      :: nnStart         ! first node ID in g_coord
  INTEGER                      :: nnEnd           ! last node ID in g_coord
  INTEGER                      :: bufsize         ! packet size for bcast data
  INTEGER                      :: ier             ! MPI error code
  INTEGER                      :: iel,i,j,k,l,m,n,p   ! loop counters
  INTEGER                      :: bitBucket
  INTEGER                      :: readSteps
  INTEGER                      :: readCount
  INTEGER                      :: readRemainder    
  INTEGER                      :: nn_in,part
  REAL(iwp), INTENT(INOUT)     :: g_coord_pp(:,:,:) 
  REAL(iwp), ALLOCATABLE       :: ord(:)    ! temporary array
  REAL(iwp)                    :: zero = 0.0_iwp
  LOGICAL                      :: verbose=.true.
! LOGICAL                      :: verbose=.false.

!------------------------------------------------------------------------------
! 1. Find READSTEPS, the number of steps in which the read will be carried
!    out, READCOUNT, the size of each read and READREMAINDER, the number
!    of entries to be read after the last READSTEP.
!------------------------------------------------------------------------------

  fname     = job_name(1:INDEX(job_name, " ")-1) // ".ensi.geo"  

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

  ALLOCATE(ord(readCount))

!------------------------------------------------------------------------------
! 3. Master process opens the data file
!------------------------------------------------------------------------------

  IF(numpe==1)THEN
    
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')

    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer

    READ(10,*)   part        ; IF(verbose) PRINT *, part

    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    
    READ(10,*)   nn_in       ; IF(verbose) PRINT *, "nn_in = ",nn_in
    
  END IF

!------------------------------------------------------------------------------
! 4. Go round READSTEPS loop, read data, broadcast and populate g_coord_pp
!------------------------------------------------------------------------------

! For Ensight Gold we need to read this as a loop on x, y and z ordinates
! separately

  DO m = 1, ndim ! Coordinates
  
    IF(ALLOCATED(ord)) DEALLOCATE(ord)
    ALLOCATE(ord(readCount))
    ord     = zero
    bufsize = readCount
    p       = 0

    DO i = 1, readSteps
      ord     = zero
      nnStart = (i-1) * readCount + 1
      nnEnd   =  i    * readCount 
      IF(numpe == 1) THEN
        READ(10,*) ord(1:readCount)
      END IF
      CALL MPI_BCAST(ord,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      DO iel = 1, nels_pp
        DO k = 1, nod
          IF(g_num_pp(k,iel) < nnStart) CYCLE
          IF(g_num_pp(k,iel) > nnEnd)   CYCLE
          l                  = g_num_pp(k,iel) - nnStart + 1
          g_coord_pp(k,m,iel)= ord(l)
        END DO
      END DO
    END DO
  
!------------------------------------------------------------------------------
! 5. If READREMAINDER > 0, collect remaining entries
!------------------------------------------------------------------------------
  
    IF(readRemainder > 0) THEN
      DEALLOCATE(ord)
      ALLOCATE(ord(readRemainder))
      ord      = zero
      bufsize  = readRemainder
      nnStart  = (readSteps * readCount) + 1
      nnEnd    = nnStart + readRemainder - 1
      IF(nnEnd > nn) THEN
        PRINT *, "Too many nodes"
        CALL shutdown()
        STOP
      END IF
      IF(numpe == 1) THEN
        READ(10,*) ord(1:readRemainder)
      END IF
      CALL MPI_BCAST(ord,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      DO iel = 1, nels_pp
        DO k = 1, nod
          IF(g_num_pp(k,iel) < nnStart) CYCLE
          IF(g_num_pp(k,iel) > nnEnd)   CYCLE
          l                  = g_num_pp(k,iel) - nnStart + 1
          g_coord_pp(k,m,iel)= ord(l)
        END DO
      END DO
    END IF

  END DO ! Coordinates
  
!------------------------------------------------------------------------------
! 6. Deallocate global arrays and close data file
!------------------------------------------------------------------------------

  DEALLOCATE(ord)    
  IF(numpe==1) CLOSE(10)

  RETURN
  
  END SUBROUTINE READ_G_COORD_PP_E

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE READ_G_COORD_PP_BE(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)

  !/****f* input/read_g_coord_pp_be
  !*  NAME
  !*    SUBROUTINE: read_g_coord_pp_be
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_coord_pp_be(job_name,g_num_pp,nn,npes,numpe,  &
  !*                                        g_coord_pp)
  !*  FUNCTION
  !*    
  !*  INPUTS
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    22 September 2014
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2014
  !******
  !*  THIS SUBROUTINE IS WORK IN PROGRESS
  !*/

  USE, INTRINSIC :: ISO_C_BINDING

  IMPLICIT NONE

  CHARACTER(LEN=50),INTENT(IN) :: job_name
  CHARACTER(LEN=50)            :: fname
  CHARACTER(LEN=80,KIND=C_CHAR)            :: cbuffer
  INTEGER, INTENT(IN)          :: nn, npes, numpe
  INTEGER, INTENT(IN)          :: g_num_pp(:,:)
  INTEGER                      :: ndim
  INTEGER                      :: nels_pp
  INTEGER                      :: nod
  INTEGER                      :: nnStart         ! first node ID in g_coord
  INTEGER                      :: nnEnd           ! last node ID in g_coord
  INTEGER                      :: bufsize         ! packet size for bcast data
  INTEGER                      :: ier             ! MPI error code
  INTEGER                      :: iel,i,j,k,l,m,n,p   ! loop counters
  INTEGER                      :: bitBucket
  INTEGER                      :: readSteps
  INTEGER                      :: readCount
  INTEGER                      :: readRemainder    
  INTEGER(KIND=C_INT)          :: nn_in,part
  REAL(iwp), INTENT(INOUT)     :: g_coord_pp(:,:,:) 
  REAL(KIND=C_FLOAT), ALLOCATABLE :: ord(:)    ! temporary array
  REAL(iwp)                    :: zero = 0.0_iwp
  LOGICAL                      :: verbose=.true.
! LOGICAL                      :: verbose=.false.

!------------------------------------------------------------------------------
! 1. Find READSTEPS, the number of steps in which the read will be carried
!    out, READCOUNT, the size of each read and READREMAINDER, the number
!    of entries to be read after the last READSTEP.
!------------------------------------------------------------------------------

  fname     = job_name(1:INDEX(job_name, " ")-1) // ".bin.ensi.geo"  

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

  ALLOCATE(ord(readCount))

!------------------------------------------------------------------------------
! 3. Master process opens the data file
!------------------------------------------------------------------------------

  IF(numpe==1)THEN
    
    OPEN(10,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',         &
                       ACCESS='STREAM')

    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer

    READ(10)   part    ; IF(verbose) PRINT *, part

    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    
    READ(10)   nn_in   ; IF(verbose) PRINT *, "nn_in = ",nn_in
    
  END IF

!------------------------------------------------------------------------------
! 4. Go round READSTEPS loop, read data, broadcast and populate g_coord_pp
!------------------------------------------------------------------------------

! For Ensight Gold we need to read this as a loop on x, y and z ordinates
! separately

  DO m = 1, ndim ! Coordinates
  
    IF(ALLOCATED(ord)) DEALLOCATE(ord)
    ALLOCATE(ord(readCount))
    ord     = zero
    bufsize = readCount
    p       = 0

    DO i = 1, readSteps
      ord     = zero
      nnStart = (i-1) * readCount + 1
      nnEnd   =  i    * readCount 
      IF(numpe == 1) THEN
        READ(10) ord(1:readCount)
      END IF
      CALL MPI_BCAST(ord,bufsize,MPI_REAL4,0,MPI_COMM_WORLD,ier)
      DO iel = 1, nels_pp
        DO k = 1, nod
          IF(g_num_pp(k,iel) < nnStart) CYCLE
          IF(g_num_pp(k,iel) > nnEnd)   CYCLE
          l                  = g_num_pp(k,iel) - nnStart + 1
          g_coord_pp(k,m,iel)= ord(l)
        END DO
      END DO
    END DO
  
!------------------------------------------------------------------------------
! 5. If READREMAINDER > 0, collect remaining entries
!------------------------------------------------------------------------------
  
    IF(readRemainder > 0) THEN
      DEALLOCATE(ord)
      ALLOCATE(ord(readRemainder))
      ord      = zero
      bufsize  = readRemainder
      nnStart  = (readSteps * readCount) + 1
      nnEnd    = nnStart + readRemainder - 1
      IF(nnEnd > nn) THEN
        PRINT *, "Too many nodes"
        CALL shutdown()
        STOP
      END IF
      IF(numpe == 1) THEN
        READ(10) ord(1:readRemainder)
      END IF
      CALL MPI_BCAST(ord,bufsize,MPI_REAL4,0,MPI_COMM_WORLD,ier)
      DO iel = 1, nels_pp
        DO k = 1, nod
          IF(g_num_pp(k,iel) < nnStart) CYCLE
          IF(g_num_pp(k,iel) > nnEnd)   CYCLE
          l                  = g_num_pp(k,iel) - nnStart + 1
          g_coord_pp(k,m,iel)= ord(l)
        END DO
      END DO
    END IF

  END DO ! Coordinates
  
!------------------------------------------------------------------------------
! 6. Deallocate global arrays and close data file
!------------------------------------------------------------------------------

  DEALLOCATE(ord)    
  IF(numpe==1) CLOSE(10)

  RETURN
  
  END SUBROUTINE READ_G_COORD_PP_BE

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE READ_NODES(fname,nn,nn_start,numpe,g_coord_pp)

    !/****f* input_output/read_nodes
    !*  NAME
    !*    SUBROUTINE: read_nodes
    !*  SYNOPSIS
    !*    Usage:      CALL read_nodes(fname,nn,nn_start,numpe,g_coord_pp)
    !*  FUNCTION
    !*    Master process reads the global array of nodal coordinates and
    !*    broadcasts to slave processes.
    !*    Processes record only its local part of nodal coordinates.
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    fname                  : Character
    !*                           : File name to read
    !*
    !*    nn                     : Integer
    !*                           : Total number of nodes 
    !*
    !*    nn_start               : Integer
    !*                           : First node number in a process
    !*
    !*    numpe                  : Integer
    !*                           : Process number
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    g_coord_pp(ndim,nn_pp) : Real
    !*                           : Nodal coordinates
    !*
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    01.06.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  The method based on building global arrays.
    !*  
    !*  An improvement would be to avoid allocating the global array
    !*  
    !*  THIS SUBROUTINE IS ONLY USED IN PROGRAM XX7 AND NEEDS TO BE REMOVED
    !*  FROM PARAFEM
    !*/
  
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: fname
    INTEGER,      INTENT(IN)  :: nn, nn_start, numpe
    REAL(iwp),    INTENT(OUT) :: g_coord_pp(:,:)
    INTEGER :: ndim, nn_pp, i, k, bufsize, inpe, ier
    REAL(iwp), ALLOCATABLE :: g_coord(:,:)

    !----------------------------------------------------------------------
    ! 1. Allocate global arrays
    !----------------------------------------------------------------------

    ndim  = UBOUND(g_coord_pp,1)

    ALLOCATE(g_coord(ndim,nn))

    !----------------------------------------------------------------------
    ! 2. Master process reads the data
    !----------------------------------------------------------------------

    IF(numpe==1)THEN
      OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
      READ(10,*)   !headers
      READ(10,*)   !headers
      DO i = 1,nn
        READ(10,*)k,g_coord(:,i)
      END DO
      CLOSE(10)
    END IF

    !----------------------------------------------------------------------
    ! 3. Master process broadcasts the data to slave processes
    !----------------------------------------------------------------------
   
    bufsize = ndim*nn

    CALL MPI_BCAST(g_coord,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

    !----------------------------------------------------------------------
    ! 4. Each process records only its local range of nodal coordinates
    !----------------------------------------------------------------------

    nn_pp = UBOUND(g_coord_pp,2) 
    inpe  = nn_start 

    DO i = 1,nn_pp
      g_coord_pp(:,i) = g_coord(:,inpe)
      inpe = inpe + 1
    END DO

    !----------------------------------------------------------------------
    ! 5. Deallocate global arrays
    !----------------------------------------------------------------------

    DEALLOCATE(g_coord)

    RETURN

  END SUBROUTINE READ_NODES
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_G_NUM_PP(job_name,iel_start,nn,npes,numpe,g_num_pp)

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

  END SUBROUTINE READ_G_NUM_PP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_G_NUM_PP_E(job_name,iel_start,nn,npes,numpe,g_num_pp)

  !/****f* input/read_g_num_pp_e
  !*  NAME
  !*    SUBROUTINE: read_g_num_pp_e
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_num_pp_e(job_name,iel_start,nn,npes,numpe,    &
  !*                                      g_num_pp)
  !*  FUNCTION
  !*    Master process reads the global array of elements from a ensight
  !*    gold "geo" file and broadcasts to slave processes.
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
  !*  AUTHOR(s)
  !*    L. Margetts
  !*    Ll. Evans
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2014
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  This subroutine is under development 22.09.2014
  !*
  !*/
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname
  CHARACTER(LEN=80)             :: cbuffer
  REAL                          :: rdummy
  INTEGER, INTENT(IN)           :: iel_start, nn, npes, numpe
  INTEGER, INTENT(INOUT)        :: g_num_pp(:,:)
  INTEGER                       :: nn_in,nels_in,part
  INTEGER                       :: nod, nels_pp, iel, i, j, k, ipos
  INTEGER                       :: bufsize, ielpe, ier, ndim=3
  INTEGER                       :: readSteps,max_nels_pp
  INTEGER                       :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: g_num(:,:),localCount(:),readCount(:)
! LOGICAL                       :: verbose=.false.
  LOGICAL                       :: verbose=.true.

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
    fname     = job_name(1:INDEX(job_name, " ")-1) // ".ensi.geo"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')

    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer

    READ(10,*)   part    ; IF(verbose) PRINT *, part

    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10,'(A)')   cbuffer ; IF(verbose) PRINT *, cbuffer
    
    READ(10,*)   nn_in   ; IF(verbose) PRINT *, "nn_in = ",nn_in
    
  END IF

!------------------------------------------------------------------------------
! 5. Go around READSTEPS loop, read data, and send to appropriate processor
!------------------------------------------------------------------------------

  DO i=1,npes
    IF(i == 1) THEN  ! local data
      IF(numpe == 1) THEN
	  
        DO j = 1,nn*ndim  
          READ(10,*) !skip nodes until reaching the elements
        END DO
		
        READ(10,'(A)') cbuffer     ; IF(verbose) PRINT *, cbuffer
        READ(10,*) nels_in         ; IF(verbose) PRINT *, nels_in
        iel = readCount(i)
        READ(10,*) g_num_pp(:,1:iel)
      END IF
    ELSE
      bufsize = readCount(i)*nod
      IF(numpe == 1) THEN
        iel = readCount(i)
        READ(10,*) g_num(:,1:iel)
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

  END SUBROUTINE READ_G_NUM_PP_E
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_G_NUM_PP_BE(job_name,iel_start,nn,npes,numpe,g_num_pp)

  !/****f* input/read_g_num_pp2_be
  !*  NAME
  !*    SUBROUTINE: read_g_num_pp2_be
  !*  SYNOPSIS
  !*    Usage:      CALL read_g_num_pp2_be(job_name,iel_start,nn,npes,numpe,  &
  !*                                    g_num_pp)
  !*  FUNCTION
  !*    Master process reads the global array of elements from a binary ensight
  !*    gold "geo" file and broadcasts to slave processes.
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
  !*    (c) University of Manchester 2007-2014
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  This subroutine is under development 22.09.2014
  !*
  !*/

  USE, INTRINSIC :: ISO_C_BINDING
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname
  CHARACTER(LEN=80)             :: cbuffer
  REAL(KIND=C_FLOAT)            :: rdummy
  INTEGER, INTENT(IN)           :: iel_start, nn, npes, numpe
  INTEGER, INTENT(INOUT)        :: g_num_pp(:,:)
  INTEGER(KIND=C_INT)           :: nn_in,nels_in,part
  INTEGER                       :: nod, nels_pp, iel, i, j, k, ipos
  INTEGER                       :: bufsize, ielpe, ier, ndim=3
  INTEGER                       :: readSteps,max_nels_pp
  INTEGER                       :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: g_num(:,:),localCount(:),readCount(:)
! LOGICAL                       :: verbose=.false.
  LOGICAL                       :: verbose=.true.

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
    fname     = job_name(1:INDEX(job_name, " ")-1) // ".bin.ensi.geo"
    OPEN(10,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',         &
                       ACCESS='STREAM')

    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer

    READ(10)   part    ; IF(verbose) PRINT *, part

    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    READ(10)   cbuffer ; IF(verbose) PRINT *, cbuffer
    
    READ(10)   nn_in   ; IF(verbose) PRINT *, "nn_in = ",nn_in
    
  END IF

!------------------------------------------------------------------------------
! 5. Go around READSTEPS loop, read data, and send to appropriate processor
!
!    Use INQUIRE to find out file position then skip coordinate data nn*ndim 
!    entries. Assumes single precision 4 byte reals for each value skipped.
!------------------------------------------------------------------------------

  DO i=1,npes
    IF(i == 1) THEN  ! local data
      IF(numpe == 1) THEN
        ipos = 0
        INQUIRE(UNIT=10,POS=ipos)
        ipos = ipos + (4 * nn * ndim) ; IF(verbose) PRINT *, "POS = ", ipos
        READ(10,POS=ipos) cbuffer     ; IF(verbose) PRINT *, cbuffer
        READ(10) nels_in              ; IF(verbose) PRINT *, nels_in
        iel = readCount(i)
        READ(10) g_num_pp(:,1:iel)
      END IF
    ELSE
      bufsize = readCount(i)*nod
      IF(numpe == 1) THEN
        iel = readCount(i)
        READ(10) g_num(:,1:iel)
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

  END SUBROUTINE READ_G_NUM_PP_BE
    
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

  SUBROUTINE READ_ELEMENTS_2(fname,npes,nn,numpe,g_num_pp)

    !/****f* input_output/read_elements
    !*  NAME
    !*    SUBROUTINE: read_elements
    !*  SYNOPSIS
    !*    Usage:      CALL read_elements(fname,npes,nn,numpe,g_num_pp)
    !*  FUNCTION
    !*    Master process reads the global array of elements and broadcasts
    !*    to slave processes.
    !*    Processes record only its local part of elements.
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    fname                 : Character
    !*                          : File name to read
    !*
    !*    npes                  : Integer
    !*                          : Number of processes
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
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    01.06.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  The method based on building global arrays.
    !*  
    !*  An improvement would be to avoid allocating the global array
    !* 
    !*  From an old version of ParaFEM. Works with program xx7.f90. Needs 
    !*  removing. Compare with READ_ELEMENTS
    !*/
  
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: fname
    INTEGER,      INTENT(IN)  :: npes, nn, numpe
    INTEGER,      INTENT(OUT) :: g_num_pp(:,:)
    INTEGER :: i, j, k, nod, nels_pp, iel, bufsize, ielpe, ier, buf2, &
               statu(MPI_STATUS_SIZE)
    INTEGER, ALLOCATABLE :: num_elem_pp(:), g_num(:,:)
  
    !----------------------------------------------------------------------
    ! 1. Allocate element information arrays
    !----------------------------------------------------------------------

    ALLOCATE(num_elem_pp(npes))

    !----------------------------------------------------------------------
    ! 2. Master process receives the integer nels_pp from slave processes
    !----------------------------------------------------------------------
    
    nels_pp = UBOUND(g_num_pp,2)

    IF (numpe==1) THEN
      num_elem_pp    = 0
      num_elem_pp(1) = nels_pp
    END IF
    bufsize = 1

    DO i = 2,npes
      IF(numpe==i) THEN
        CALL MPI_SEND(nels_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
      END IF
      IF(numpe==1) THEN
        CALL MPI_RECV(j,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
        num_elem_pp(i) = j
      END IF
    END DO 

    !----------------------------------------------------------------------
    ! 3. Allocate array to read elements
    !----------------------------------------------------------------------

    nod     = UBOUND(g_num_pp,1)

    ALLOCATE(g_num(nod,nels_pp))

    !----------------------------------------------------------------------
    ! 4. Master process reads the elements and sends them to slave processes
    !----------------------------------------------------------------------

    !----------------------------------------------------------------------
    ! 4.1 Master process reads its own elements
    !----------------------------------------------------------------------

    IF (numpe==1) THEN
      OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
      READ(10,*)    !header
      READ(10,*)    !header
      DO i = 1,nn   !read nodes until reaching the elements
        READ(10,*)
      END DO
      READ(10,*)    !keyword element
      DO iel = 1,num_elem_pp(1)
        READ(10,*)k,k,k,k,g_num_pp(:,iel),k
      END DO
    END IF

    !----------------------------------------------------------------------
    ! 4.2 Master process reads slave processes' elements and sends them
    !----------------------------------------------------------------------

    buf2 = nod*nels_pp

    DO i = 2,npes
      IF (numpe==1) THEN
        DO iel = 1,num_elem_pp(i)
          READ(10,*)k,k,k,k,g_num(:,iel),k
        END DO
        bufsize = nod*num_elem_pp(i)
        CALL MPI_SEND(g_num,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,ier)
      END IF      
      IF (numpe==i) THEN
        CALL MPI_RECV(g_num_pp,buf2,MPI_INTEGER,0,i,MPI_COMM_WORLD,statu,ier)
      END IF
    END DO
    CLOSE(10)
 
    !----------------------------------------------------------------------
    ! 5. Deallocate arrays
    !----------------------------------------------------------------------

    DEALLOCATE(g_num, num_elem_pp)
   
    RETURN

  END SUBROUTINE READ_ELEMENTS_2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE READ_ETYPE_PP(job_name,npes,numpe,etype_pp)
  
  !/****f* input/read_etype_pp
  !*  NAME
  !*    SUBROUTINE: read_etype_pp
  !*  SYNOPSIS
  !*    Usage:      CALL read_etype_pp(job_name,npes,numpe,etype_pp)
  !*  FUNCTION
  !*    The master process reads the element property type for each element 
  !*    and broadcasts that data to the slave processes. Each process only 
  !*    records its local data.
  !*  INPUTS
  !*    The following character argument has the INTENT(IN) attribute:
  !*
  !*    job_name              : Used to create file name to read
  !*
  !*    The following scalar integer arguments have the INTENT(IN) attribute:
  !*
  !*    npes                  : Integer
  !*                          : Total number of processors
  !*
  !*    numpe                 : Integer
  !*                          : Process number
  !*
  !*    The following scalar array argument has the INTENT(INOUT) attribute:
  !*
  !*    etype_pp(nels_pp)     : Element property type vector
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*    Ll.M. Evans
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2015
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  ParaView includes a bug which prevents the 'Material' section of the 
  !*  ENSIGHT Gold to be read and therefore integer MATID values
  !*  http://www.paraview.org/Bug/view.php?id=15151
  !*  http://www3.ensight.com/EnSight10_Docs/UserManual.pdf pp.713
  !*  Workaround - use MATIDs as reals and convert into int for ParaFEM
  !*/
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)     :: job_name
  CHARACTER(LEN=50)                 :: fname
  CHARACTER(LEN=80)                 :: cbuffer
  INTEGER, INTENT(IN)               :: npes, numpe
  INTEGER, INTENT(INOUT)            :: etype_pp(:)
  INTEGER                           :: nels_pp,i,ier,bufsize
  INTEGER                           :: int_in,etype_int,iel
  INTEGER                           :: readSteps,max_nels_pp
  INTEGER                           :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE              :: localCount(:),readCount(:)
  REAL(iwp), ALLOCATABLE            :: etype(:)
  REAL(iwp)                         :: etype_r

!------------------------------------------------------------------------------
! 1. Initiallize variables
!------------------------------------------------------------------------------

  nels_pp   = UBOUND(etype_pp,1)

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
  
  ALLOCATE(etype(max_nels_pp))
  
  etype    = 0.0_iwp
  etype_pp = 0
  
!------------------------------------------------------------------------------
! 4. Master processor opens the data file and advances to the start of the 
!    element steering array
!------------------------------------------------------------------------------
  
  IF (numpe==1) THEN
    fname     = job_name(1:INDEX(job_name, " ")-1) // ".ensi.MATID"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*)!   cbuffer !header
    READ(10,*)!   cbuffer !header
    READ(10,*)!   int_in  !header
    READ(10,*)!   cbuffer !header
  END IF
  
!------------------------------------------------------------------------------
! 5. Go around READSTEPS loop, read data, and send to appropriate processor
!------------------------------------------------------------------------------
  
  DO i=1,npes
    IF(i == 1) THEN  ! local data
      IF(numpe == 1) READ(10,*) etype(1:readCount(i))
      etype_pp(1:readCount(i))=int(etype(1:readCount(i)))
    ELSE
      bufsize = readCount(i)
      IF(numpe == 1) THEN
        etype = 0.0_iwp
        READ(10,*) etype(1:readCount(i))
!       etype_pp(1:readCount(i))=int(etype(1:readCount(i)))
        CALL MPI_SEND(etype(1:readCount(i)),bufsize,MPI_REAL8,i-1,i,           &
                      MPI_COMM_WORLD,status,ier)
      END IF
      IF(numpe == i) THEN
        CALL MPI_RECV(etype(1:readCount(i)),bufsize,MPI_REAL8,0,i,             &
                      MPI_COMM_WORLD,status,ier)
        etype_pp(1:readCount(i))=int(etype(1:readCount(i)))
      END IF
    END IF
  END DO
  
!------------------------------------------------------------------------------
! 6. Close file and deallocate global arrays
!------------------------------------------------------------------------------

  IF(numpe==1) CLOSE(10)

  DEALLOCATE(readCount,localCount,etype)

  RETURN

  END SUBROUTINE READ_ETYPE_PP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_ETYPE_PP_BE(job_name,npes,numpe,etype_pp)

  !/****f* input/read_etype_pp_be
  !*  NAME
  !*    SUBROUTINE: read_etype_pp_be
  !*  SYNOPSIS
  !*    Usage:      CALL read_etype_pp_be(job_name,npes,numpe,etype_pp)
  !*  FUNCTION
  !*    The master process reads the element property type for each element 
  !*    and broadcasts that data to the slave processes. Each process only 
  !*    records its local data.
  !*  INPUTS
  !*    The following character argument has the INTENT(IN) attribute:
  !*
  !*    job_name              : Used to create file name to read
  !*
  !*    The following scalar integer arguments have the INTENT(IN) attribute:
  !*
  !*    npes                  : Integer
  !*                          : Total number of processors
  !*
  !*    numpe                 : Integer
  !*                          : Process number
  !*
  !*    The following scalar array argument has the INTENT(INOUT) attribute:
  !*
  !*    etype_pp(nels_pp)     : Element property type vector     
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*    Ll. Evans
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2014
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  This subroutine can be modified to use MPI I/O
  !*
  !*  ParaView includes a bug which prevents the 'Material' section of the 
  !*  ENSIGHT Gold to be read and therefore integer MATID values
  !*  http://www.paraview.org/Bug/view.php?id=15151
  !*  http://www3.ensight.com/EnSight10_Docs/UserManual.pdf pp.713
  !*  Workaround - use MATIDs as reals and convert into int for ParaFEM
  !*/
  
  USE, INTRINSIC :: ISO_C_BINDING

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)     :: job_name
  CHARACTER(LEN=50)                 :: fname
  CHARACTER(LEN=80,KIND=C_CHAR)     :: cbuffer
  INTEGER, INTENT(IN)               :: npes, numpe
  INTEGER, INTENT(INOUT)            :: etype_pp(:)
  INTEGER                           :: nels_pp,i,ier,bufsize
  INTEGER(KIND=C_INT)               :: int_in
  INTEGER                           :: readSteps,max_nels_pp
  INTEGER                           :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE              :: localCount(:),readCount(:)
  REAL(KIND=C_FLOAT), ALLOCATABLE   :: etype(:)

!------------------------------------------------------------------------------
! 1. Initiallize variables
!------------------------------------------------------------------------------

  nels_pp   = UBOUND(etype_pp,1)

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
! 3. Allocate the array etype into which the element property type vector is  
!    to be read.
!------------------------------------------------------------------------------

  max_nels_pp = MAXVAL(readCount,1)
  
  ALLOCATE(etype(max_nels_pp))
  
  etype    = 0.0_iwp
  etype_pp = 0
  
!------------------------------------------------------------------------------
! 4. The master process opens the data file
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname     = job_name(1:INDEX(job_name, " ")-1) // ".bin.ensi.MATID"
    OPEN(10,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',         &
                       ACCESS='STREAM')
    READ(10)   cbuffer !header
    READ(10)   cbuffer !header
    READ(10)   int_in  !header
    READ(10)   cbuffer !header
  END IF

!------------------------------------------------------------------------------
! 5. Go around READSTEPS loop, read data, and send to appropriate process
!------------------------------------------------------------------------------

  DO i=1,npes
    IF(i == 1) THEN  ! local data
      IF(numpe == 1) READ(10) etype(1:readCount(i))
      etype_pp(1:readCount(i))=int(etype(1:readCount(i)))
    ELSE
      bufsize = readCount(i)
      IF(numpe == 1) THEN
        etype = 0.0_iwp
        READ(10) etype(1:readCount(i))
!       etype_pp(1:readCount(i))=int(etype(1:readCount(i)))
        CALL MPI_SEND(etype(1:readCount(i)),bufsize,MPI_REAL4,i-1,i,           &
                      MPI_COMM_WORLD,status,ier)
      END IF
      IF(numpe == i) THEN
        CALL MPI_RECV(etype(1:readCount(i)),bufsize,MPI_REAL4,0,i,             &
                      MPI_COMM_WORLD,status,ier)
        etype_pp(1:readCount(i))=int(etype(1:readCount(i)))
      END IF
    END IF
  END DO
  
!------------------------------------------------------------------------------
! 6. Close file and deallocate global arrays
!------------------------------------------------------------------------------

  IF(numpe==1) CLOSE(10)

  DEALLOCATE(readCount,localCount,etype)

  RETURN

  END SUBROUTINE READ_ETYPE_PP_BE

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
    
  SUBROUTINE READ_G_NUM_PP2(job_name,iel_start,nels,nn,numpe,g_num_pp)

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

  END SUBROUTINE READ_G_NUM_PP2
                
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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE READ_RESTRAINTS(fname,numpe,rest)

    !/****f* input/read_restraints
    !*  NAME
    !*    SUBROUTINE: read_restraints
    !*  SYNOPSIS
    !*    Usage:      CALL read_restraints(fname,numpe,rest)
    !*  FUNCTION
    !*    Master process reads the global array of nodes with restrained 
    !*    degrees of freedom.
    !*    Master process broadcasts to slave processes.
    !*
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    fname                  : Character
    !*                           : File name to read
    !*
    !*    numpe                  : Integer
    !*                           : Process number
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    rest(nr,nodof+1)       : Integer
    !*                           : Nodes (column 1) and degrees of freedom
    !*                             restrained (columns 2,3,4). 
    !*                             The criterion is:  1 fixed, 0 free
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    01.06.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  The method is based on reading a global array, although it's not a
    !*  worry because the nodes with constraints are a small proportion of the 
    !*  whole mesh.
    !*
    !*  THIS SUBROUTINE IS ONLY USED BY PROGRAM XX7 AND NEEDS REMOVING FROM 
    !*  PARAFEM
    !*/

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: fname
    INTEGER,      INTENT(IN)  :: numpe
    INTEGER,      INTENT(OUT) :: rest(:,:)
    INTEGER :: i, nr, nodof1, bufsize, ier

    !----------------------------------------------------------------------
    ! 1. Master process reads the data
    !----------------------------------------------------------------------

    nr = UBOUND(rest,1)

    IF(numpe==1)THEN
      OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
      DO i = 1,nr
        READ(10,*)rest(i,:)
      END DO
      CLOSE(10)
    END IF

    !----------------------------------------------------------------------
    ! 2. Master process broadcasts the data to slave processes
    !----------------------------------------------------------------------

    nodof1  = UBOUND(rest,2)
    bufsize = nr*nodof1

    CALL MPI_BCAST(rest,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    RETURN

  END SUBROUTINE READ_RESTRAINTS
  
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
    INTEGER                       :: i,k,nmats,nvals,bufsize,ier,ielpe
    INTEGER,INTENT(IN)            :: numpe,npes
    REAL(iwp),INTENT(INOUT)       :: materialValues(:,:)
    CHARACTER(LEN=50), INTENT(in) :: fname
    CHARACTER(LEN=10)             :: keyword

    IF(numpe==1)THEN
     OPEN(21,FILE=fname, STATUS='OLD', ACTION='READ')
     !Read the *MATERIAL keyword, num of materials in file and num values per material line (2)
     READ(21,*) keyword, nmats, nvals
     !Read the material labels line
     READ(21,*)                  ! skip line

     PRINT *, "nmats =", nmats
     DO i = 1,nmats
       READ(21,*)k, materialValues(:,i)
 !      PRINT *, "materialValues = ", materialValues(:,i)
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

  SUBROUTINE READ_MATERIAL(JOB_NAME,MATERIALVALUES,NUMPE,NPES)

    !Master process reads the global array from a single file and 
    !sends a copy to each slave. (I/O method 1: refer to notes)
    !k is the material number in the file, it is read and discarded

    IMPLICIT NONE
    INTEGER                       :: i,k,nmats,nvals,bufsize,ier,ielpe
    INTEGER,INTENT(IN)            :: numpe,npes
    REAL(iwp),INTENT(INOUT)       :: materialValues(:,:)
    CHARACTER(LEN=50), INTENT(IN) :: job_name
    CHARACTER(LEN=50)             :: fname    
    CHARACTER(LEN=10)             :: keyword

    IF(numpe==1)THEN
     fname=job_name(1:INDEX(job_name," ")-1) // ".mat"
     OPEN(21,FILE=fname, STATUS='OLD', ACTION='READ')
     !Read the *MATERIAL keyword, num of materials in file and num values per material line (2)
     READ(21,*) keyword, nmats, nvals
     !Read the material labels line
     READ(21,*)                  ! skip line

     PRINT *, "nmats =", nmats
     DO i = 1,nmats
       READ(21,*)k, materialValues(:,i)
 !      PRINT *, "materialValues = ", materialValues(:,i)
     END DO
     CLOSE(21)
    END IF
    
    bufsize       = ubound(materialValues,1)*ubound(materialValues,2)
    CALL MPI_BCAST(materialValues,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
   
    RETURN

  END SUBROUTINE READ_MATERIAL

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_NELS_PP(job_name,iel_start,nels_pp,npes,numpe)

  !/****f* input/read_nels_pp
  !*  NAME
  !*    SUBROUTINE: read_nels_pp
  !*  SYNOPSIS
  !*    Usage:      CALL read_nels_pp(job_name,iel_start,nels_pp,npes,numpe)
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
  INTEGER, INTENT(OUT)             :: nels_pp,iel_start
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

  nels_pp   = psize(numpe)
  iel_start = 0

  IF(numpe==1) THEN
    iel_start = 1
  ELSE
    DO i = 2, numpe
      iel_start = iel_start + psize(i-1)
    END DO
    iel_start   = iel_start + 1
  END IF

  PRINT *, "iel_start = ", iel_start, " on process ", numpe 

  DEALLOCATE(psize)

  RETURN
  END SUBROUTINE READ_NELS_PP


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  SUBROUTINE READ_LF(job_name,numpe,lf)
  
  IMPLICIT NONE
  
  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  INTEGER, INTENT(IN)              :: numpe
  INTEGER                          :: i,incs
  REAL(iwp), INTENT(INOUT)         :: lf(:,:)
  
  CHARACTER(LEN=50)                :: fname

  lf=0.0_iwp
  incs= UBOUND(lf,2)
  
  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    ! skip 2 lines
    READ(10,*); READ(10,*)
    READ(10,*) lf
    CLOSE(10)
  END IF
  
  CALL MPI_BCAST(lf,incs*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
  
  END SUBROUTINE READ_LF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P121(job_name,numpe,e,element,limit,loaded_nodes,mesh,nels, &
                       nip,nn,nod,nr,partition,tol,v)

  !/****f* input/read_p121
  !*  NAME
  !*    SUBROUTINE: read_p121
  !*  SYNOPSIS
  !*    Usage:      CALL read_p121(job_name,numpe,e,element,limit,            &
  !*                               loaded_nodes,mesh,nels,nip,nn,nod,nr,
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
  !*    (c) University of Manchester 2010-2013
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: limit,mesh,partition 
  REAL(iwp), INTENT(INOUT)         :: e,v,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(9)
  REAL(iwp)                        :: real_store(3)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,nels,nn,nr,nip,nod,loaded_nodes,        &
               e,v,tol,limit
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
    integer_store(9)   = partition

    real_store         = 0.0_iwp

    real_store(1)      = e  
    real_store(2)      = v  
    real_store(3)      = tol  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 9
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
    limit           = integer_store(8)
    partition       = integer_store(9)

    e               = real_store(1)
    v               = real_store(2)
    tol             = real_store(3)

  END IF

  RETURN
  END SUBROUTINE READ_P121

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P122(job_name,numpe,c,cjits,cjtol,e,element,                &
    fixed_freedoms,loaded_nodes,incs,mesh,nels,nip,nn,nod,nr,phi,partition,   &
    plasits,plastol,psi,v)

  !/****f* input/read_p122
  !*  NAME
  !*    SUBROUTINE: read_p122
  !*  SYNOPSIS
  !*    Usage:      CALL read_p122(job_name,numpe,c,cjits,cjtol,e,            &
  !*                  element,loaded_nodes,incs,mesh,nels,nip,nn,nr,phi,      &
  !*                  partition,plasits,plastol,psi,v)
  !*
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
  !*    element                : Character
  !*                           : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*    fixed_freedoms         : Integer
  !*                           : Number of fixed displacements
  !*
  !*    kx                     : Real
  !*
  !*    ky                     : Real
  !*
  !*    kz                     : Real
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
  !*    partition              : Integer
  !*                           : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    08.03.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nod,nr,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: mesh,partition 
  INTEGER, INTENT(INOUT)           :: incs,plasits,cjits,fixed_freedoms
  REAL(iwp), INTENT(INOUT)         :: phi,c,psi,e,v,plastol,cjtol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(12)
  REAL(iwp)                        :: real_store(7)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,nels,nn,nr,nip,nod,fixed_freedoms, &
               loaded_nodes,phi,c,psi,e,v,incs,plasits,cjits,            &
               plastol,cjtol
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = partition
    integer_store(3)   = nels
    integer_store(4)   = nn
    integer_store(5)   = nr 
    integer_store(6)   = nip
    integer_store(7)   = nod
    integer_store(8)   = loaded_nodes
    integer_store(9)   = incs
    integer_store(10)  = plasits
    integer_store(11)  = cjits
    integer_store(12)  = fixed_freedoms

    real_store         = 0.0_iwp

    real_store(1)      = phi  
    real_store(2)      = c  
    real_store(3)      = psi  
    real_store(4)      = e  
    real_store(5)      = v  
    real_store(6)      = plastol
    real_store(7)      = cjtol    

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 12
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 7
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  bufsize = 15
  CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN
    mesh            = integer_store(1)
    partition       = integer_store(2)
    nels            = integer_store(3)
    nn              = integer_store(4)
    nr              = integer_store(5)
    nip             = integer_store(6)
    nod             = integer_store(7)
    loaded_nodes    = integer_store(8)
    incs            = integer_store(9)
    plasits         = integer_store(10)
    cjits           = integer_store(11)
    fixed_freedoms  = integer_store(12)

    phi             = real_store(1)
    c               = real_store(2)
    psi             = real_store(3)
    e               = real_store(4)
    v               = real_store(5)
    plastol         = real_store(6)
    cjtol           = real_store(7)

  END IF
  
  RETURN
  END SUBROUTINE READ_P122

  SUBROUTINE READ_QINC(job_name,numpe,qinc)
  
  IMPLICIT NONE
  
  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  INTEGER, INTENT(IN)              :: numpe
  INTEGER                          :: i,incs
  REAL(iwp), INTENT(INOUT)         :: qinc(:)
  
  CHARACTER(LEN=50)                :: fname

  qinc=0.0_iwp
  incs= UBOUND(qinc,1)
  
  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    ! skip 6 lines
    READ(10,*); READ(10,*); READ(10,*); READ(10,*); READ(10,*); READ(10,*)
    DO i=1,incs
      READ(10,*) qinc(i)
    END DO
    CLOSE(10)
  END IF
  
  CALL MPI_BCAST(qinc,incs,MPI_REAL8,0,MPI_COMM_WORLD,ier)
  
  END SUBROUTINE READ_QINC
  


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


  SUBROUTINE READ_P123(job_name,numpe,element,fixed_freedoms,kx,ky,kz,limit,  &
                       loaded_nodes,mesh,nels,nip,nn,nod,nr,nres,partition,tol)

  !/****f* input/read_p123
  !*  NAME
  !*    SUBROUTINE: read_p123
  !*  SYNOPSIS
  !*    Usage:      CALL read_p123(job_name,numpe,element,fixed_freedoms,
  !*                               kx,ky,kz,limit,loaded_nodes,mesh,nels,nip,
  !*                               nn,nod,nr,nres,partition,tol,v)
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
  !*    element                : Character
  !*                           : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*    fixed_freedoms         : Integer
  !*                           : Number of fixed displacements
  !*
  !*    kx                     : Real
  !*
  !*    ky                     : Real
  !*
  !*    kz                     : Real
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
  !*    nres                   : Integer
  !*                           : Equation to output summary data
  !*
  !*    partition              : Integer
  !*                           : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*
  !*    tol                    : Real
  !*                           : Tolerance for PCG
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    08.03.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012-2013
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nres,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: limit,mesh,fixed_freedoms,partition 
  REAL(iwp), INTENT(INOUT)         :: kx,ky,kz,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(11)
  REAL(iwp)                        :: real_store(4)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,nels,nn,nr,nip,nod,loaded_nodes,        &
               fixed_freedoms,kx,ky,kz,tol,limit,nres
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
    integer_store(11)  = nres

    real_store         = 0.0_iwp

    real_store(1)      = kx  
    real_store(2)      = ky  
    real_store(3)      = kz  
    real_store(4)      = tol  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 11
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 4
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
    nres            = integer_store(11)

    kx              = real_store(1)
    ky              = real_store(2)
    kz              = real_store(3)
    tol             = real_store(4)

  END IF

  !Commented out below lines to test whether loded nodes and fixed freedoms work together
  
  !IF(fixed_freedoms > 0 .AND. loaded_nodes > 0) THEN
  !  PRINT *
  !  PRINT *, "Error - model has", fixed_freedoms, " fixed freedoms and"
  !  PRINT *, loaded_nodes, " loaded nodes"
  !  PRINT *, "Mixed displacement and load control not supported"
  !  PRINT *
  !  CALL shutdown()
  !END IF

  RETURN
  END SUBROUTINE READ_P123
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P124_4(job_name,numpe,dtim,element,fixed_freedoms,kx,ky,kz, &
                       limit,loaded_nodes,mesh,nels,nip,nn,nod,npri,nr,       &
                       nstep,partition,theta,tol)

  !/****f* input/read_p124_4
  !*  NAME
  !*    SUBROUTINE: read_p124_4
  !*  SYNOPSIS
  !*    Usage:      CALL read_p124_4(job_name,numpe,dtim,element,fixed_freedoms,
  !*                               kx,ky,kz,limit,loaded_nodes,mesh,nels,nip,
  !*                               nn,nod,npri,nr,nstep,partition,theta,tol)
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer has the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor number
  !*
  !*    The following scalar character has the INTENT(INOUT) attribute:
  !*
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*    The following scalar integers have the INTENT(INOUT) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    limit                  : Maximum number of PCG iterations allowed
  !*    loaded_nodes           : Number of nodes with applied forces
  !*    mesh                   : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    nels                   : Total number of elements
  !*    nip                    : Number of Gauss integration points
  !*    nn                     : Total number of nodes in the mesh
  !*    nod                    : Number of nodes per element
  !*    npri                   : Number of timesteps to skip before printing
  !*    nr                     : Number of nodes with restrained degrees of
  !*                             freedom 
  !*    nstep                  : Number of steps to complete in the simulation  
  !*    partition              : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*
  !*    The following scalar reals have the INTENT(INOUT) attribute:
  !*
  !*    dtim                   : Time step
  !*    kx                     : Conductivity in x-direction
  !*    ky                     : Conductivity in y-direction
  !*    kz                     : Conductivity in z-direction
  !*    theta                  : Parameter in theta integrator
  !*    tol                    : Tolerance for PCG
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    25.06.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012-2013
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*  This was called READ_P124 
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: limit,mesh,fixed_freedoms,partition 
  INTEGER, INTENT(INOUT)           :: npri,nstep
  REAL(iwp), INTENT(INOUT)         :: kx,ky,kz,tol
  REAL(iwp), INTENT(INOUT)         :: dtim,theta

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(12)
  REAL(iwp)                        :: real_store(6)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,nels,nn,nr,nip,nod,loaded_nodes,        &
               fixed_freedoms,kx,ky,kz,tol,limit,                             &
               dtim,nstep,npri,theta
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
    integer_store(11)  = npri
    integer_store(12)  = nstep

    real_store         = 0.0_iwp

    real_store(1)      = kx  
    real_store(2)      = ky  
    real_store(3)      = kz  
    real_store(4)      = tol  
    real_store(5)      = dtim  
    real_store(6)      = theta

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 12
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 6
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
    npri            = integer_store(11)
    nstep           = integer_store(12)

    kx              = real_store(1)
    ky              = real_store(2)
    kz              = real_store(3)
    tol             = real_store(4)
    dtim            = real_store(5)
    theta           = real_store(6)

  END IF

  !Commented out below lines to test whether loded nodes and fixed freedoms work together
  
  !IF(fixed_freedoms > 0 .AND. loaded_nodes > 0) THEN
  !  PRINT *
  !  PRINT *, "Error - model has", fixed_freedoms, " fixed freedoms and"
  !  PRINT *, loaded_nodes, " loaded nodes"
  !  PRINT *, "Mixed displacement and load control not supported"
  !  PRINT *
  !  CALL shutdown()
  !END IF

  RETURN
  END SUBROUTINE READ_P124_4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P124(job_name,numpe,dtim,element,fixed_freedoms,            &
                       limit,loaded_nodes,mesh,nels,nip,nn,nod,npri,nr,nres,  &
                       nstep,partition,theta,tol,np_types,val0)

  !/****f* input/read_p124
  !*  NAME
  !*    SUBROUTINE: read_p124
  !*  SYNOPSIS
  !*    Usage:      CALL read_p124(job_name,numpe,dtim,element,fixed_freedoms,
  !*                               kx,ky,kz,limit,loaded_nodes,mesh,nels,nip,
  !*                               nn,nod,npri,nr,nres,nstep,partition,theta,
  !*                               tol,np_types,rho,cp,val0)
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer has the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor number
  !*
  !*    The following scalar character has the INTENT(INOUT) attribute:
  !*
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*    The following scalar integers have the INTENT(INOUT) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    limit                  : Maximum number of PCG iterations allowed
  !*    loaded_nodes           : Number of nodes with applied forces
  !*    mesh                   : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    nels                   : Total number of elements
  !*    nip                    : Number of Gauss integration points
  !*    nn                     : Total number of nodes in the mesh
  !*    nod                    : Number of nodes per element
  !*    npri                   : Number of timesteps to skip before printing
  !*    nr                     : Number of nodes with restrained degrees of
  !*                             freedom 
  !*    nres                   : Equation number for summary results 
  !*    nstep                  : Number of steps to complete in the simulation  
  !*    partition              : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*    np_types               : Number of property types
  !*
  !*    The following scalar reals have the INTENT(INOUT) attribute:
  !*
  !*    dtim                   : Time step
  !*    theta                  : Parameter in theta integrator
  !*    tol                    : Tolerance for PCG
  !*    val0                   : Initial temperature of whole model
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Llion Marc Evans
  !*  CREATION DATE
  !*    21.08.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012-2013
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*  This is a modified copy of READ_XX12
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nres,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: limit,mesh,fixed_freedoms,partition 
  INTEGER, INTENT(INOUT)           :: npri,nstep
  INTEGER, INTENT(INOUT)           :: np_types
  REAL(iwp), INTENT(INOUT)         :: tol
  REAL(iwp), INTENT(INOUT)         :: dtim,theta
  REAL(iwp), INTENT(INOUT)         :: val0

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(14)
  REAL(iwp)                        :: real_store(4)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,np_types,nels,nn,nr,nip,nod,            &
               loaded_nodes,fixed_freedoms,val0,                              &
               dtim,nstep,npri,theta,tol,limit,nres
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
    integer_store(11)  = npri
    integer_store(12)  = nstep
    integer_store(13)  = np_types
    integer_store(14)  = nres

    real_store         = 0.0_iwp

    real_store(1)      = tol  
    real_store(2)      = dtim  
    real_store(3)      = theta
    real_store(4)      = val0

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 14
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 4
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
    npri            = integer_store(11)
    nstep           = integer_store(12)
    np_types        = integer_store(13)
    nres            = integer_store(14)

    tol             = real_store(1)
    dtim            = real_store(2)
    theta           = real_store(3)
    val0            = real_store(4)

  END IF

  RETURN
  END SUBROUTINE READ_P124

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P125(job_name,numpe,dtim,element,fixed_freedoms,kx,ky,kz,   &
                       loaded_nodes,mesh,nels,nip,nn,nod,npri,nr,nres,nstep,  &
                       partition,val0)

  !/****f* input/read_p125
  !*  NAME
  !*    SUBROUTINE: read_p125
  !*  SYNOPSIS
  !*    Usage:      CALL read_p125(job_name,numpe,dtim,element,fixed_freedoms,
  !*                               kx,ky,kz,loaded_nodes,mesh,nels,nip,nn,nod,
  !*                               npri,nr,nres,nstep,partition,val0)
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer has the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor number
  !*
  !*    The following scalar character has the INTENT(INOUT) attribute:
  !*
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*    The following scalar integers have the INTENT(INOUT) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    loaded_nodes           : Number of nodes with applied forces
  !*    mesh                   : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    nels                   : Total number of elements
  !*    nip                    : Number of Gauss integration points
  !*    nn                     : Total number of nodes in the mesh
  !*    nod                    : Number of nodes per element
  !*    npri                   : Number of timesteps to skip before printing
  !*    nres                   : Equation number for print out
  !*    nr                     : Number of nodes with restrained degrees of
  !*                             freedom 
  !*    nstep                  : Number of steps to complete in the simulation  
  !*    partition              : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*
  !*    The following scalar reals have the INTENT(INOUT) attribute:
  !*
  !*    dtim                   : Time step
  !*    kx                     : Conductivity in x-direction
  !*    ky                     : Conductivity in y-direction
  !*    kz                     : Conductivity in z-direction
  !*    val0                   : Initial value
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    25.06.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: mesh,fixed_freedoms,partition 
  INTEGER, INTENT(INOUT)           :: npri,nstep,nres
  REAL(iwp), INTENT(INOUT)         :: kx,ky,kz,dtim,val0

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(12)
  REAL(iwp)                        :: real_store(5)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,nels,nn,nr,nip,nod,loaded_nodes,        &
               fixed_freedoms,kx,ky,kz,dtim,nstep,npri,nres,val0
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
    integer_store(9)   = partition
    integer_store(10)  = npri
    integer_store(11)  = nstep
    integer_store(12)  = nres

    real_store         = 0.0_iwp

    real_store(1)      = kx  
    real_store(2)      = ky  
    real_store(3)      = kz  
    real_store(4)      = dtim  
    real_store(5)      = val0 

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 12
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 5
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
    partition       = integer_store(9)
    npri            = integer_store(10)
    nstep           = integer_store(11)
    nres            = integer_store(12)

    kx              = real_store(1)
    ky              = real_store(2)
    kz              = real_store(3)
    dtim            = real_store(4)
    val0            = real_store(5)

  END IF

  !Commented out below lines to test whether loded nodes and fixed freedoms work together
  
  !IF(fixed_freedoms > 0 .AND. loaded_nodes > 0) THEN
  !  PRINT *
  !  PRINT *, "Error - model has", fixed_freedoms, " fixed freedoms and"
  !  PRINT *, loaded_nodes, " loaded nodes"
  !  PRINT *, "Mixed displacement and load control not supported"
  !  PRINT *
  !  CALL shutdown()
  !END IF

  RETURN
  END SUBROUTINE READ_P125


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_P126(job_name,numpe,cjits,cjtol,ell,fixed_equations,        &
                       kappa,limit,mesh,nels,nip,nn,nr,nres,partitioner,      &
                       penalty,rho,tol,x0,visc) 

  !/****f* input/read_p126
  !*  NAME
  !*    SUBROUTINE: read_p126
  !*  SYNOPSIS
  !*    Usage:      CALL read_p126(job_name,numpe,cjits,cjtol,ell,
  !*                               fixed_equations,kappa,limit,mesh,          &
  !*                               nels,nip,nn,nr,nres,partitioner,penalty,   &
  !*                               rho,tol,x0,visc)
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
  !*    nres                   : Integer
  !*                           : Equation number for output
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
  !*    (c) University of Manchester 2010-13
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nres,nip,cjits,fixed_equations
  INTEGER, INTENT(INOUT)           :: limit,mesh,ell,partitioner 
  REAL(iwp), INTENT(INOUT)         :: cjtol,kappa,penalty,rho,tol,x0,visc

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(11)
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
    READ(10,*) mesh,partitioner,nels,nn,nres,nr,fixed_equations,nip,visc,rho, &
               tol,limit,cjtol,cjits,penalty,x0,ell,kappa
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
    integer_store(11)  = nres

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

  bufsize = 11 
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
    nres            = integer_store(11)

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

  SUBROUTINE READ_P127(job_name,numpe,nels,nxe,nze,aa,bb,cc,nip,kx,ky,kz,e,v, &
    dtim,nstep,theta,cjits,cjtol,nlfp)

  !/****f* input/read_p127
  !*  NAME
  !*    SUBROUTINE: read_p127
  !*  SYNOPSIS
  !*    Usage:      CALL read_p127(job_name,numpe,nels,nxe,nze,aa,bb,cc,nip,  &
  !*                               kx,ky,kz,e,v,dtim,nstep,theta,cjits,cjtol, &
  !*                               nlfp)
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer has the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor number
  !*
  !*    The following scalar character has the INTENT(INOUT) attribute:
  !*
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*    The following scalar integers have the INTENT(INOUT) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    loaded_nodes           : Number of nodes with applied forces
  !*    mesh                   : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    nels                   : Total number of elements
  !*    nip                    : Number of Gauss integration points
  !*    nn                     : Total number of nodes in the mesh
  !*    nod                    : Number of nodes per element
  !*    npri                   : Number of timesteps to skip before printing
  !*    nres                   : Equation number for print out
  !*    nr                     : Number of nodes with restrained degrees of
  !*                             freedom 
  !*    nstep                  : Number of steps to complete in the simulation  
  !*    partition              : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*    nlfp                   :
  !*
  !*    The following scalar reals have the INTENT(INOUT) attribute:
  !*
  !*    dtim                   : Time step
  !*    kx                     : Conductivity in x-direction
  !*    ky                     : Conductivity in y-direction
  !*    kz                     : Conductivity in z-direction
  !*    val0                   : Initial value
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    01.02.2013
  !*  COPYRIGHT
  !*    (c) University of Manchester 2013
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nxe,nze,nip,nstep,cjits,nlfp
  REAL(iwp), INTENT(INOUT)         :: kx,ky,kz,e,v,dtim,theta,cjtol,aa,bb,cc

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(7)
  REAL(iwp)                        :: real_store(11)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) nels,nxe,nze,aa,bb,cc,nip,kx,ky,kz,e,v,dtim,nstep,theta,      &
               cjits,cjtol,nlfp
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = nels
    integer_store(2)   = nxe
    integer_store(3)   = nze
    integer_store(4)   = nip 
    integer_store(5)   = nstep
    integer_store(6)   = cjits
    integer_store(7)   = nlfp

    real_store         = 0.0_iwp

    real_store(1)      = kx  
    real_store(2)      = ky  
    real_store(3)      = kz  
    real_store(4)      = e  
    real_store(5)      = v 
    real_store(6)      = dtim  
    real_store(7)      = theta  
    real_store(8)      = cjtol  
    real_store(9)      = aa  
    real_store(10)     = bb 
    real_store(11)     = cc 

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 7
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 11
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    nels            = integer_store(1)
    nxe             = integer_store(2)
    nze             = integer_store(3)
    nip             = integer_store(4)
    nstep           = integer_store(5)
    cjits           = integer_store(6)
    nlfp            = integer_store(7)

    kx              = real_store(1)
    ky              = real_store(2)
    kz              = real_store(3)
    e               = real_store(4)
    v               = real_store(5)
    dtim            = real_store(6)
    theta           = real_store(7)
    cjtol           = real_store(8)
    aa              = real_store(9)
    bb              = real_store(10)
    cc              = real_store(11)

  END IF

  RETURN
  END SUBROUTINE READ_P127
  

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE READ_P128(job_name,numpe,acc,e,el,er,lalfa,leig,lx,lz,mesh,      &
    nels,nip,nmodes,nn,nr,partitioner,rho,v)

  !/****f* input/read_p128
  !*  NAME
  !*    SUBROUTINE: read_p128
  !*  SYNOPSIS
  !*    Usage:      CALL read_p128(job_name,numpe,acc,e,el,er,lalfa,leig,lx,  &
  !*                               lz,mesh,nels,nip,nmodes,nn,nr,             &
  !*                               partitioner,rho,v)
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
  !*    04.02.2013
  !*  COPYRIGHT
  !*    (c) University of Manchester 2013
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=50)                :: fname
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nip,lalfa,leig,lx,lz
  INTEGER, INTENT(INOUT)           :: mesh,partitioner,nmodes
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,el,er,acc

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(11)
  REAL(iwp)                        :: real_store(6)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) mesh,partitioner,nels,nn,nr,nip,rho,e,v,nmodes,el,er,          &
      lalfa,leig,lx,lz,acc
    CLOSE(10)
   
    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nr 
    integer_store(5)   = nip
    integer_store(6)   = nmodes
    integer_store(7)   = lalfa
    integer_store(8)   = leig
    integer_store(9)   = lx
    integer_store(10)  = lz
    integer_store(11)  = partitioner

    real_store         = 0.0_iwp

    real_store(1)      = rho  
    real_store(2)      = e  
    real_store(3)      = v  
    real_store(4)      = acc  
    real_store(5)      = el  
    real_store(6)      = er   

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 11
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 6
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    mesh         = integer_store(1)
    nels         = integer_store(2)
    nn           = integer_store(3)
    nr           = integer_store(4)
    nip          = integer_store(5)
    nmodes       = integer_store(6)
    lalfa        = integer_store(7)
    leig         = integer_store(8)
    lx           = integer_store(9)
    lz           = integer_store(10)
    partitioner  = integer_store(11)
    
    rho          = real_store(1)
    e            = real_store(2)
    v            = real_store(3)
    acc          = real_store(4)
    el           = real_store(5)
    er           = real_store(6)

  END IF

  RETURN
  END SUBROUTINE READ_P128

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  
  SUBROUTINE READ_P129(job_name,numpe,alpha1,beta1,e,element,                 &
                       limit,loaded_nodes,mesh,nels,nip,nn,nod,               &
                       npri,nr,nres,nstep,omega,partitioner,rho,theta,tol,v)

  !/****f* input/read_p129
  !*  NAME
  !*    SUBROUTINE: read_p129
  !*  SYNOPSIS
  !*    Usage:      CALL read_p129(job_name,numpe,alpha1,beta1,e,element,     &
  !*                               limit,loaded_nodes,mesh,nels,nip,nn,nod,   &
  !*                               npri,nr,nres,nstep,omega,partitioner,rho,  &
  !*                               theta,tol,v)
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
  !*    (c) University of Manchester 2010-13
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  CHARACTER(LEN=50)                :: fname
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nres,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: nstep,npri,limit,mesh,partitioner 
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,alpha1,beta1,theta,omega,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(12)
  REAL(iwp)                        :: real_store(8)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partitioner,nels,nn,nr,nip,nod,loaded_nodes,nres, &
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
    integer_store(12)  = nres

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

  bufsize = 12
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
    nres         = integer_store(12)

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
                        nels,nip,nn,nod,npri,nr,nres,nstep,partitioner,pload, &
                        rho,sbary,v)

  !/****f* input/read_p1210
  !*  NAME
  !*    SUBROUTINE: read_p1210
  !*  SYNOPSIS
  !*    Usage:      CALL read_p1210(job_name,numpe,dtim,e,element,            &
  !*                                loaded_nodes,meshgen,nels,nip,nn,nod,npri,&
  !*                                nr,nres,nstep,partitioner,pload,rho,      &
  !*                                sbary,v)
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
  INTEGER, INTENT(INOUT)           :: nels,nip,nn,nod,npri,nr,nres,nstep
  INTEGER, INTENT(INOUT)           :: loaded_nodes,meshgen,partitioner 
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,sbary,dtim,pload

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(11)
  REAL(iwp)                        :: real_store(6)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,meshgen,partitioner,nels,nip,nn,nr,nod,loaded_nodes,   &
               nres,rho,e,v,sbary,dtim,nstep,npri,pload
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
    integer_store(11)  = nres

    real_store         = 0.0_iwp

    real_store(1)      = dtim
    real_store(2)      = e  
    real_store(3)      = rho  
    real_store(4)      = sbary
    real_store(5)      = v  
    real_store(6)      = pload  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 11
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
    nres          = integer_store(11)

    dtim          = real_store(1)
    e             = real_store(2)
    rho           = real_store(3)
    sbary         = real_store(4)
    v             = real_store(5)
    pload         = real_store(6)
    
  END IF

  RETURN
  END SUBROUTINE READ_P1210



!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_XX1(job_name,numpe,e,element,fixed_freedoms,limit,          &
                      loaded_nodes,mesh,nels,nip,nn,nod,nr,partition,tol,v)

  !/****f* input/read_xx1
  !*  NAME
  !*    SUBROUTINE: read_xx1
  !*  SYNOPSIS
  !*    Usage:      CALL read_xx1(job_name,numpe,e,element,fixed_freedoms,    &
  !*                              limit,loaded_nodes,mesh,nels,nip,nn,nod,    &
  !*                              nr,partition,tol,v)
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
  !*                           : Number of fixed freedoms
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
  !*    10.06.2013
  !*  COPYRIGHT
  !*    (c) University of Manchester 2013
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: limit,mesh,partition,fixed_freedoms 
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
    READ(10,*) element,mesh,partition,nels,nn,nr,nip,nod,fixed_freedoms,      &
               loaded_nodes,e,v,tol,limit
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
    integer_store(9)   = partition
    integer_store(10)  = fixed_freedoms

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
    limit           = integer_store(8)
    partition       = integer_store(9)
    fixed_freedoms  = integer_store(10)

    e               = real_store(1)
    v               = real_store(2)
    tol             = real_store(3)

  END IF

  RETURN
  END SUBROUTINE READ_XX1


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

    SUBROUTINE BCAST_INPUTDATA_XX5(numpe,npes,element,meshgen,nels,nn,nr,     &
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
    INTEGER, PARAMETER              :: ilength=4, rlength=8
!   INTEGER, PARAMETER              :: ilength=8, rlength=8
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
    END SUBROUTINE BCAST_INPUTDATA_XX5

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    SUBROUTINE CHECK_INPUTDATA_XX5(numpe,npes,element,meshgen,nels,nn,nr,   &
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
    END SUBROUTINE CHECK_INPUTDATA_XX5
  

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_XX6(job_name,numpe,bmat,e,element,maxitr,mesh,ncv,nels,     &
                      nev,nip,nn,nod,nr,rho,tol,v,which)

  !/****f* input/read_xx6
  !*  NAME
  !*    SUBROUTINE: read_xx6
  !*  SYNOPSIS
  !*    Usage:      CALL read_xx6(job_name,numpe,bmat,e,element,maxitr,       &
  !*                              mesh,ncv,nels,nev,nip,nn,nod,nr,rho,tol,    &
  !*                              v,which) 
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
  END SUBROUTINE READ_XX6


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    SUBROUTINE READ_DATA_XX7(fname,numpe,nels,nn,nr,loaded_nodes,fixed_nodes, &
                             nip,limit,tol,e,v,nod,num_load_steps,jump,tol2)

    !/****f* input_output/read_data_xx7
    !*  NAME
    !*    SUBROUTINE: read_data_xx7
    !*  SYNOPSIS
    !*    Usage:      CALL read_data_xx7(fname,numpe,nels,nn,nr,              &
    !*                                   loaded_nodes,fixed_nodes,nip,        &
    !*                                   limit,tol,e,v,nod,num_load_steps,    &
    !*                                   jump,tol2)
    !*  FUNCTION
    !*    Master process reads the general data of the problem
    !*    Master process broadcasts to slave processes.
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    fname                  : Character
    !*                           : File name to read
    !*
    !*    numpe                  : Integer
    !*                           : Process number
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    nels                   : Integer
    !*                           : Total number of elements
    !*
    !*    nn                     : Integer
    !*                           : Total number of nodes 
    !*
    !*    nr                     : Integer
    !*                           : Number of nodes with restrained degrees of
    !*                             freedom 
    !*
    !*    loaded_nodes           : Integer
    !*                           : Number of nodes with applied forces
    !*
    !*    fixed_nodes            : Integer
    !*                           : Number of restrained degrees of freedom 
    !*                             with a non-zero applied value
    !*
    !*    nip                    : Integer
    !*                           : Number of Gauss integration points
    !*
    !*    limit                  : Integer
    !*                           : Maximum number of PCG iterations allowed
    !*
    !*    tol                    : Real
    !*                           : Tolerance for PCG
    !*
    !*    e                      : Real
    !*                           : Young's modulus
    !*
    !*    v                      : Real
    !*                           : Poisson coefficient
    !*
    !*    nod                    : Integer
    !*                           : Number of nodes per element
    !*
    !*    num_load_steps         : Integer
    !*                           : Number of load steps
    !*
    !*    jump                   : Integer
    !*                           : Number of load steps to skip before writing
    !*                             results (periodically)
    !*
    !*    tol2                   : Real
    !*                           : Tolerance for Newton-Raphson loop
    !*
    !*  AUTHOR
    !*    Francisco Calvo
    !*    L. Margetts
    !*  CREATION DATE
    !*    01.06.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
  
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)  :: fname
    INTEGER,      INTENT(IN)  :: numpe
    INTEGER,      INTENT(OUT) :: nels, nn, nr, loaded_nodes, fixed_nodes, nip,&
                                 limit, nod, num_load_steps, jump
    REAL(iwp),    INTENT(OUT) :: tol, e, v, tol2
    INTEGER                   :: bufsize, ier, vec_integer(10)
    REAL(iwp)                 :: vec_real(4)

    !----------------------------------------------------------------------
    ! 1. Master process reads the data and builds the integer and real
    !    vectors with the data
    !----------------------------------------------------------------------

    IF (numpe==1) THEN
      OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
      READ(10,*)nels,nn,nr,loaded_nodes,fixed_nodes,nip
      READ(10,*)limit,tol,e,v
      READ(10,*)nod
      READ(10,*)num_load_steps,jump
      READ(10,*)tol2
      CLOSE(10)
      
      vec_integer(1)  = nels
      vec_integer(2)  = nn
      vec_integer(3)  = nr
      vec_integer(4)  = loaded_nodes
      vec_integer(5)  = fixed_nodes
      vec_integer(6)  = nip
      vec_integer(7)  = limit
      vec_real(1)     = tol
      vec_real(2)     = e
      vec_real(3)     = v
      vec_integer(8)  = nod
      vec_integer(9)  = num_load_steps
      vec_integer(10) = jump
      vec_real(4)     = tol2
      
    END IF

    !----------------------------------------------------------------------
    ! 2. Master process broadcasts the data to slave processes
    !----------------------------------------------------------------------

    bufsize = 10
    CALL MPI_BCAST(vec_integer,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    bufsize = 4
    CALL MPI_BCAST(vec_real,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

    !----------------------------------------------------------------------
    ! 3. Slave processes extract the variables back from the vectors
    !----------------------------------------------------------------------

    IF (numpe/=1) THEN
      nels           = vec_integer(1)
      nn             = vec_integer(2)
      nr             = vec_integer(3)
      loaded_nodes   = vec_integer(4)
      fixed_nodes    = vec_integer(5)
      nip            = vec_integer(6)
      limit          = vec_integer(7)
      tol            = vec_real(1)
      e              = vec_real(2)
      v              = vec_real(3)
      nod            = vec_integer(8)
      num_load_steps = vec_integer(9)
      jump           = vec_integer(10)
      tol2           = vec_real(4)
    END IF

    RETURN

  END SUBROUTINE READ_DATA_XX7


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_XX12(job_name,numpe,element,fixed_freedoms,limit,            &
                       loaded_nodes,mesh,nels,nip,nn,nod,ntime,nr,partition,   &
                       theta,tol,np_types,chk,val0,el_print,                   &
					   xx12_input,xx12_output)

  !/****f* input/read_xx12
  !*  NAME
  !*    SUBROUTINE: read_xx12
  !*  SYNOPSIS
  !*    Usage:      CALL read_xx12(job_name,numpe,dtim,element,fixed_freedoms,
  !*                               kx,ky,kz,limit,loaded_nodes,mesh,nels,nip,
  !*                               nn,nod,npri,npri_chk,nr,nstep,partition,
  !*                               theta,tol,np_types,rho,cp,val0,el_print,
  !*                               xx12_input,xx12_output)
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts 
  !*    it to the slave processors.
  !*  INPUTS
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer has the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor number
  !*
  !*    The following scalar character has the INTENT(INOUT) attribute:
  !*
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*    chk                    : Checkpoint flag
  !*                           : Values: 'initialise' or 'restart'
  !*
  !*    The following scalar integers have the INTENT(INOUT) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    limit                  : Maximum number of PCG iterations allowed
  !*    loaded_nodes           : Number of nodes with applied forces
  !*    mesh                   : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    nels                   : Total number of elements
  !*    nip                    : Number of Gauss integration points
  !*    nn                     : Total number of nodes in the mesh
  !*    nod                    : Number of nodes per element
  !*    npri                   : Number of timesteps to skip before printing
  !*    nr                     : Number of nodes with restrained degrees of
  !*                             freedom 
  !*    nstep                  : Number of steps to complete in the simulation  
  !*    partition              : Type of partitioning 
  !*                           : 1 = internal partitioning
  !*                           : 2 = external partitioning with .psize file
  !*    np_types               : Number of property types
  !*    el_print               : Element number to print its values to .ttr2
  !*    xx12_input             : 1 = input in binary
  !*                           : 2 = input in plain text
  !*    xx12_output            : 1 = output in binary
  !*                           : 2 = output in plain text
  !*
  !*    The following scalar reals have the INTENT(INOUT) attribute:
  !*
  !*    dtim                   : Time step
  !*    theta                  : Parameter in theta integrator
  !*    tol                    : Tolerance for PCG
  !*    val0                   : Initial temperature of whole model
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Llion Marc Evans
  !*  CREATION DATE
  !*    21.08.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)       :: job_name
  CHARACTER(LEN=15), INTENT(INOUT)    :: element,chk
  INTEGER, INTENT(IN)                 :: numpe
  INTEGER, INTENT(INOUT)              :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)              :: limit,mesh,fixed_freedoms,partition 
  INTEGER, INTENT(INOUT)              :: el_print,xx12_input,xx12_output
  INTEGER, INTENT(INOUT)              :: np_types,ntime
  REAL(iwp), INTENT(INOUT)            :: tol,theta,val0
  
!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(15),i
  REAL(iwp)                        :: real_store(3)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,np_types,nels,nn,nr,nip,nod,loaded_nodes, &
               fixed_freedoms,chk,val0,ntime,theta,tol,limit,el_print,          &
			   xx12_input,xx12_output
    
!    ALLOCATE(dtim(nsection),npri(nsection),npri_chk(nsection),nstep(nsection))    
!    dtim     = 0.0_iwp
!    nstep    = 0
!    npri     = 0
!    npri_chk = 0
!    DO i=1,nsection
!      READ(10,*) dtim(i),nstep(i),npri(i),npri_chk(i)
!    END DO    
!    READ(10,*) theta,tol,limit,el_print,i_o

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
    integer_store(12)  = el_print
    integer_store(13)  = xx12_input
    integer_store(14)  = xx12_output
    integer_store(15)  = ntime

!    integer_store(11)  = npri
!    integer_store(12)  = npri_chk
!    integer_store(13)  = nstep
!    real_store(2)      = dtim

    real_store         = 0.0_iwp

    real_store(1)      = tol
    real_store(2)      = theta
    real_store(3)      = val0

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------
  
!  bufsize = nsection
!  CALL MPI_BCAST(npri,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!  CALL MPI_BCAST(npri_chk,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!  CALL MPI_BCAST(nstep,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!  CALL MPI_BCAST(dtim,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
  
  bufsize = 15
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 3
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  bufsize = 15
  CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  
  bufsize = 15
  CALL MPI_BCAST(chk,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  
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
    el_print        = integer_store(12)
    xx12_input      = integer_store(13)
	xx12_output     = integer_store(14)
    ntime           = integer_store(15)

!    npri            = integer_store(11)
!    npri_chk        = integer_store(12)
!    nstep           = integer_store(13)
!    dtim            = real_store(2)

    tol             = real_store(1)
    theta           = real_store(2)
    val0            = real_store(3)

  END IF

  !Commented out below lines to test whether loded nodes and fixed freedoms work together
  
  !IF(fixed_freedoms > 0 .AND. loaded_nodes > 0) THEN
  !  PRINT *
  !  PRINT *, "Error - model has", fixed_freedoms, " fixed freedoms and"
  !  PRINT *, loaded_nodes, " loaded nodes"
  !  PRINT *, "Mixed displacement and load control not supported"
  !  PRINT *
  !  CALL shutdown()
  !END IF

  RETURN
  END SUBROUTINE READ_xx12
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_TIMESTEPS(TIMESTEPS_REAL,TIMESTEPS_INT,FNAME,NUMPE,NPES)

    !Master process reads the global array from a single file and 
    !sends a copy to each slave. (I/O method 1: refer to notes)
    !k is the material number in the file, it is read and discarded

    IMPLICIT NONE
    INTEGER                       :: i,j,k,ntime,n_int,nreal,bufsize,ier,ielpe
    INTEGER,INTENT(IN)            :: numpe,npes
!    INTEGER,ALLOCATABLE,INTENT(INOUT)         :: timesteps_int(:,:)
!    REAL(iwp),ALLOCATABLE,INTENT(INOUT)       :: timesteps_real(:,:)
    INTEGER,INTENT(INOUT)         :: timesteps_int(:,:)
    REAL(iwp),INTENT(INOUT)       :: timesteps_real(:,:)
    CHARACTER(LEN=50), INTENT(in) :: fname
    CHARACTER(LEN=10)             :: keyword

    IF(numpe==1)THEN
     OPEN(21,FILE=fname, STATUS='OLD', ACTION='READ')
     !Read the *MATERIAL keyword, num of materials in file and num values per material line (2)
     READ(21,*) keyword, ntime, nreal, n_int
     !Read the material labels line
     READ(21,*)                  ! skip line

     PRINT *, "ntime =", ntime
     PRINT *, "nreal =", nreal
     PRINT *, "n_int =", n_int
     
!     ALLOCATE(timesteps_real(ntime,nreal),timesteps_int(ntime,n_int))
     
     DO i = 1,ntime
       READ(21,*)k, timesteps_real(i,:), timesteps_int(i,:)
       PRINT *, "k              =", k
       PRINT *, "timesteps_real =", timesteps_real(i,:)
       PRINT *, "timesteps_int  =", timesteps_int(i,:)
     END DO
     
!     DO i = 1,ntime
!       READ(21,'(I1)',ADVANCE='no')k
!       PRINT *, "k =", k
!       DO j = 1,nreal
!         READ(21,'(F5.3)',ADVANCE='no')timesteps_real(i,j)
!         PRINT *, "timesteps_real =", timesteps_real(i,j)
!       END DO
!       DO j = 1,n_int
!         READ(21,*)timesteps_int(i,j)
!         PRINT *, "timesteps_int =", timesteps_int(i,j)
!       END DO
!     END DO
     CLOSE(21)
    END IF
    
    bufsize       = ubound(timesteps_real,1)*ubound(timesteps_real,2)
    CALL MPI_BCAST(timesteps_real,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
    bufsize       = ubound(timesteps_int,1)*ubound(timesteps_int,2)
    CALL MPI_BCAST(timesteps_int,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
   
    RETURN

  END SUBROUTINE READ_TIMESTEPS
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  
  SUBROUTINE READ_X_PP(job_name,npes,numpe,j_chk,x_pp)

  !/****f* input/read_x_pp
  !*  NAME
  !*    SUBROUTINE: read_x_pp
  !*  SYNOPSIS
  !*    Usage:      CALL read_x_pp(job_name,ieq_start,nn,npes,numpe,nod,    &
  !*                               j_chk,x_pp)
  !*  FUNCTION
  !*    Master process reads the global array of elements IDs and broadcasts
  !*    to slave processes.
  !*    Processes record only its local part of elements.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    job_name              : Character
  !*                          : Used to create file name to read
  !*
  !*    ieq_start             : Integer
  !*                          : First equation number in a process
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
  !*    nod                   : Integer
  !*                          : Number of nodes per element
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    The following scalar array arguments have the INTENT(OUT) attribute:
  !*
  !*    x_pp(neq_pp)          : Initial temperatures vector
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*    Ll.M. Evans
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2014
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*
  !*/

  IMPLICIT NONE
  
  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=50)                :: fname,stepnum,cbuffer
  INTEGER, INTENT(IN)              :: npes, numpe
!  INTEGER, INTENT(IN)              :: ieq_start, nn, nod
  INTEGER                          :: neq_pp, ieq, i, k
  INTEGER                          :: bufsize, ier, int_in
  INTEGER                          :: readSteps,max_neq_pp
  INTEGER                          :: status(MPI_STATUS_SIZE)
  INTEGER, INTENT(INOUT)           :: j_chk
  REAL(iwp), INTENT(INOUT)         :: x_pp(:) 
  REAL(iwp), ALLOCATABLE           :: x_temp(:)
  INTEGER, ALLOCATABLE             :: localCount(:),readCount(:)
  LOGICAL                          :: verbose=.true.
! LOGICAL                           :: verbose=.false.
  
!------------------------------------------------------------------------------
! 1. Initialize variables
!------------------------------------------------------------------------------
  
  neq_pp = UBOUND(x_pp,1)
  
!------------------------------------------------------------------------------
! 2. Find READSTEPS, the number of steps in which the read will be carried
!    out and READCOUNT, the size of each read.
!------------------------------------------------------------------------------

  ALLOCATE(readCount(npes))
  ALLOCATE(localCount(npes))
  
  readCount         = 0
  localCount        = 0
  readSteps         = npes
  localCount(numpe) = neq_pp
  
  CALL MPI_ALLREDUCE(localCount,readCount,npes,MPI_INTEGER,MPI_SUM,           &
  MPI_COMM_WORLD,ier)
  
!------------------------------------------------------------------------------
! 3. Allocate the array x_temp into which the steering array for each processor
!    is to be read
!------------------------------------------------------------------------------

  max_neq_pp = MAXVAL(readCount,1)
  
  ALLOCATE(x_temp(max_neq_pp))
  
  x_temp = 0             ! different value for each processor
  
!------------------------------------------------------------------------------
! 4. Master processor opens the data file and advances to the start of the 
!    element steering array
!------------------------------------------------------------------------------
  
  fname     = job_name(1:INDEX(job_name, " ")-1) // ".chk"

  IF (numpe==1) THEN
    OPEN(10,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',         &
                       ACCESS='STREAM')
    READ(10)   j_chk
    READ(10)   cbuffer
    READ(10)   int_in
  END IF
  
!------------------------------------------------------------------------------
! 5. Go around READSTEPS loop, read data, and send to appropriate processor
!------------------------------------------------------------------------------

  DO i=1,npes
    IF(i == 1) THEN  ! local data
      IF(numpe == 1) THEN
        DO ieq = 1,readCount(i)
          READ(10)x_temp(ieq)
!          PRINT*,x_temp(ieq)
          x_pp(ieq) = x_temp(ieq)
        END DO
      END IF
    ELSE
      bufsize = readCount(i)
      IF(numpe == 1) THEN
        DO ieq = 1,readCount(i)
          READ(10)x_temp(ieq)
!          PRINT*,x_temp(ieq)
        END DO
        CALL MPI_SEND(x_temp(1:readCount(i)),bufsize,MPI_REAL8,i-1,i,      &
                      MPI_COMM_WORLD,status,ier)
      END IF
      IF(numpe == i) THEN
        x_temp = 0
        CALL MPI_RECV(x_temp,bufsize,MPI_REAL8,0,i,MPI_COMM_WORLD,status,ier)
        x_pp = x_temp
      END IF
    END IF
  END DO
  
!------------------------------------------------------------------------------
! 6. Close file and deallocate global arrays
!------------------------------------------------------------------------------

  IF(numpe==1) CLOSE(10)

  DEALLOCATE(x_temp,readCount,localCount)

  RETURN

  END SUBROUTINE READ_X_PP
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE READ_AMPLITUDE(JOB_NAME,NUMPE,NSTEP,VALUE)
  
  !/****f* input/read_amplitude
  !*  NAME
  !*    SUBROUTINE: read_amplitude
  !*  SYNOPSIS
  !*    Usage:      CALL read_amplitude(job_name,numpe,nstep,value)
  !*  FUNCTION
  !*    Master processor reads the global array of nodal force amplitudes 
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
  !*    nstep                    : Integer
  !*                             : Total number of time steps
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    value(load_amp)          : Real
  !*                             : Amplitude of force applied for given time step
  !*                             : Factor value to be multiplied with loaded_freedoms
  !*  AUTHOR
  !*    Lee Margetts
  !*    Llion Marc Evans
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010
  !******
  !*  The method is modified version of the subroutine read_loads
  !*
  !*/
  
  IMPLICIT NONE
  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname  
  INTEGER                       :: i,bufsize,ier,ielpe
  INTEGER,INTENT(IN)            :: numpe,nstep
  REAL(iwp),INTENT(INOUT)       :: value(:)
  
  IF(numpe==1)THEN
    fname = job_name(1:INDEX(job_name, " ")-1) // ".amp"
    OPEN(24, FILE=fname, STATUS='OLD', ACTION='READ')
    DO i = 1,nstep
      READ(24,*) value(i)
    END DO
    CLOSE(24)
  END IF
  
  bufsize = nstep
  CALL MPI_BCAST(value,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
  
  RETURN
  
  END SUBROUTINE READ_AMPLITUDE

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE bcast_inputdata_p123 (numpe,npes,nels,nxe,nze,nip,              &
                aa,bb,cc,kx,ky,kz,tol,limit,loaded_freedoms,fixed_freedoms)
  
  IMPLICIT NONE
  
  
  integer, INTENT(INOUT)   :: numpe,npes,nels,nxe,nze,nip,limit,             &
                              loaded_freedoms,fixed_freedoms
  INTEGER, PARAMETER       :: ilength=4, rlength=8
  INTEGER                  :: bufsizer,bufdecl
  REAL(iwp), INTENT(INOUT) :: aa,bb,cc,tol,kx,ky,kz
  !
  ! Assign temporary buffer for broadcast data
  !
  
  bufsizer=7*ilength + 7*rlength
  
  !------------------------------------------------------------------------------
  !---- Broadcast buffersize to allow slaves to allocate tempbuf
  !------------------------------------------------------------------------------
  
  call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)
  
  bufdecl=bufsizer/4
  allocate(tempbuf(bufdecl))
  
  !------------------------------------------------------------------------------
  !---- Pack all data ready for broadcast, broadcast, receive and unpack
  !------------------------------------------------------------------------------
  
  if (numpe==npes) then
  
    position = 0
  
  !---- integers
  
  CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
        &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (limit,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (loaded_freedoms,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (fixed_freedoms,1,MPI_INTEGER,tempbuf,bufsizer,position,        &
       &      MPI_COMM_WORLD,ier) !
  
  !----- reals
  
  CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (tol,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (kx,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (ky,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_PACK (kz,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
       &      MPI_COMM_WORLD,ier) !
  end if
  
  call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)
  
  if (numpe/=npes) then
  
  position=0
  
  !---- integers
  
  CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,limit,1,MPI_INTEGER,               &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,loaded_freedoms,1,MPI_INTEGER,     &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,fixed_freedoms,1,MPI_INTEGER,      &
       &      MPI_COMM_WORLD,ier) !
  
  !---- reals
  
  CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,tol,1,MPI_REAL8,                &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,kx,1,MPI_REAL8,                 &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,ky,1,MPI_REAL8,                 &
       &      MPI_COMM_WORLD,ier) !
  CALL MPI_UNPACK (tempbuf,bufsizer,position,kz,1,MPI_REAL8,                  &
       &      MPI_COMM_WORLD,ier) !
  
  end if
  
  END SUBROUTINE bcast_inputdata_p123

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE bcast_inputdata_p127 (numpe,npes,nels,nxe,nze, aa,bb,cc,       &
                        nip,kx,ky,kz,e,v,dtim,nstep,theta,cjits,cjtol)

IMPLICIT NONE


INTEGER,INTENT(INOUT)    :: numpe,npes,nels,nxe,nze,nip,nstep,cjits  
INTEGER,PARAMETER        :: ilength=8, rlength=8
INTEGER                  :: bufsizer,bufdecl
REAL(iwp), INTENT(INOUT) :: aa,bb,cc,kx,ky,kz,e,v,dtim,theta,cjtol

!
! Assign temporary buffer for broadcast data
!

bufsizer=6*ilength + 11*rlength

!------------------------------------------------------------------------------
!---- Broadcast buffersize to allow slaves to allocate tempbuf
!------------------------------------------------------------------------------

call MPI_BCAST(bufsizer,1,MPI_INTEGER,npes-1,MPI_COMM_WORLD,ier)

bufdecl=bufsizer/4
allocate(tempbuf(bufdecl))

!------------------------------------------------------------------------------
!---- Pack all data ready for broadcast, broadcast, receive and unpack
!------------------------------------------------------------------------------

if (numpe==npes) then

  position = 0

!---- integers

CALL MPI_PACK (nels,1,MPI_INTEGER,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nxe,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
      &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nze,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nip,1,MPI_INTEGER,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (nstep,1,MPI_INTEGER,tempbuf,bufsizer,position,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjits,1,MPI_INTEGER,tempbuf,bufsizer,position,       &
     &      MPI_COMM_WORLD,ier) !

!----- reals

CALL MPI_PACK (aa,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (bb,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cc,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kx,1,MPI_REAL8,tempbuf,bufsizer,position,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (ky,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (kz,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (e,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (v,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (dtim,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (theta,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_PACK (cjtol,1,MPI_REAL8,tempbuf,bufsizer,position,                   &
     &      MPI_COMM_WORLD,ier) !

end if

call MPI_BCAST(tempbuf,bufsizer,MPI_BYTE,npes-1,MPI_COMM_WORLD,ier)

if (numpe/=npes) then

position=0

!---- integers

CALL MPI_UNPACK (tempbuf,bufsizer,position,nels,1,MPI_INTEGER,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nxe,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nze,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nip,1,MPI_INTEGER,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,nstep,1,MPI_INTEGER,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjits,1,MPI_INTEGER,     &
     &      MPI_COMM_WORLD,ier) !

!---- reals

CALL MPI_UNPACK (tempbuf,bufsizer,position,aa,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,bb,1,MPI_REAL8,                 &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cc,1,MPI_REAL8,                  &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kx,1,MPI_REAL8,                &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,ky,1,MPI_REAL8,              &   
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,kz,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,e,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,v,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,dtim,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,theta,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !
CALL MPI_UNPACK (tempbuf,bufsizer,position,cjtol,1,MPI_REAL8,               &
     &      MPI_COMM_WORLD,ier) !

end if

END SUBROUTINE bcast_inputdata_p127

  

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_RFEMSOLVE(job_name,numpe,element,fixed_freedoms,limit,      &
                      loaded_nodes,mesh,mises,nels,nip,nn,nod,np_types,nr,    &
                      partition,tol)

  !/****f* input/read_rfemsolve
  !*  NAME
  !*    SUBROUTINE: read_rfemsolve
  !*  SYNOPSIS
  !*    Usage:      CALL read_rfemsolve(job_name,numpe,element,fixed_freedoms,
  !*                              limit,loaded_nodes,mesh,mises,nels,nip,nn,nod,
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
  !*    mises                  : Real
  !*                           : Threshold value for von Mises stress
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    12.06.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012
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
  REAL(iwp), INTENT(INOUT)         :: mises,tol

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(11)
  REAL(iwp)                        :: real_store(2)
  CHARACTER(LEN=50)                :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partition,np_types,nels,nn,nr,nip,nod,            &
               loaded_nodes,fixed_freedoms,tol,limit,mises
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

    real_store(1)      = mises
    real_store(2)      = tol  

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 11
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 2
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

    mises           = real_store(1)
    tol             = real_store(2)

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
  END SUBROUTINE READ_RFEMSOLVE
  

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


  SUBROUTINE MESH_ENSI(argv,nlen,g_coord,g_num,element,etype,nf,loads,        &
                       nstep,npri,dtim,solid)

   !/****f* input/mesh_ensi
   !*  NAME
   !*    SUBROUTINE: mesh_ensi
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf, &
   !*                               loads,nstep,npri,dtim,solid)
   !*  FUNCTION
   !*    This subroutine outputs a set of files in the Ensight gold format.
   !*    Models in this format can be viewed in ParaView.
   !* 
   !*    Element types supported:                Tested with:
   !*
   !*    2-node bar
   !*    3-node triangle                         p51  (4th edition p51_1.dat)
   !*    6-node triangle
   !*    4-node quadrilateral                    p115 (4th edition)
   !*    8-node quadrilateral                    p116 (4th edition)
   !*    4-node tetrahedron                      p54  (4th edition p54_2.dat)
   !*    8-node hexahedron                       p86  (4th edition)       
   !*    20-node hexahedron                      p55  (4th edition)
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    npri             : print interval
   !*	 nstep            : number of time steps in analysis
   !*
   !*    Scalar reals
   !*    dtim             : time step
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*	 element          : element type
   !*
   !*    Scalar logicals
   !*    solid            : type of analysis solid if .true. fluid if .false.
   !*
   !*    Dynamic scalar arrays
   !*    g_num            : global element node numbers vector
   !*    etype            : element property type vector
   !*    nf               : nodal freedom matrix
   !* 
   !*    Dynamic real arrays
   !* 	 g_coord          : global nodal coordinates
   !*	 oldlds           : initial loads vector
   !*
   !*  OUTPUTS
   !*  AUTHOR
   !*    I.M. Smith
   !*    D.V. Griffiths
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*  Subroutine required by the 5th Edition of "Programming the Finite
   !*  Element Method". Take care when modifying
   !*/

    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen,nstep,npri
    INTEGER,   INTENT(IN)         :: g_num(:,:),etype(:),nf(:,:)
    INTEGER                       :: i,j,k,l,m,n,nfe,nod,nels,ndim,nn
    INTEGER                       :: prnwidth,remainder
    REAL(iwp), INTENT(IN)         :: g_coord(:,:),loads(:),dtim
    CHARACTER(LEN=15), INTENT(IN) :: argv,element  
    LOGICAL, INTENT(IN)           :: solid
    
  !------------------------------------------------------------------------------
  ! 1. Initialisation
  !------------------------------------------------------------------------------
  
    nn   = UBOUND(g_coord,2) ; ndim = UBOUND(g_coord,1)
    nels = UBOUND(g_num,2)   ; nod  = UBOUND(g_num,1)
  
  !------------------------------------------------------------------------------
  ! 2. Write case file
  !------------------------------------------------------------------------------
  
    OPEN(12,FILE=argv(1:nlen)//'.ensi.case')
  
    WRITE(12,'(A/A)')    "#", "# Post-processing file generated by subroutine &
                               &WRITE_ENSI in "
    WRITE(12,'(A,A,/A)') "#", " Smith, Griffiths and Margetts, 'Programming the &
                               &Finite Element Method',","# Wiley, 2013."        
    WRITE(12,'(A/A/A)')  "#","# Ensight Gold Format","#"
    WRITE(12,'(2A/A)')   "# Problem name: ",argv(1:nlen),"#"
    WRITE(12,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
    WRITE(12,'(2A/A)')   "model: 1  ",argv(1:nlen)//'.ensi.geo',"VARIABLE"
    WRITE(12,'(2A)')     "scalar per element:  material      ",                &
                          argv(1:nlen)//'.ensi.MATID'
    IF(solid) THEN
      WRITE(12,'(2A)')   "scalar per node:     restraint     ",                &
                          argv(1:nlen)//'.ensi.NDBND'
      WRITE(12,'(2A)')   "vector per node:     displacement  ",                &
                          argv(1:nlen)//'.ensi.DISPL-******'
    ELSE
      WRITE(12,'(2A)')   "scalar per node:     pressure      ",                &
                          argv(1:nlen)//'.ensi.PRESSURE-******'
    END IF
    WRITE(12,'(2A)')     "vector per node:     load          ",                &
                          argv(1:nlen)//'.ensi.NDLDS'
    WRITE(12,'(A/A)')     "TIME","time set:     1"
    WRITE(12,'(A,I5)')    "number of steps:",nstep/npri
    WRITE(12,'(A,I5)')    "filename start number:",npri
    WRITE(12,'(A,I5)')    "filename increment:",npri
    WRITE(12,'(A)')       "time values:"
    prnwidth  = 5
    remainder = mod(nstep/npri,prnwidth)
    n         = ((nstep/npri) - remainder)/prnwidth
    IF(nstep/npri<=prnwidth) THEN
      DO i=1,nstep,npri
        IF(i==nstep) THEN
          WRITE(12,'(E12.5)') i*dtim
        ELSE
          WRITE(12,'(E12.5)',ADVANCE='no') i*dtim
        END IF
      END DO
    ELSE
      IF(remainder==0) THEN
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
      ELSE
  !     DO j=1,n-1
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
        m = (n*prnwidth)+1
        l = (n*prnwidth)+remainder
        DO i=m,l
          IF(i==l) THEN
            WRITE(12,'(E12.5)') dtim*i
          ELSE
            WRITE(12,'(E12.5)',ADVANCE='no') dtim*i
          END IF
        END DO
      END IF
    END IF
   
    CLOSE(12)
  
  !------------------------------------------------------------------------------
  ! 3. Write geometry file
  !------------------------------------------------------------------------------
  
    OPEN(13,FILE=argv(1:nlen)//'.ensi.geo')
    WRITE(13,'(/2A)')   "Problem name: ", argv(1:nlen)
    WRITE(13,'(A/A/A)') "Geometry files","node id given","element id given"
    WRITE(13,'(A/A)')   "part","      1"
    IF(ndim==2) WRITE(13,'(A)') "2d-mesh"
    IF(ndim==3) WRITE(13,'(A)') "Volume Mesh"
    WRITE(13,'(A)')     "coordinates"
    
    WRITE(13,'(I10)') nn
    DO j=1,ndim
      DO i=1,nn  
        WRITE(13,'(E12.5)') g_coord(j,i)
      END DO
    END DO
  
    IF(ndim==2) THEN ! ensight requires zeros for the z-ordinate
      DO i=1,nn
        WRITE(13,'(A)') " 0.00000E+00"
      END DO
    END IF
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod)
          CASE(3)
            WRITE(13,'(A/I10)') "tria3", nels
            DO i = 1,nels
              WRITE(13,'(3I10)')g_num(3,i),g_num(2,i),g_num(1,i)
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod)
          CASE(4)
            WRITE(13,'(A/I10)') "quad4", nels
            DO i = 1,nels
              WRITE(13,'(4I10)')g_num(1,i),g_num(4,i),g_num(3,i),g_num(2,i)
            END DO
          CASE(8)
            WRITE(13,'(A/I10)') "quad8", nels
            DO i = 1,nels
              WRITE(13,'(8I10)')g_num(1,i),g_num(7,i),g_num(5,i),g_num(3,i),    &
                                g_num(8,i),g_num(6,i),g_num(4,i),g_num(2,i)
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod)
          CASE(8)
            WRITE(13,'(A/I10)') "hexa8", nels
            DO i = 1,nels
              WRITE(13,'(8I10)') g_num(1,i),g_num(4,i),g_num(8,i),g_num(5,i),   &
                                 g_num(2,i),g_num(3,i),g_num(7,i),g_num(6,i)
            END DO
          CASE(20)
            WRITE(13,'(A/I10)') "hexa20", nels
            DO i = 1,nels
              WRITE(13,'(20I10)')                                               &
                g_num(1,i), g_num(7,i), g_num(19,i),g_num(13,i),g_num(3,i),     &
                g_num(5,i), g_num(17,i),g_num(15,i),g_num(8,i), g_num(12,i),    &
                g_num(20,i),g_num(9,i), g_num(4,i), g_num(11,i),g_num(16,i),    &
                g_num(10,i),g_num(2,i), g_num(6,i), g_num(18,i),g_num(14,i) 
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            WRITE(13,'(A/I10)') "tetra4", nels
            DO i = 1,nels
              WRITE(13,'(4I10)') g_num(1,i),g_num(3,i),g_num(2,i),g_num(4,i)
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE DEFAULT
        WRITE(13,'(A)')       "# Element type not recognised"
    END SELECT
  
    CLOSE(13)
  
  !------------------------------------------------------------------------------
  ! 4. Write file containing material IDs
  !------------------------------------------------------------------------------
  
    OPEN(14,FILE=argv(1:nlen)//'.ensi.MATID')
    WRITE(14,'(A)') "Alya Ensight Gold --- Scalar per-element variable file"
    WRITE(14,'(A/A)') "part", "      1"
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod) 
          CASE(3)
            WRITE(14,'(A)') "tria3"
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod) 
          CASE(4)
            WRITE(14,'(A)') "quad4"
          CASE(8)
            WRITE(14,'(A)') "quad8"
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod) 
          CASE(8)
            WRITE(14,'(A)') "hexa8"
          CASE(20)
            WRITE(14,'(A)') "hexa20"
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            WRITE(14,'(A)') "tetra4"
          CASE DEFAULT
          WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE DEFAULT
        WRITE(14,'(A)')   "# Element type not recognised"
    END SELECT
   
    DO i=1,nels; WRITE(14,'(I10)') etype(i); END DO
  
    WRITE(14,'(A)')
  
    CLOSE(14)
  
  !------------------------------------------------------------------------------
  ! 5. Write boundary conditions. Encoded using formula: 4z + 2y + 1x
  !
  !    110 = 1   010 = 2   100 = 3   011 = 4   101 = 5   001 = 6   000 = 7
  !------------------------------------------------------------------------------
  
    IF(solid) THEN
      OPEN(15,FILE=argv(1:nlen)//'.ensi.NDBND')
      WRITE(15,'(A)')     "Alya Ensight Gold --- Scalar per-node variable file"
      WRITE(15,'(A/A/A)') "part", "      1","coordinates"
      IF(ndim==3) THEN
        DO i=1,UBOUND(g_coord,2) 
          nfe=0
          IF(nf(1,i)==0) nfe=nfe+1
          IF(nf(2,i)==0) nfe=nfe+2
          IF(nf(3,i)==0) nfe=nfe+4
          WRITE(15,'(I2)') nfe
        END DO
      ELSE IF(ndim==2) THEN
        DO i=1,nn
          nfe=0
          IF(nf(1,i)==0) nfe=nfe+1
          IF(nf(2,i)==0) nfe=nfe+2
          WRITE(15,'(I2)') nfe
        END DO
      ELSE
        PRINT *, "Wrong number of dimensions in mesh_ensi"
      END IF   
    END IF
  
    CLOSE(15)
  
  !------------------------------------------------------------------------------
  ! 6. Write loaded nodes
  !------------------------------------------------------------------------------
  
    OPEN(16,FILE=argv(1:nlen)//'.ensi.NDLDS')
    WRITE(16,'(A)')     "Alya Ensight Gold --- Vector per-node variable file"
    WRITE(16,'(A/A/A)') "part", "      1","coordinates"
    DO j=1,UBOUND(nf,1)
      DO i=1, UBOUND(nf,2)
        WRITE(16,'(E12.5)') loads(nf(j,i))
      END DO
    END DO
    CLOSE(16)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MESH_ENSI_BIN(argv,nlen,g_coord,g_num,element,etype,nf,loads,    &
                           nstep,npri,dtim,solid)

   !/****f* input/mesh_ensi_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_bin(argv,nlen,g_coord,g_num,element,      &
   !*                                   etype,nf,loads,nstep,npri,dtim,solid)
   !*  FUNCTION
   !*    This subroutine outputs a set of files in the C binary version of the 
   !*    Ensight gold format. Models in this format can be viewed in ParaView.
   !* 
   !*    Element types supported:                Tested with:
   !*
   !*    2-node bar
   !*    3-node triangle                         p51  (4th edition p51_1.dat)
   !*    6-node triangle
   !*    4-node quadrilateral                    p115 (4th edition)
   !*    8-node quadrilateral                    p116 (4th edition)
   !*    4-node tetrahedron                      p54  (4th edition p54_2.dat)
   !*    8-node hexahedron                       p86  (4th edition)       
   !*    20-node hexahedron                      p55  (4th edition)
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    npri             : print interval
   !*	 nstep            : number of time steps in analysis
   !*
   !*    Scalar reals
   !*    dtim             : time step
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*	 element          : element type
   !*
   !*    Scalar logicals
   !*    solid            : type of analysis solid if .true. fluid if .false.
   !*
   !*    Dynamic scalar arrays
   !*    g_num            : global element node numbers vector
   !*    nf               : nodal freedom matrix
   !* 
   !*    Dynamic real arrays
   !* 	 g_coord          : global nodal coordinates
   !*	 oldlds           : initial loads vector
   !*    etype            : element property type vector
   !*
   !*  OUTPUTS
   !*  AUTHOR
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen,nstep,npri
    INTEGER,   INTENT(IN)         :: g_num(:,:),nf(:,:)
    INTEGER                       :: i,j,k,l,m,n,nfe,nod,nels,ndim,nn
    INTEGER                       :: prnwidth,remainder
    REAL(iwp), INTENT(IN)         :: g_coord(:,:),loads(:),dtim
    REAL(iwp), INTENT(IN)         :: etype(:)
    CHARACTER(LEN=15), INTENT(IN) :: argv,element  
    CHARACTER(LEN=80)             :: cbuffer
    LOGICAL, INTENT(IN)           :: solid
    
!------------------------------------------------------------------------------
! 1. Initialisation
!------------------------------------------------------------------------------
  
    nn   = UBOUND(g_coord,2) ; ndim = UBOUND(g_coord,1)
    nels = UBOUND(g_num,2)   ; nod  = UBOUND(g_num,1)
  
!------------------------------------------------------------------------------
! 2. Write case file
!------------------------------------------------------------------------------
  
    OPEN(12,FILE=argv(1:nlen)//'.bin.ensi.case')
  
    WRITE(12,'(A/A)')    "#", "# Post-processing file generated by subroutine &
                               &WRITE_ENSI in "
    WRITE(12,'(A,A,/A)') "#"," Smith, Griffiths and Margetts, 'Programming the &
                               &Finite Element Method',","# Wiley, 2013."        
    WRITE(12,'(A/A/A)')  "#","# Ensight Gold Format","#"
    WRITE(12,'(2A/A)')   "# Problem name: ",argv(1:nlen),"#"
    WRITE(12,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
    WRITE(12,'(2A/A)')   "model: 1  ",argv(1:nlen)//'.bin.ensi.geo',"VARIABLE"
    WRITE(12,'(2A)')     "scalar per element:  material      ",                &
                          argv(1:nlen)//'.bin.ensi.MATID'
    IF(solid) THEN
      WRITE(12,'(2A)')   "scalar per node:     restraint     ",                &
                          argv(1:nlen)//'.bin.ensi.NDBND'
      WRITE(12,'(2A)')   "vector per node:     displacement  ",                &
                          argv(1:nlen)//'.bin.ensi.DISPL-******'
    ELSE
      WRITE(12,'(2A)')   "scalar per node:     pressure      ",                &
                          argv(1:nlen)//'.bin.ensi.PRESSURE-******'
    END IF
    WRITE(12,'(2A)')     "vector per node:     load          ",                &
                          argv(1:nlen)//'.bin.ensi.NDLDS'
    WRITE(12,'(A/A)')     "TIME","time set:     1"
    WRITE(12,'(A,I5)')    "number of steps:",nstep/npri
    WRITE(12,'(A,I5)')    "filename start number:",npri
    WRITE(12,'(A,I5)')    "filename increment:",npri
    WRITE(12,'(A)')       "time values:"
    prnwidth  = 5
    remainder = mod(nstep/npri,prnwidth)
    n         = ((nstep/npri) - remainder)/prnwidth
    IF(nstep/npri<=prnwidth) THEN
      DO i=1,nstep,npri
        IF(i==nstep) THEN
          WRITE(12,'(E12.5)') i*dtim
        ELSE
          WRITE(12,'(E12.5)',ADVANCE='no') i*dtim
        END IF
      END DO
    ELSE
      IF(remainder==0) THEN
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
      ELSE
  !     DO j=1,n-1
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
        m = (n*prnwidth)+1
        l = (n*prnwidth)+remainder
        DO i=m,l
          IF(i==l) THEN
            WRITE(12,'(E12.5)') dtim*i
          ELSE
            WRITE(12,'(E12.5)',ADVANCE='no') dtim*i
          END IF
        END DO
      END IF
    END IF
   
    CLOSE(12)
  
  !----------------------------------------------------------------------------
  ! 3. Write geometry file
  !
  !    Only 8 node bricks tested
  !----------------------------------------------------------------------------
  
    OPEN(13,FILE=argv(1:nlen)//'.bin.ensi.geo',STATUS="REPLACE",              &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "C Binary"                     ; WRITE(13) cbuffer
    cbuffer = "Problem name: "//argv(1:nlen) ; WRITE(13) cbuffer
    cbuffer = "Geometry files"               ; WRITE(13) cbuffer
    cbuffer = "node id off"                  ; WRITE(13) cbuffer
    cbuffer = "element id off"               ; WRITE(13) cbuffer
    cbuffer = "part"                         ; WRITE(13) cbuffer
    WRITE(13) int(1,kind=c_int)
    IF(ndim==2) THEN 
       cbuffer = "2d-mesh"                   ; WRITE(13) cbuffer
    END IF
    IF(ndim==3) THEN
       cbuffer = "Volume"                    ; WRITE(13) cbuffer
    END IF
    cbuffer = "coordinates"                  ; WRITE(13) cbuffer
    
    WRITE(13) int(nn,kind=c_int)
    DO j=1,ndim
      DO i=1,nn  
        WRITE(13) real(g_coord(j,i),kind=c_float)
      END DO
    END DO
  
    IF(ndim==2) THEN ! ensight requires zeros for the z-ordinate
      DO i=1,nn
        WRITE(13,'(A)') " 0.00000E+00" ! needs fixing for binary
      END DO
    END IF
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod)
!         CASE(3)
!           WRITE(13,'(A/I10)') "tria3", nels
!           DO i = 1,nels
!             WRITE(13,'(3I10)')g_num(3,i),g_num(2,i),g_num(1,i)
!           END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod)
          CASE(4)
            WRITE(13,'(A/I10)') "quad4", nels
            DO i = 1,nels
              WRITE(13,'(4I10)')g_num(1,i),g_num(4,i),g_num(3,i),g_num(2,i)
            END DO
          CASE(8)
            WRITE(13,'(A/I10)') "quad8", nels
            DO i = 1,nels
              WRITE(13,'(8I10)')g_num(1,i),g_num(7,i),g_num(5,i),g_num(3,i),   &
                                g_num(8,i),g_num(6,i),g_num(4,i),g_num(2,i)
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod)
          CASE(8)
            cbuffer = "hexa8"       ; WRITE(13) cbuffer
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13) int(g_num(1,i),kind=c_int),int(g_num(4,i),kind=c_int),&
                        int(g_num(8,i),kind=c_int),int(g_num(5,i),kind=c_int),&
                        int(g_num(2,i),kind=c_int),int(g_num(3,i),kind=c_int),&
                        int(g_num(7,i),kind=c_int),int(g_num(6,i),kind=c_int)
            END DO
          CASE(20)
            cbuffer = "hexa20"       ; WRITE(13) cbuffer
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13)                                                       &
                int(g_num(1,i),kind=c_int), int(g_num(7,i),kind=c_int),       &
                int(g_num(19,i),kind=c_int),int(g_num(13,i),kind=c_int),      &
                int(g_num(3,i),kind=c_int),int(g_num(5,i),kind=c_int),        &
                int(g_num(17,i),kind=c_int),int(g_num(15,i),kind=c_int),      &
                int(g_num(8,i),kind=c_int),int(g_num(12,i),kind=c_int),       &
                int(g_num(20,i),kind=c_int),int(g_num(9,i),kind=c_int),       &
                int(g_num(4,i),kind=c_int),int(g_num(11,i),kind=c_int),       &
                int(g_num(16,i),kind=c_int),int(g_num(10,i),kind=c_int),      &
                int(g_num(2,i),kind=c_int),int(g_num(6,i),kind=c_int),        &
                int(g_num(18,i),kind=c_int),int(g_num(14,i),kind=c_int) 
            END DO
          CASE DEFAULT
            cbuffer = "# Element type not recognised" ; WRITE(13) cbuffer
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            cbuffer = "tetra4" ; WRITE(13)
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13) int(g_num(1,i),kind=c_int),int(g_num(3,i),kind=c_int), &
                        int(g_num(2,i),kind=c_int),int(g_num(4,i),kind=c_int)
            END DO
          CASE DEFAULT
            cbuffer = "# Element type not recognised" ; WRITE(13)
        END SELECT
      CASE DEFAULT
        cbuffer = "# Element type not recognised" ; WRITE(13)
    END SELECT
  
    CLOSE(13)
  
  !-----------------------------------------------------------------------------
  ! 4. Write file containing material IDs
  !-----------------------------------------------------------------------------
  
    OPEN(14,FILE=argv(1:nlen)//'.bin.ensi.MATID',STATUS="REPLACE",             &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
    WRITE(14) cbuffer
    cbuffer = "part"
    WRITE(14) cbuffer
    WRITE(14) int(1,kind=c_int)
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod) 
          CASE(3)
            cbuffer = "tria3" ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod) 
          CASE(4)
            cbuffer = "quad4" ; WRITE(14) cbuffer
          CASE(8)
            cbuffer = "quad8" ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod) 
          CASE(8)
            cbuffer = "hexa8"   ; WRITE(14) cbuffer
          CASE(20)
            cbuffer = "hexa20"  ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            WRITE(14,'(A)') "tetra4"
          CASE DEFAULT
          WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE DEFAULT
        WRITE(14,'(A)')   "# Element type not recognised"
    END SELECT
   
    WRITE(14) real(etype(:),kind=c_float) 
  
    CLOSE(14)
  
  !-----------------------------------------------------------------------------
  ! 5. Write boundary conditions. Encoded using formula: 4z + 2y + 1x
  !
  !    110 = 1   010 = 2   100 = 3   011 = 4   101 = 5   001 = 6   000 = 7
  !-----------------------------------------------------------------------------
  
    IF(solid) THEN

      OPEN(15,FILE=argv(1:nlen)//'.bin.ensi.NDBND',STATUS="REPLACE",           &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

      cbuffer = "Alya Ensight Gold --- Scalar per-node variable file"
      WRITE(15)
      cbuffer = "part"         ; WRITE(15)
      WRITE(15) int(1,kind=c_int)
      cbuffer = "coordinates"  ; WRITE(15)

      IF(ndim==3) THEN
        DO i=1,UBOUND(g_coord,2) 
          nfe=0
          IF(nf(1,i)==0) nfe=nfe+1
          IF(nf(2,i)==0) nfe=nfe+2
          IF(nf(3,i)==0) nfe=nfe+4
          WRITE(15) int(nfe,kind=c_int)
        END DO
      ELSE IF(ndim==2) THEN
        DO i=1,nn
          nfe=0
          IF(nf(1,i)==0) nfe=nfe+1
          IF(nf(2,i)==0) nfe=nfe+2
          WRITE(15) int(nfe,kind=c_int)
        END DO
      ELSE
        PRINT *, "Wrong number of dimensions in mesh_ensi"
      END IF   
    END IF
  
    CLOSE(15)
  
  !-----------------------------------------------------------------------------
  ! 6. Write loaded nodes
  !-----------------------------------------------------------------------------
  
    OPEN(16,FILE=argv(1:nlen)//'.bin.ensi.NDLDS',STATUS="REPLACE",             &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "Alya Ensight Gold --- Vector per-node variable file"
    WRITE(16)
    cbuffer = "part"        ; WRITE(16)
    WRITE(16) int(1,kind=c_int)
    cbuffer = "coordinates" ; WRITE(16)

    DO j=1,UBOUND(nf,1)
      DO i=1, UBOUND(nf,2)
        WRITE(16) real(loads(nf(j,i)),kind=c_float)
      END DO
    END DO
    CLOSE(16)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI_BIN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

 SUBROUTINE MESH_ENSI_GEO_BIN(argv,nlen,g_coord,g_num,element)

   !/****f* input/mesh_ensi_geo_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_geo_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_geo_bin(argv,nlen,g_coord,g_num,element)
   !*  FUNCTION
   !*    This subroutine outputs the "geo" file in the C binary version of the 
   !*    Ensight gold format. Models in this format can be viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*	 element          : element type
   !*
   !*    Scalar logicals
   !*
   !*    Dynamic scalar arrays
   !*    g_num            : global element node numbers vector
   !* 
   !*    Dynamic real arrays
   !* 	 g_coord          : global nodal coordinates
   !*
   !*  OUTPUTS
   !*  AUTHOR
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !* 
   !*  Used in program p12meshgenbin
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen
    INTEGER,   INTENT(INOUT)      :: g_num(:,:)
    INTEGER                       :: i,j
    INTEGER                       :: nod,nels,ndim,nn
    REAL(iwp), INTENT(INOUT)      :: g_coord(:,:)
    CHARACTER(LEN=15), INTENT(IN) :: argv,element  
    CHARACTER(LEN=80)             :: cbuffer
    
  !----------------------------------------------------------------------------
  ! 1. Initialisation
  !----------------------------------------------------------------------------
  
    nn   = UBOUND(g_coord,2) ; ndim = UBOUND(g_coord,1)
    nels = UBOUND(g_num,2)   ; nod  = UBOUND(g_num,1)
  
  !----------------------------------------------------------------------------
  ! 2. Write geometry file
  !
  !    Only 8 node bricks tested
  !----------------------------------------------------------------------------
  
    OPEN(13,FILE=argv(1:nlen)//'.bin.ensi.geo',STATUS="REPLACE",              &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "C Binary"                     ; WRITE(13) cbuffer
    cbuffer = "Problem name: "//argv(1:nlen) ; WRITE(13) cbuffer
    cbuffer = "Geometry files"               ; WRITE(13) cbuffer
    cbuffer = "node id off"                  ; WRITE(13) cbuffer
    cbuffer = "element id off"               ; WRITE(13) cbuffer
    cbuffer = "part"                         ; WRITE(13) cbuffer
    WRITE(13) int(1,kind=c_int)
    IF(ndim==2) THEN 
       cbuffer = "2d-mesh"                   ; WRITE(13) cbuffer
    END IF
    IF(ndim==3) THEN
       cbuffer = "Volume"                    ; WRITE(13) cbuffer
    END IF
    cbuffer = "coordinates"                  ; WRITE(13) cbuffer
    
    WRITE(13) int(nn,kind=c_int)
    DO j=1,ndim
      DO i=1,nn  
        
        WRITE(13) real(g_coord(j,i),kind=c_float)
	
      END DO
    END DO


    IF(ndim==2) THEN ! ensight requires zeros for the z-ordinate
      DO i=1,nn
        WRITE(13,'(A)') " 0.00000E+00" ! needs fixing for binary
      END DO
    END IF
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod)
!         CASE(3)
!           WRITE(13,'(A/I10)') "tria3", nels
!           DO i = 1,nels
!             WRITE(13,'(3I10)')g_num(3,i),g_num(2,i),g_num(1,i)
!           END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod)
          CASE(4)
            WRITE(13,'(A/I10)') "quad4", nels
            DO i = 1,nels
              WRITE(13,'(4I10)')g_num(1,i),g_num(4,i),g_num(3,i),g_num(2,i)
            END DO
          CASE(8)
            WRITE(13,'(A/I10)') "quad8", nels
            DO i = 1,nels
              WRITE(13,'(8I10)')g_num(1,i),g_num(7,i),g_num(5,i),g_num(3,i),  &
                                g_num(8,i),g_num(6,i),g_num(4,i),g_num(2,i)
             

            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod)
          CASE(8)
            cbuffer = "hexa8"       ; WRITE(13) cbuffer
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13) int(g_num(1,i),kind=c_int),int(g_num(4,i),kind=c_int),&
                        int(g_num(8,i),kind=c_int),int(g_num(5,i),kind=c_int),&
                        int(g_num(2,i),kind=c_int),int(g_num(3,i),kind=c_int),&
                        int(g_num(7,i),kind=c_int),int(g_num(6,i),kind=c_int)
		
            END DO 
          CASE(20)
            cbuffer = "hexa20"       ; WRITE(13) cbuffer
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13)                                                       &
                int(g_num(1,i),kind=c_int), int(g_num(7,i),kind=c_int),       &
                int(g_num(19,i),kind=c_int),int(g_num(13,i),kind=c_int),      &
                int(g_num(3,i),kind=c_int),int(g_num(5,i),kind=c_int),        &
                int(g_num(17,i),kind=c_int),int(g_num(15,i),kind=c_int),      &
                int(g_num(8,i),kind=c_int),int(g_num(12,i),kind=c_int),       &
                int(g_num(20,i),kind=c_int),int(g_num(9,i),kind=c_int),       &
                int(g_num(4,i),kind=c_int),int(g_num(11,i),kind=c_int),       &
                int(g_num(16,i),kind=c_int),int(g_num(10,i),kind=c_int),      &
                int(g_num(2,i),kind=c_int),int(g_num(6,i),kind=c_int),        &
                int(g_num(18,i),kind=c_int),int(g_num(14,i),kind=c_int) 
	
            END DO		
          CASE DEFAULT
            cbuffer = "# Element type not recognised" ; WRITE(13) cbuffer
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            cbuffer = "tetra4" ; WRITE(13)
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13) int(g_num(1,i),kind=c_int),int(g_num(3,i),kind=c_int),&
                        int(g_num(2,i),kind=c_int),int(g_num(4,i),kind=c_int)
            END DO
          CASE DEFAULT
            cbuffer = "# Element type not recognised" ; WRITE(13)
        END SELECT
      CASE DEFAULT
        cbuffer = "# Element type not recognised" ; WRITE(13)
    END SELECT
  
    CLOSE(13)

    RETURN
  
  END SUBROUTINE MESH_ENSI_GEO_BIN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

 SUBROUTINE MESH_ENSI_MATID_BIN(argv,nlen,nod,element,etype)

   !/****f* input/mesh_ensi_matid_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_matid_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_matid_bin(argv,nlen,nod,element,etype)
   !*  FUNCTION
   !*    This subroutine outputs material type for each element in the mesh
   !*    in the C binary version of the Ensight gold format. Models in this 
   !*    format can be viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    nod              : number of nodes in the element
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*	 element          : element type
   !*
   !*    Dynamic scalar arrays
   !*    etype            : element property type vector
   !* 
   !*  OUTPUTS
   !*  AUTHOR
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*  ParaView has a bug which prevents the 'Material' section of the 
   !*  ENSIGHT Gold to be read and therefore integer MATID values
   !*  http://www.paraview.org/Bug/view.php?id=15151
   !*  http://www3.ensight.com/EnSight10_Docs/UserManual.pdf pp.713
   !*  Workaround - use MATIDs as reals and convert into int for ParaFEM
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen,nod
    REAL(iwp), INTENT(IN)         :: etype(:)
    CHARACTER(LEN=15), INTENT(IN) :: argv,element  
    CHARACTER(LEN=80)             :: cbuffer
  
!------------------------------------------------------------------------------
! 1. Write file containing material IDs
!------------------------------------------------------------------------------
  
    OPEN(14,FILE=argv(1:nlen)//'.bin.ensi.MATID',STATUS="REPLACE",            &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
    WRITE(14) cbuffer
    cbuffer = "part" ;  WRITE(14)
    WRITE(14) cbuffer
    WRITE(14) int(1,kind=c_int)
    cbuffer = "coordinates" ; WRITE(14)
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod) 
          CASE(3)
            cbuffer = "tria3" ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod) 
          CASE(4)
            cbuffer = "quad4" ; WRITE(14) cbuffer
          CASE(8)
            cbuffer = "quad8" ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod) 
          CASE(8)
            cbuffer = "hexa8"   ; WRITE(14) cbuffer
          CASE(20)
            cbuffer = "hexa20"  ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            WRITE(14,'(A)') "tetra4"
          CASE DEFAULT
          WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE DEFAULT
        WRITE(14,'(A)')   "# Element type not recognised"
    END SELECT
   
    WRITE(14) real(etype(:),kind=c_float) 
    
    CLOSE(14)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI_MATID_BIN
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

 SUBROUTINE MESH_ENSI_NDBND_BIN(argv,nf,nlen,nod,solid)

   !/****f* input/mesh_ensi_ndbnd_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_ndbnd_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_ndbnd_bin(argv,nf,nlen,nod,solid)
   !*  FUNCTION
   !*    This subroutine outputs a file of restrained nodes in the C binary 
   !*    version of the Ensight gold format. Models in this format can be 
   !*    viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    nod              : number of nodes per element
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*
   !*    Scalar logicals
   !*    solid            : type of analysis solid if .true. fluid if .false.
   !*
   !*    Dynamic scalar arrays
   !*    nf               : nodal freedom matrix
   !*  OUTPUTS
   !*    File: <job_name>.bin.ensi.NDBND 
   !*  AUTHOR
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen,nod
    INTEGER,   INTENT(IN)         :: nf(:,:)
    INTEGER                       :: i,nn,ndim
!   INTEGER                       :: nfe
    REAL(iwp)                     :: nfe
    CHARACTER(LEN=15), INTENT(IN) :: argv  
    CHARACTER(LEN=80)             :: cbuffer
    LOGICAL, INTENT(IN)           :: solid
    
  !-----------------------------------------------------------------------------
  ! 1. Initialisation
  !-----------------------------------------------------------------------------
  
    ndim = UBOUND(nf,1)-1  
    nn   = UBOUND(nf,2)
 
    PRINT *, "ndim=",ndim
    PRINT *, "nn=",nn
    PRINT *, "solid=",solid
 
  !-----------------------------------------------------------------------------
  ! 2. Write boundary conditions. Encoded using formula: 4z + 2y + 1x
  !
  !    110 = 1   010 = 2   100 = 3   011 = 4   101 = 5   001 = 6   000 = 7
  !-----------------------------------------------------------------------------
  
    IF(solid) THEN

      OPEN(15,FILE=argv(1:nlen)//'.bin.ensi.NDBND',STATUS="REPLACE",           &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

      cbuffer = "Alya Ensight Gold --- Scalar per-node variable file"
      WRITE(15) cbuffer
      cbuffer = "part"         ; WRITE(15) cbuffer
      WRITE(15) int(1,kind=c_int)
      cbuffer = "coordinates"  ; WRITE(15) cbuffer

      IF(ndim==3) THEN
        DO i=1,nn 
          nfe=0
          IF(nf(2,i)==0) nfe=nfe+1
          IF(nf(3,i)==0) nfe=nfe+2
          IF(nf(4,i)==0) nfe=nfe+4
          WRITE(15) real(nfe,kind=c_float) 
         !WRITE(15) int(nfe,kind=c_int)
        END DO
      ELSE IF(ndim==2) THEN
        DO i=1,nn
          nfe=0
          IF(nf(1,i)==0) nfe=nfe+1
          IF(nf(2,i)==0) nfe=nfe+2
          WRITE(15) int(nfe,kind=c_int)
        END DO
      ELSE
        PRINT *, "Wrong number of dimensions in mesh_ensi"
      END IF   
    END IF
  
	
    CLOSE(15)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI_NDBND_BIN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE MESH_ENSI_NDLDS_BIN(argv,nlen,nn,val,node)

   !/****f* input/mesh_ensi_ndlds_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_ndlds_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_ndlds_bin(argv,nlen,nf,loads)
   !*  FUNCTION
   !*    This subroutine outputs a file of loads in the C binary version of the 
   !*    Ensight gold format. Models in this format can be viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*
   !*    Dynamic scalar arrays
   !*    nf               : nodal freedom matrix
   !* 
   !*    Dynamic real arrays
   !*	 xxxxx           : initial loads vector
   !*
   !*  OUTPUTS
   !*  AUTHORS
   !*    J.D Arregui-Mena &L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2016
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER, INTENT(IN)           :: nlen,nn
    INTEGER, INTENT(IN)           :: node(:)  !(loaded_nodes)
    REAL(iwp), INTENT(IN)         :: val(:,:) !(ndim,loaded_nodes)
    INTEGER                       :: i,j,k,ndim
    CHARACTER(LEN=15), INTENT(IN) :: argv
    CHARACTER(LEN=80)             :: cbuffer
    REAL(iwp),PARAMETER           :: zero=0.0_iwp

    ndim = UBOUND(val,1)  

  !-----------------------------------------------------------------------------
  ! 1. Write loaded nodes
  !-----------------------------------------------------------------------------
  
    OPEN(16,FILE=argv(1:nlen)//'.bin.ensi.NDLDS',STATUS="REPLACE",             &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "Alya Ensight Gold --- Vector per-node variable file"
    WRITE(16) cbuffer
    cbuffer = "part"        
    WRITE(16) cbuffer
    WRITE(16) int(1,kind=c_int)
    cbuffer = "coordinates" 
    WRITE(16) cbuffer

   
    DO i=1,ndim
     k=1
      DO j=1,nn
        IF(j==node(k)) THEN
          k=k+1
          WRITE(16) REAL(val(i,k), kind=c_float)
        ELSE
          WRITE(16) REAL(zero, kind=c_float)
        END IF
      END DO
    END DO
       
    CLOSE(16)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI_NDLDS_BIN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MESH_ENSI_CASE(argv,nlen,nstep,npri,dtim,solid)

   !/****f* input/mesh_ensi_case
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_case
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_case(argv,nlen,nstep,npri,dtim,solid)
   !*  FUNCTION
   !*    This subroutine outputs the "case" file required for the C binary 
   !*    version of the Ensight gold format. Models in this format can be 
   !*    viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    npri             : print interval
   !*	 nstep            : number of time steps in analysis
   !*
   !*    Scalar reals
   !*    dtim             : time step
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*
   !*    Scalar logicals
   !*    solid            : type of analysis solid if .true. fluid if .false.
   !*
   !*  OUTPUTS
   !*  AUTHOR
   !*    I.M. Smith
   !*    D.V. Griffiths
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen,nstep,npri
    INTEGER                       :: i,j,k,l,m,n,nfe,nod,nels,ndim,nn
    INTEGER                       :: prnwidth,remainder
    REAL(iwp), INTENT(IN)         :: dtim
    CHARACTER(LEN=15), INTENT(IN) :: argv
    CHARACTER(LEN=80)             :: cbuffer
    LOGICAL, INTENT(IN)           :: solid
    
  !-----------------------------------------------------------------------------
  ! 1. Write case file
  !-----------------------------------------------------------------------------
  
    OPEN(12,FILE=argv(1:nlen)//'.bin.ensi.case')
  
    WRITE(12,'(A/A)')    "#", "# Post-processing file generated by subroutine  &
                               &WRITE_ENSI in "
    WRITE(12,'(A,A,/A)') "#"," Smith, Griffiths and Margetts, 'Programming the &
                               &Finite Element Method',","# Wiley, 2013."        
    WRITE(12,'(A/A/A)')  "#","# Ensight Gold Format","#"
    WRITE(12,'(2A/A)')   "# Problem name: ",argv(1:nlen),"#"
    WRITE(12,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
    WRITE(12,'(2A/A)')   "model: 1  ",argv(1:nlen)//'.bin.ensi.geo',"VARIABLE"
    WRITE(12,'(2A)')     "scalar per element:  material      ",                &
                          argv(1:nlen)//'.bin.ensi.MATID'
    IF(solid) THEN
      WRITE(12,'(2A)')   "scalar per node:     restraint     ",                &
                          argv(1:nlen)//'.bin.ensi.NDBND'
      WRITE(12,'(2A)')   "vector per node:     displacement  ",                &
                          argv(1:nlen)//'.bin.ensi.DISPL-******'
    ELSE
      WRITE(12,'(2A)')   "scalar per node:     pressure      ",                &
                          argv(1:nlen)//'.bin.ensi.PRESSURE-******'
    END IF
    WRITE(12,'(2A)')     "vector per node:     load          ",                &
                          argv(1:nlen)//'.bin.ensi.NDLDS'
    WRITE(12,'(A/A)')     "TIME","time set:     1"
    WRITE(12,'(A,I5)')    "number of steps:",nstep/npri
    WRITE(12,'(A,I5)')    "filename start number:",npri
    WRITE(12,'(A,I5)')    "filename increment:",npri
    WRITE(12,'(A)')       "time values:"
    prnwidth  = 5
    remainder = mod(nstep/npri,prnwidth)
    n         = ((nstep/npri) - remainder)/prnwidth
    IF(nstep/npri<=prnwidth) THEN
      DO i=1,nstep,npri
        IF(i==nstep) THEN
          WRITE(12,'(E12.5)') i*dtim
        ELSE
          WRITE(12,'(E12.5)',ADVANCE='no') i*dtim
        END IF
      END DO
    ELSE
      IF(remainder==0) THEN
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
      ELSE
  !     DO j=1,n-1
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
        m = (n*prnwidth)+1
        l = (n*prnwidth)+remainder
        DO i=m,l
          IF(i==l) THEN
            WRITE(12,'(E12.5)') dtim*i
          ELSE
            WRITE(12,'(E12.5)',ADVANCE='no') dtim*i
          END IF
        END DO
      END IF
    END IF
   
    CLOSE(12)
    
    RETURN
  
  END SUBROUTINE MESH_ENSI_CASE
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MESH_ENSI_GEO(argv,nlen,g_coord,g_num,element)

   !/****f* input/mesh_ensi_geo
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_geo
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_geo(argv,nlen,g_coord,g_num,element)
   !*  FUNCTION
   !*    This subroutine outputs a set of files in the C binary version of the 
   !*    Ensight gold format. Models in this format can be viewed in ParaView.
   !* 
   !*    Element types supported:                Tested with:
   !*
   !*    2-node bar
   !*    3-node triangle                         p51  (4th edition p51_1.dat)
   !*    6-node triangle
   !*    4-node quadrilateral                    p115 (4th edition)
   !*    8-node quadrilateral                    p116 (4th edition)
   !*    4-node tetrahedron                      p54  (4th edition p54_2.dat)
   !*    8-node hexahedron                       p86  (4th edition)       
   !*    20-node hexahedron                      p55  (4th edition)
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    npri             : print interval
   !*	 nstep            : number of time steps in analysis
   !*
   !*    Scalar reals
   !*    dtim             : time step
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*	 element          : element type
   !*
   !*    Scalar logicals
   !*    solid            : type of analysis solid if .true. fluid if .false.
   !*
   !*    Dynamic scalar arrays
   !*    g_num            : global element node numbers vector
   !* 
   !*    Dynamic real arrays
   !* 	 g_coord          : global nodal coordinates
   !*
   !*  OUTPUTS
   !*  AUTHOR
   !*    I.M. Smith
   !*    D.V. Griffiths
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen
    INTEGER,   INTENT(IN)         :: g_num(:,:)
    INTEGER                       :: i,j,k,l,m,n,nod,nels,ndim,nn
    REAL(iwp), INTENT(IN)         :: g_coord(:,:)
    CHARACTER(LEN=15), INTENT(IN) :: argv,element  
    CHARACTER(LEN=80)             :: cbuffer
    
  !----------------------------------------------------------------------------
  ! 1. Initialisation
  !----------------------------------------------------------------------------
  
    nn   = UBOUND(g_coord,2) ; ndim = UBOUND(g_coord,1)
    nels = UBOUND(g_num,2)   ; nod  = UBOUND(g_num,1)
  
  !----------------------------------------------------------------------------
  ! 2. Write geometry file
  !
  !    Only 8 node bricks tested
  !----------------------------------------------------------------------------
  
    OPEN(13,FILE=argv(1:nlen)//'.bin.ensi.geo',STATUS="REPLACE",              &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "C Binary"                     ; WRITE(13) cbuffer
    cbuffer = "Problem name: "//argv(1:nlen) ; WRITE(13) cbuffer
    cbuffer = "Geometry files"               ; WRITE(13) cbuffer
    cbuffer = "node id off"                  ; WRITE(13) cbuffer
    cbuffer = "element id off"               ; WRITE(13) cbuffer
    cbuffer = "part"                         ; WRITE(13) cbuffer
    WRITE(13) int(1,kind=c_int)
    IF(ndim==2) THEN 
       cbuffer = "2d-mesh"                   ; WRITE(13) cbuffer
    END IF
    IF(ndim==3) THEN
       cbuffer = "Volume"                    ; WRITE(13) cbuffer
    END IF
    cbuffer = "coordinates"                  ; WRITE(13) cbuffer
    
    WRITE(13) int(nn,kind=c_int)
    DO j=1,ndim
      DO i=1,nn  
        WRITE(13) real(g_coord(j,i),kind=c_float)
      END DO
    END DO
  
    IF(ndim==2) THEN ! ensight requires zeros for the z-ordinate
      DO i=1,nn
        WRITE(13,'(A)') " 0.00000E+00" ! needs fixing for binary
      END DO
    END IF
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod)
!         CASE(3)
!           WRITE(13,'(A/I10)') "tria3", nels
!           DO i = 1,nels
!             WRITE(13,'(3I10)')g_num(3,i),g_num(2,i),g_num(1,i)
!           END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod)
          CASE(4)
            WRITE(13,'(A/I10)') "quad4", nels
            DO i = 1,nels
              WRITE(13,'(4I10)')g_num(1,i),g_num(4,i),g_num(3,i),g_num(2,i)
            END DO
          CASE(8)
            WRITE(13,'(A/I10)') "quad8", nels
            DO i = 1,nels
              WRITE(13,'(8I10)')g_num(1,i),g_num(7,i),g_num(5,i),g_num(3,i),   &
                                g_num(8,i),g_num(6,i),g_num(4,i),g_num(2,i)
            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod)
          CASE(8)
            cbuffer = "hexa8"       ; WRITE(13) cbuffer
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13) int(g_num(1,i),kind=c_int),int(g_num(4,i),kind=c_int),&
                        int(g_num(8,i),kind=c_int),int(g_num(5,i),kind=c_int),&
                        int(g_num(2,i),kind=c_int),int(g_num(3,i),kind=c_int),&
                        int(g_num(7,i),kind=c_int),int(g_num(6,i),kind=c_int)
            END DO
          CASE(20)
            cbuffer = "hexa20"       ; WRITE(13) cbuffer
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13)                                                       &
                int(g_num(1,i),kind=c_int), int(g_num(7,i),kind=c_int),       &
                int(g_num(19,i),kind=c_int),int(g_num(13,i),kind=c_int),      &
                int(g_num(3,i),kind=c_int),int(g_num(5,i),kind=c_int),        &
                int(g_num(17,i),kind=c_int),int(g_num(15,i),kind=c_int),      &
                int(g_num(8,i),kind=c_int),int(g_num(12,i),kind=c_int),       &
                int(g_num(20,i),kind=c_int),int(g_num(9,i),kind=c_int),       &
                int(g_num(4,i),kind=c_int),int(g_num(11,i),kind=c_int),       &
                int(g_num(16,i),kind=c_int),int(g_num(10,i),kind=c_int),      &
                int(g_num(2,i),kind=c_int),int(g_num(6,i),kind=c_int),        &
                int(g_num(18,i),kind=c_int),int(g_num(14,i),kind=c_int) 
            END DO
          CASE DEFAULT
            cbuffer = "# Element type not recognised" ; WRITE(13) cbuffer
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            cbuffer = "tetra4" ; WRITE(13)
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13) int(g_num(1,i),kind=c_int),int(g_num(3,i),kind=c_int), &
                        int(g_num(2,i),kind=c_int),int(g_num(4,i),kind=c_int)
            END DO
          CASE DEFAULT
            cbuffer = "# Element type not recognised" ; WRITE(13)
        END SELECT
      CASE DEFAULT
        cbuffer = "# Element type not recognised" ; WRITE(13)
    END SELECT
  
    CLOSE(13)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI_GEO
    
END MODULE INPUT
