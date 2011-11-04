MODULE GATHER_SCATTER

  !/****h* /gather_scatter
  !*  NAME
  !*    MODULE: gather_scatter
  !*  SYNOPSIS
  !*    Usage:      USE gather_scatter
  !*  FUNCTION
  !*    
  !*    Subroutine                 Purpose
  !*
  !*    MY_BARRIER                 Synchronises all processors
  !*    MPERROR                    Returns an error code
  !*    CALC_NELS_PP               Computes number of elements per processor
  !*    CALC_NEQ_PP                Computes number of equations per processor
  !*    CALC_NPES_PP               Computes size of arrays for data exchange
  !*    ALLOCATE_GATHER_SCATTER    Allocates arrays used in GATHER & SCATTER
  !*    DEALLOCATE_GATHER_SCATTER  Deallocates arrays used in GATHER & SCATTER
  !*    GATHER                     Performs the gather operation: pmul = p(g)
  !*    SCATTER                    Performs the operation: u(g) = u(g) + utemp
  !*    SCATTER_NOADD              Performs the operation: u(g) = utemp
  !*    MAKE_GGL                   Relates equations to elements and processors
  !*    SCATTER_NODES                         
  !*
  !*  AUTHOR
  !*    M.A. Pettipher
  !*    L. Margetts
  !*    V. Szeremi
  !*    F. Calvo
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  USE precision
  USE mp_interface
  USE global_variables
  USE input, ONLY: read_nels_pp

!------------------------------------------------------------------------------
! 1. Variables restricted to gather_scatter module:
!------------------------------------------------------------------------------
  
  INTEGER, ALLOCATABLE            :: ggl_pp(:,:)
  INTEGER, ALLOCATABLE, PRIVATE   :: lenget(:), lenput(:), lengetsum(:)
  INTEGER, ALLOCATABLE, PRIVATE   :: toget(:,:), toput(:,:)
  INTEGER, ALLOCATABLE, PRIVATE   :: pesget(:), pesput(:)
  INTEGER, ALLOCATABLE, PRIVATE   :: getpes(:), putpes(:)
  INTEGER, ALLOCATABLE, PRIVATE   :: toget_temp(:,:), toput_temp(:,:)
  REAL(iwp), ALLOCATABLE, PRIVATE :: tempget(:), tempput(:), tempput1(:,:)
  REAL(iwp), ALLOCATABLE, PRIVATE :: sumget(:,:)
  INTEGER, PRIVATE                :: numpesget, numpesput, numpesgetput
  INTEGER, PRIVATE                :: len_pl_pp

!------------------------------------------------------------------------------
! 2. Variables also used in main program:
!------------------------------------------------------------------------------

  INTEGER :: npes, neq_pp1, neq_pp2
  INTEGER :: num_neq_pp1, threshold
  INTEGER :: ibar
  INTEGER :: iel_start, ieq_start

  CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MY_BARRIER(numpe,ibar,channel,chstr)

  !/****f* gather_scatter/my_barrier
  !*  NAME
  !*    SUBROUTINE: my_barrier
  !*  SYNOPSIS
  !*    Usage:      CALL my_barrier(numpe,ibar,channel,chstr)
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !*/ 

    IMPLICIT NONE
  
    INTEGER, INTENT(IN)           :: numpe, ibar, channel
    CHARACTER (LEN=*), INTENT(IN) :: chstr
    INTEGER                       :: ier

    ier = 0
  
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  
    IF (ier .NE. MPI_SUCCESS) THEN
      CALL MPERROR(chstr,ier)
    END IF

  END SUBROUTINE MY_BARRIER

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE MPERROR(cherror,errcode)

  !/****f* gather_scatter/mperror
  !*  NAME
  !*    SUBROUTINE: mperror
  !*  SYNOPSIS
  !*    Usage:      CALL mperror(cherror,errcode)
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !* Should call SHUTDOWN rather than STOP. Best practice is to use STOP
  !* in one subroutine only.
  !*
  !*/ 

    IMPLICIT NONE
  
    CHARACTER*(*) :: cherror
    INTEGER       :: errcode
  
    WRITE(*,'(A,A,I5)') cherror, ' errcode = ', errcode
  
    STOP

  END SUBROUTINE MPERROR

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE CALC_NELS_PP(job_name,nels,npes,numpe,partitioner,nels_pp)

  !/****f* gather_scatter/calc_nels_pp
  !*  NAME
  !*    SUBROUTINE: calc_nels_pp
  !*  SYNOPSIS
  !*    Usage:      CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner, &
  !*                                  nels_pp)
  !*  FUNCTION
  !*    Calculates the number of elements, NELS_PP, assigned to each processor.
  !*    It is effectively a very naive method of mesh partitioning. The 
  !*    subroutine also computes, IEL_START, the first element number on each 
  !*    processor. IEL_START and NELS_PP are the external variables modified by
  !*    this subroutine. Note that they are global variables that are not
  !*    passed through the list of arguments.
  !*
  !*    An example follows which explains how this subroutine works: 
  !*
  !*    If we have:   nels = 103 (number of elements)     
  !*                  npes = 5  (number of processors)
  !*
  !*    The output will be:
  !*
  !*    nels_pp2     = 20 (103/5=20, so 20 elements for each processor and the 
  !*                       remaining 3 elements will be assigned to the first 
  !*                       3 processors)
  !*
  !*    num_nels_pp1 = 3  (103 - 20*5 = 3, these are the 3 remaining elements)
  !*
  !*    nels_pp2     = 21 (no. of elements assigned to the first 3 processors)
  !*
  !*    So the result is:
  !*
  !*    On processor 1 (numpe = 1) -----> nels_pp=21  iel_start=1
  !*    On processor 2 (numpe = 2) -----> nels_pp=21  iel_start=22
  !*    On processor 3 (numpe = 3) -----> nels_pp=21  iel_start=43
  !*    On processor 4 (numpe = 4) -----> nels_pp=20  iel_start=64
  !*    On processor 5 (numpe = 5) -----> nels_pp=20  iel_start=84
  !*
  !*  INPUTS
  !*    The following argument has the INTENT(IN) attribute:
  !*
  !*    nels                  : Integer
  !*                          : Total number of elements in the mesh
  !*  AUTHOR
  !*    M.A. Pettipher
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2011
  !******
  !*
  !* Would be improved if nels_pp and iel_start were subroutine arguments
  !* 
  !* 20.05.2011 Added ability to use external partitioner
  !*/
  
    IMPLICIT NONE
  
    CHARACTER(LEN=50),INTENT(IN) :: job_name 
    INTEGER, INTENT(IN)          :: nels,npes,numpe,partitioner
    INTEGER, INTENT(INOUT)       :: nels_pp
    INTEGER                      :: num_nels_pp1, nels_pp1, nels_pp2

    num_nels_pp1 = 0
    nels_pp1     = 0
    nels_pp2     = 0
 
    SELECT CASE (partitioner)
 
      CASE (1) ! S&G partitioning 

      IF (npes == 1) THEN
        nels_pp1 = 0
        nels_pp2 = nels
        nels_pp  = nels_pp2
        iel_start = 1
      ELSE
        nels_pp2     = nels/npes
        num_nels_pp1 = nels - nels_pp2*npes
        IF (num_nels_pp1 == 0) THEN
          nels_pp1   = nels_pp2
        ELSE
          nels_pp1   = nels_pp2 + 1
        ENDIF
        IF (numpe <= num_nels_pp1 .OR. num_nels_pp1 == 0) THEN
          nels_pp    = nels_pp1
          iel_start  = (numpe - 1)*nels_pp1 + 1
        ELSE
          nels_pp    = nels_pp2
          iel_start  = num_nels_pp1*nels_pp1+                                 &
                       (numpe-num_nels_pp1-1)*(nels_pp1-1)+1
        ENDIF
      ENDIF
  
    CASE (2)

      CALL read_nels_pp(job_name,nels_pp,npes,numpe)

    CASE DEFAULT
     
      PRINT*
      PRINT*, "Incorrect value for variable PARTITIONER. Permitted values are:"
      PRINT*, "  1 - Smith and Griffiths partitioning"
      PRINT*, "  2 - External partitioning (requires .psize file)"

      CALL shutdown()

    END SELECT

    WRITE(details,'(A,I8,A,I8)') 'PE no: ', numpe, ' nels_pp: ', nels_pp

  END SUBROUTINE CALC_NELS_PP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE CALC_NEQ_PP
  
  !/****f* gather_scatter/calc_neq_pp
  !*  NAME
  !*    SUBROUTINE: calc_neq_pp
  !*  SYNOPSIS
  !*    Usage:      CALL calc_neq_pp(nels)
  !*  FUNCTION
  !*    Calculates the number of equations, NEQ_PP, assigned to each processor.
  !*    It is effectively a very naive method of mesh partitioning. The 
  !*    subroutine also computes, IEQ_START, the first equation number on each 
  !*    processor. IEQ_START and NEQ_PP are external variables modified by this
  !*    subroutine. Note that they are global variables that are not passed
  !*    through the list of arguments.
  !*
  !*    An example follows which explains how this subroutine works: 
  !*
  !*    If we have:   nels = 103 (number of equations)     
  !*                  npes = 5   (number of processors)
  !*
  !*    The output will be:
  !*
  !*    neq_pp2      = 20 (103/5=20, so 20 equations for each processor and the 
  !*                       remaining 3 equations will be assigned to the first 
  !*                       3 processors)
  !*
  !*    num_neq_pp1  = 3  (103 - 20*5 = 3, these are the 3 remaining equations)
  !*
  !*    nels_pp2     = 21 (number of elements assigned to first 3 processors)
  !*
  !*    So the result is:
  !*
  !*    On processor 1 (numpe = 1) -----> neq_pp=21    ieq_start=1
  !*    On processor 2 (numpe = 2) -----> neq_pp=21    ieq_start=22
  !*    On processor 3 (numpe = 3) -----> neq_pp=21    ieq_start=43
  !*    On processor 4 (numpe = 4) -----> neq_pp=20    ieq_start=64
  !*    On processor 5 (numpe = 5) -----> neq_pp=20    ieq_start=84
  !*
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !* Would be improved if neq, neq_pp and ieq_start were subroutine arguments
  !*
  !*/
  
    IMPLICIT NONE

    neq_pp2     = 0
    neq_pp1     = 0
    neq_pp      = 0
    ieq_start   = 0
    num_neq_pp1 = 0

    IF (npes == 1) THEN
      neq_pp2   = neq
      neq_pp1   = neq
      neq_pp    = neq_pp1
      ieq_start = 1
    ELSE
      neq_pp2     = neq/npes
      num_neq_pp1 = neq - neq_pp2*npes
      IF (num_neq_pp1 == 0) THEN
        neq_pp1   = neq_pp2
      ELSE
        neq_pp1   = neq_pp2 + 1
      ENDIF
      IF (numpe <= num_neq_pp1 .OR. num_neq_pp1 == 0) THEN
        neq_pp    = neq_pp1
        ieq_start = (numpe - 1)*neq_pp1 + 1
      ELSE
        neq_pp    = neq_pp2
        ieq_start = num_neq_pp1*neq_pp1 + (numpe-num_neq_pp1-1)*(neq_pp1-1) + 1
      ENDIF
    ENDIF
  
    WRITE(details,'(A,I8,A,I8)') 'PE no: ', numpe, ' neq_pp: ', neq_pp

  END SUBROUTINE CALC_NEQ_PP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE CALC_NPES_PP(npes,npes_pp)
  
  !/****f* gather_scatter/calc_npes_pp
  !*  NAME
  !*    SUBROUTINE: calc_npes_pp
  !*  SYNOPSIS
  !*    Usage:      CALL calc_npes_pp(npes,npes_pp)
  !*  FUNCTION
  !*    Sets value to be used in initial dimensions of toget and toput. We do
  !*    not know the exact number of processors needed by each processor until
  !*    MAKE_GGL is called, but we must declare these arrays (or at least
  !*    temporary versions of them) before the exact number is known. NPES
  !*    is definitely enough, but wasteful of memory, so make rough
  !*    overestimate based on number of processors, NPES.
  !*  AUTHOR
  !*    M.A. Pettipher
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !* This subroutine needs to be rewritten. It causes the most execution 
  !* failures, particularly for pathological cases.
  !*/ 
 
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: npes_pp
    INTEGER, INTENT(IN)  :: npes
 
    SELECT CASE (npes)

      CASE (1:15)
        npes_pp = npes
      CASE (16:32)  
        npes_pp = npes/2
      CASE (33:256)
        npes_pp = npes/4
      CASE (257:1024)
        npes_pp = npes/7
      CASE DEFAULT
        npes_pp = npes/12

    END SELECT

  END SUBROUTINE CALC_NPES_PP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE ALLOCATE_GATHER_SCATTER(npes_pp,npes)

  !/****f* gather_scatter/allocate_gather_scatter
  !*  NAME
  !*    SUBROUTINE: allocate_gather_scatter
  !*  SYNOPSIS
  !*    Usage:      CALL allocate_gather_scatter(npes_pp,npes)
  !*  FUNCTION
  !*    Allocates the arrays used in the subroutines GATHER and SCATTER.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    npes_pp               : Integer
  !*                          : Number of processors each processor needs to 
  !*                          : communicate with.
  !*
  !*    npes                  : Integer
  !*                          : Total number of processors used by the program.
  !*                         
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !* 
  !* This subroutine could be improved by computing the actual storage
  !* requirements
  !*
  !*/ 
 
    INTEGER, INTENT(IN) :: npes_pp,npes

    ALLOCATE(ggl_pp(ntot,nels_pp))
    ALLOCATE(lenget(npes), lenput(npes), lengetsum(0:npes))
    ALLOCATE(pesget(npes), pesput(npes))
    ALLOCATE(getpes(npes), putpes(npes))

    ggl_pp = 0 ; lenget = 0; lenput = 0 ; lengetsum = 0
    pesget = 0 ; pesput = 0; getpes = 0 ; putpes    = 0

!------------------------------------------------------------------------------
! 1. Not obvious what are reasonable values for array dimensions. With 1000 
!    element problem and 32 PEs, neq_pp1 = 413, but PE 31 requires 832 remote 
!    accesses, so array size must be more than 2*neq_pp1 - so set to 3*neq_pp1.
!------------------------------------------------------------------------------

    ALLOCATE(tempget(3*neq_pp1), tempput(3*neq_pp1), tempput1(3*neq_pp1,npes))

    tempget = 0.0_iwp ; tempput = 0.0_iwp ; tempput1 = 0.0_iwp

!------------------------------------------------------------------------------
! 2. Allocate temporary arrays for toget and toput, until actual size is known. 
!    Reallocate later.
!------------------------------------------------------------------------------

    ALLOCATE(toget_temp(neq_pp1,npes_pp), toput_temp(neq_pp1,npes_pp))
    ALLOCATE(sumget(neq_pp1,npes_pp))

    toget_temp = 0 ; toput_temp = 0 ; sumget = 0.0_iwp

  END SUBROUTINE ALLOCATE_GATHER_SCATTER

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE DEALLOCATE_GATHER_SCATTER

  !/****f* gather_scatter/deallocate_gather_scatter
  !*  NAME
  !*    SUBROUTINE: deallocate_gather_scatter
  !*  SYNOPSIS
  !*    Usage:      CALL deallocate_gather_scatter(npes_pp,npes)
  !*  FUNCTION
  !*    Deallocates the arrays used in the subroutines GATHER and SCATTER.
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !*/ 

    DEALLOCATE(ggl_pp)
    DEALLOCATE(lenget, lenput, lengetsum)
    DEALLOCATE(pesget, pesput)
    DEALLOCATE(getpes, putpes)
    DEALLOCATE(tempget, tempput, tempput1)

  END SUBROUTINE DEALLOCATE_GATHER_SCATTER

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE GATHER(p_pp,pmul_pp)

  !/****f* gather_scatter/gather
  !*  NAME
  !*    SUBROUTINE: gather
  !*  SYNOPSIS
  !*    Usage:      CALL gather(p_pp,pmul_pp)
  !*  FUNCTION
  !*    Performs the gather operation: pmul = p(g).
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    p_pp(:)               : Real
  !*                          : Distributed vector of dimension NEQ_PP
  !*
  !*  OUTPUTS
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    pmul_pp(:,:)          : Real
  !*                          : Distributed array(NTOT,NELS_PP) that is 
  !*                          : to be populated with the values in P_PP, 
  !*                          : arranged in an element-by-element form.
  !*                  
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !*/ 

    IMPLICIT NONE
   
    REAL(iwp), INTENT(IN)  :: p_pp(neq_pp)
    REAL(iwp), INTENT(OUT) :: pmul_pp(ntot,nels_pp)
    LOGICAL                :: lflag
    INTEGER                :: i, ii, j, k, ier, nstart, pe_number
    INTEGER                :: recbufsize, status(MPI_STATUS_SIZE)
    INTEGER,   ALLOCATABLE :: vrequest(:), vstatus(:,:)
    REAL(iwp), ALLOCATABLE :: vtempput(:,:), pl_pp(:)

    ALLOCATE(vrequest(numpesput))
    ALLOCATE(vstatus(MPI_STATUS_SIZE,numpesgetput))
    ALLOCATE(vtempput(neq_pp1,numpesput))
    ALLOCATE(pl_pp(0:len_pl_pp))

    vrequest   = 0 ; vstatus = 0 ; vtempput = 0 ; pl_pp(0:) = 0.0_iwp
    recbufsize = 0 ; ier     = 0 ; nstart   = 0 ; pe_number = 0

!------------------------------------------------------------------------------
! 1. Barrier to synchronise before gathering data (not always required, but 
!    safer and removes necessity for some barriers in main program).
!------------------------------------------------------------------------------

    ibar = 20
    CALL MY_BARRIER(numpe,ibar,details,'First barrier in gather')

!------------------------------------------------------------------------------
! 2. Calculate local p vector, pl_pp
!------------------------------------------------------------------------------

    pl_pp = 0.0

    DO ii = 1, numpesput
      i = pesput(ii)
      DO j = 1, lenput(i)
        vtempput(j,ii) = p_pp(toput(j,ii))
      END DO
      CALL MPI_ISEND(vtempput(1,ii),lenput(i),MPI_REAL8,i-1,numpe,            &
                     MPI_COMM_WORLD,vrequest(ii),ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B1) send',ier)
      END IF
    END DO

!------------------------------------------------------------------------------
! 3. Each PE must obtain required values of p
!    For data already available, just make a local copy:
!------------------------------------------------------------------------------

    nstart = lengetsum(numpe-1)
    DO j = 1,lenget(numpe)
      pl_pp(nstart + j) = p_pp(toget(j,getpes(numpe)))
    END DO

!------------------------------------------------------------------------------
! 4. For non local data must receive it:
!------------------------------------------------------------------------------

    DO i = 1,numpesget
      CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus(1,i),ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B1) probe',ier)
      END IF
      CALL MPI_GET_COUNT(vstatus(1,i),MPI_REAL8,recbufsize,ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B1) get_count',ier)
      END IF
      pe_number = vstatus(MPI_tag,i)
      nstart = lengetsum(pe_number-1)
      CALL MPI_RECV(pl_pp(nstart + 1),lenget(pe_number),MPI_REAL8,       &
       vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B1) receive',ier)
      END IF
    END DO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    IF (ier .NE. MPI_SUCCESS) THEN
      CALL MPERROR('Error in (B1) barrier after receives',ier)
    END IF

!------------------------------------------------------------------------------
! 5. Generate local pmul, pmul_pp
!------------------------------------------------------------------------------

    DO j = 1,nels_pp
      pmul_pp(:,j) = pl_pp(ggl_pp(:,j))
    END DO

!------------------------------------------------------------------------------
! 6. Make sure all ISENDs have completed and free up internal book-keeping 
!    of requests.
!------------------------------------------------------------------------------

    CALL MPI_WAITALL(numpesput,vrequest,vstatus,ier)

    IF (ier /= MPI_SUCCESS) THEN
      CALL MPERROR ('Error in MPI_WAITALL', ier)
    END IF

!------------------------------------------------------------------------------
! 7. Barrier to synchronise after gathering data (not always required, but 
!    safer and removes necessity for some barriers in main program).
!------------------------------------------------------------------------------

    ibar = 21
    CALL MY_BARRIER(numpe,ibar,details,'Last barrier in gather')

    DEALLOCATE(vrequest,vstatus,vtempput,pl_pp)

  END SUBROUTINE GATHER

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE SCATTER(u_pp,utemp_pp)
  
  !/****f* gather_scatter/scatter
  !*  NAME
  !*    SUBROUTINE: scatter
  !*  SYNOPSIS
  !*    Usage:      CALL scatter(u_pp,utemp_pp)
  !*  FUNCTION
  !*    Performs the scatter operation: u(g) = u(g) + utemp.
  !*  INPUTS
  !*    The following argument has the INTENT(IN) attribute:
  !*
  !*    utemp_pp(:)           : Real
  !*                          : Distributed array(NTOT,NELS_PP) that contains
  !*                          : data in an element-by-element form that is to 
  !*                          : be scattered to a vector of dimension NEQ_PP
  !*
  !*    The following argument has the INTENT(INOUT) attribute:
  !*
  !*    u_pp(:,:)             : Real
  !*                          : Distributed vector of dimension NEQ_PP.
  !*                  
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !*/ 

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)    :: utemp_pp(:,:)
    REAL(iwp), INTENT(INOUT) :: u_pp(:)
    INTEGER                  :: i, j, k, ii, ier, nstart, pe_number, ibar
    INTEGER                  :: recbufsize
    INTEGER                  :: status(MPI_STATUS_SIZE)
    INTEGER                  :: temp_status(MPI_STATUS_SIZE)
    LOGICAL                  :: lflag
    INTEGER,   ALLOCATABLE   :: vrequest(:), vstatus(:,:)
    REAL(iwp), ALLOCATABLE   :: ul_pp(:)

    status     = 0 ; temp_status = 0
    recbufsize = 0 ; ier         = 0 ; nstart = 0
    pe_number  = 0 ; ibar        = 0

    ALLOCATE(vrequest(numpesget))
    ALLOCATE(vstatus(MPI_STATUS_SIZE,numpesgetput))
    ALLOCATE(ul_pp(0:len_pl_pp))

    vrequest   = 0 ; vstatus     = 0 ; ul_pp  = 0.0_iwp

!------------------------------------------------------------------------------
! 1. Barrier to synchronise before scattering data (not always required, but 
!    safer and removes necessity for some barriers in main program).
!------------------------------------------------------------------------------

    ibar = 23
    CALL MY_BARRIER(numpe,ibar,details,'First barrier in scatter')

!------------------------------------------------------------------------------
! 2. Generate ul_pp vector
!------------------------------------------------------------------------------

    ul_pp = 0.0_iwp
    DO j = 1,nels_pp
      ul_pp(ggl_pp(:,j)) = ul_pp(ggl_pp(:,j)) + utemp_pp(:,j)
    END DO

!------------------------------------------------------------------------------
! 3. Generate u_pp vector - using ul_pp vectors from all appropriate PEs
!------------------------------------------------------------------------------

! Local data:

    nstart = lengetsum(numpe-1)
    DO j = 1, lenget(numpe)
      u_pp(toget(j,getpes(numpe))) = u_pp(toget(j,getpes(numpe))) +           &
                                     ul_pp(nstart+j)
    END DO

! Remote data:

    DO ii = 1, numpesget
      i = pesget(ii)
      nstart = lengetsum(i-1)
      CALL MPI_ISEND(ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,           &
                     MPI_COMM_WORLD,vrequest(ii),ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B5) isend',ier)
      END IF
    END DO

!------------------------------------------------------------------------------
! 4. Now receive data
!------------------------------------------------------------------------------

    DO i = 1, numpesput
      CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,               &
!                    vstatus(1,i),ier)
                     vstatus(:,i),ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B5) probe',ier)
      END IF
      CALL MPI_GET_COUNT(vstatus(1,i),MPI_REAL8,recbufsize,ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B5) get_count',ier)
      END IF
      pe_number = vstatus(MPI_TAG,i)

      IF (recbufsize .NE. lenput(pe_number)) THEN
        WRITE(*,*)'PE: ', numpe, ' Scatter error: recbufsize: ', recbufsize,  &
                  ' /= lenput(pe_number): ', lenput(pe_number)
        WRITE(*,*)' pe_number: ', pe_number, ' i: ', i
        WRITE(*,*)'source: ', vstatus(MPI_SOURCE,i),' tag: ',vstatus(MPI_TAG,i)
        CALL MPI_GET_COUNT (temp_status,MPI_REAL8,recbufsize,ier)
        IF (ier .NE. MPI_SUCCESS) THEN
          CALL MPERROR('Error in (B5) get_count',ier)
        END IF
        WRITE(*,*)'Temp probe. recbufsize: ', recbufsize, ' source: ',        &
                  temp_status(MPI_SOURCE),' tag: ', temp_status(MPI_TAG)
        STOP
      ENDIF  
      CALL MPI_RECV (tempput(1),lenput(pe_number),MPI_REAL8,                  &
       vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B5) receive',ier)
      END IF

!------------------------------------------------------------------------------
! 5. Note - need to use tempput to receive data, as must be added to
!    u_pp in loop below. No problem using the same temporary buffer
!    for the receives as these are blocking receives.
!------------------------------------------------------------------------------

      DO j = 1,lenput(pe_number)
        u_pp(toput(j,putpes(pe_number))) = u_pp(toput(j,putpes(pe_number))) + &
                                           tempput(j)
      END DO
    END DO

!------------------------------------------------------------------------------
! 6. Make sure all ISENDs have completed and free up internal book-keeping 
!    of requests.
!------------------------------------------------------------------------------

    CALL MPI_WAITALL(numpesget,vrequest,vstatus,ier)
    IF (ier /= MPI_SUCCESS) THEN
      CALL MPERROR('Error in MPI_WAITALL', ier)
    END IF

    ibar = 24
    CALL MY_BARRIER(numpe,ibar,details,'Last barrier in scatter')

    DEALLOCATE(vrequest,vstatus,ul_pp)

  END SUBROUTINE SCATTER

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE SCATTER_NOADD(utemp_pp,u_pp)

  !/****f* gather_scatter/scatter_noadd
  !*  NAME
  !*    SUBROUTINE: scatter_noadd
  !*  SYNOPSIS
  !*    Usage:      CALL scatter_noadd(utemp_pp,u_pp)
  !*  FUNCTION
  !*    Performs the scatter operation: u(g) = utemp.
  !*  INPUTS
  !*    The following argument has the INTENT(IN) attribute:
  !*
  !*    utemp_pp(:)           : Real
  !*                          : Distributed array(NTOT,NELS_PP) that contains
  !*                          : data in an element-by-element form that is to 
  !*                          : be scattered to a vector of dimension NEQ_PP
  !*
  !*    The following argument has the INTENT(INOUT) attribute:
  !*
  !*    u_pp(:,:)             : Real
  !*                          : Distributed vector of dimension NEQ_PP.
  !*                  
  !*  AUTHOR
  !*    M.A. Pettipher
  !*    V. Szeremi
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !*/ 

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)  :: utemp_pp(:,:)
    REAL(iwp), INTENT(OUT) :: u_pp(:)
    LOGICAL                :: lflag
    INTEGER                :: i, j, k, ii, ier, nstart, pe_number, l, ibar
    INTEGER                :: recbufsize
    INTEGER                :: status(MPI_STATUS_SIZE)
    INTEGER                :: temp_status(MPI_STATUS_SIZE)
    INTEGER,   ALLOCATABLE :: vrequest(:), vstatus(:,:)
    REAL(iwp), ALLOCATABLE :: ul_pp(:)

    ALLOCATE(vrequest(numpesget))
    ALLOCATE(vstatus(MPI_STATUS_SIZE,numpesgetput))
    ALLOCATE(ul_pp(0:len_pl_pp))

    vrequest   = 0 ; vstatus = 0 ; ul_pp  = 0.0_iwp
    recbufsize = 0 ; ier     = 0 ; nstart = 0
    pe_number  = 0 ; ibar    = 0

!------------------------------------------------------------------------------
! 1. Generate ul_pp vector
!------------------------------------------------------------------------------

    ul_pp = 0.0_iwp
    DO j = 1,nels_pp
      ul_pp(ggl_pp(:,j)) = utemp_pp(:,j)
    END DO

!------------------------------------------------------------------------------
! 2. Local data
!------------------------------------------------------------------------------

    nstart = lengetsum(numpe-1)
    DO j = 1, lenget(numpe)
      u_pp(toget(j,getpes(numpe))) = ul_pp(nstart + j)
    END DO

!------------------------------------------------------------------------------
! 3. Remote data
!------------------------------------------------------------------------------

    DO ii = 1,numpesget
      i = pesget(ii)
      nstart = lengetsum(i-1)
      CALL MPI_ISEND(ul_pp(nstart+1),lenget(i),MPI_REAL8,i-1,numpe,    &
                     MPI_COMM_WORLD,vrequest(ii),ier)
      IF (ier.NE.MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B5) isend',ier)
      END IF
    END DO

!------------------------------------------------------------------------------
! 4. Receive data
!------------------------------------------------------------------------------

    DO i = 1,numpesput
      CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,     &
                     vstatus(1,i),ier)
      IF (ier.NE.MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B5) probe',ier)
      END IF
      CALL MPI_GET_COUNT(vstatus(1,i),MPI_REAL8,recbufsize,ier)
      IF (ier.NE.MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B5) get_count',ier)
      END IF
      pe_number = vstatus(MPI_TAG,i)

      !#ifdef debug
      IF (recbufsize.NE.lenput(pe_number)) THEN
        WRITE(*,*)'PE: ', numpe, ' Scatter error: recbufsize: ', &
                  recbufsize, ' /= lenput(pe_number): ',lenput(pe_number)
        WRITE(*,*)' pe_number: ', pe_number, ' i: ',i
        WRITE(*,*)'source: ', vstatus(MPI_SOURCE,i),' tag: ',vstatus(MPI_TAG,i)
        CALL MPI_GET_COUNT (temp_status,MPI_REAL8,recbufsize,ier)
        IF (ier.NE.MPI_SUCCESS) THEN
          CALL MPERROR('Error in (B5) get_count',ier)
        END IF
        WRITE(*,*)'Temp probe. recbufsize: ', recbufsize, ' source: ', &
                  temp_status(MPI_SOURCE), ' tag: ', temp_status(MPI_TAG)
        STOP
      ENDIF
      !#endif

      CALL MPI_RECV(tempput(1),lenput(pe_number),MPI_REAL8,     &
      vstatus(MPI_SOURCE,i),vstatus(MPI_TAG,i),MPI_COMM_WORLD,vstatus(1,i),ier)
      IF (ier.NE.MPI_SUCCESS) THEN
        CALL MPERROR('Error in (B5) receive',ier)
      END IF
      DO j = 1,lenput(pe_number)
        u_pp(toput(j,putpes(pe_number))) = tempput(j)
      END DO
    END DO

!------------------------------------------------------------------------------
! 5. Make sure all ISENDs have completed and free up internal book-keeping
!    of requests
!------------------------------------------------------------------------------

    CALL MPI_WAITALL(numpesget,vrequest,vstatus,ier)
    IF (ier /= MPI_SUCCESS) THEN
      CALL MPERROR('Error in MPI_WAITALL',ier)
    END IF

    DEALLOCATE(vrequest,vstatus,ul_pp)

  END SUBROUTINE SCATTER_NOADD

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE MAKE_GGL(npes_pp,npes,gg_pp)

  !/****f* gather_scatter/make_ggl
  !*  NAME
  !*    SUBROUTINE: make_ggl
  !*  SYNOPSIS
  !*    Usage:      CALL make_ggl(npes_pp,npes,gg_pp)
  !*  FUNCTION
  !*    Generates ggl_pp and associated arrays.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*
  !*                  
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*
  !*     Version 1a, 25-11-96, M.A. Pettipher
  !*     Version 1b, 11-12-96, Replaces arguments with global variables
  !*     Version 1c, 27-01-97, Allows more than one PE to have neq_pp2 and 
  !*                           nels_pp2
  !*     Version 1d, 26-03-97, Uses vrequest and vstatus
  !*     Version 1e, 16-06-97, Reduces size of toget and toput by declaring
  !*                           second dimension as npes/4, and assigning from
  !*                           1 to numpesget (or numpesput) instead of to npes.
  !*                           Requires trap to ensure npes/4 is big enough.
  !*                           Subsequently reallocates with numpesget and
  !*                           numpesput.
  !*                           Creates getpes and putpes to provide reverse
  !*                           indexing of pesget and pesput.
  !*/
  
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: npes_pp, npes, gg_pp(ntot,nels_pp)
    INTEGER              :: i, j, k, ii, iii, ier, pe_number, bufid, count
    INTEGER              :: position, recbufsize
    INTEGER              :: vstatus(MPI_STATUS_SIZE,npes) 
                          ! numpesget not yet known, use npes
    INTEGER              :: vrequest(npes)      
                          ! numpesget not yet known, so use npes
    INTEGER              :: rem_acc, loc_acc, sum_rem_acc, sum_numpesget
    REAL(iwp)            :: sumtemp1, sumtemp2, rem_loc, sum_rem_loc
    LOGICAL              :: lflag, local, flag, newpe
    INTEGER, ALLOCATABLE :: preq_pp(:,:)

    ALLOCATE(preq_pp(neq_pp1,npes))
    preq_pp = 0

    recbufsize    = 0       ; ier           = 0      ; pe_number = 0
    bufid         = 0       ; count         = 0      
    position      = 0       ; rem_acc       = 0      ; loc_acc   = 0 
    sum_rem_acc   = 0       ; sum_numpesget = 0
    sum_numpesget = 0       
 
    sumtemp1      = 0.0_iwp ; sumtemp2      = 0.0_iwp
    rem_loc       = 0.0_iwp ; sum_rem_loc   = 0.0_iwp

!------------------------------------------------------------------------------
! 1. Call routine to allocate arrays for gather_scatter related operations:
!------------------------------------------------------------------------------

    CALL ALLOCATE_GATHER_SCATTER(npes_pp,npes)

    preq_pp   = 0
    threshold = num_neq_pp1*neq_pp1
    DO j = 1,nels_pp
      DO i = 1,ntot

!------------------------------------------------------------------------------
! 2. Convert index location, g, into location on specific PE. preq_pp(i,j) 
!    refers to location j on ith PE. Do not bother if index = 0 
!    (also prevents problem with assigning preq_pp(0, ))
!------------------------------------------------------------------------------

        IF (gg_pp(i,j).NE.0) THEN
          IF (gg_pp(i,j) <= threshold .OR. threshold == 0) THEN
            preq_pp(MOD((gg_pp(i,j)-1),neq_pp1)+1,(gg_pp(i,j)-1)/neq_pp1+1 ) = 1
          ELSE
            preq_pp(MOD((gg_pp(i,j)-threshold-1),(neq_pp1-1)) + 1,         &
                        (gg_pp(i,j)-threshold-1)/(neq_pp1-1)+1+num_neq_pp1 ) = 1
          ENDIF
        ENDIF
      END DO
    END DO

!------------------------------------------------------------------------------
! 3. Find number of elements required from each processor (lenget(i))
!    Also pesget and getpes - relating actual PE number, i,  to index, ii,
!    of toget(:,ii).
!    And numpesget - number of remote PEs from which data is required.
!------------------------------------------------------------------------------

    toget_temp = 0
    lenget     = 0
    lengetsum  = 0
    numpesget  = 0
    getpes     = 0
    ii         = 0
    local      = .FALSE.

    DO_find_data_outer: DO i = 1,npes
      newpe = .TRUE.
      k = 1
      DO_find_data_inner1: DO j = 1,neq_pp1
        IF (preq_pp(j,i).EQ.1) THEN
          IF (i == numpe) THEN             ! Save local until after all remote
            local = .TRUE.
            CYCLE DO_find_data_inner1
          END IF
          IF (newpe) THEN
            ii = ii + 1
          END IF
          IF (ii > npes_pp) THEN
            print *, 'PE: ',numpe,                                            &
                     '1 Number of PEs to get exceeds dimension of toget:',    &
                      ii, '>', npes_pp
            WRITE(details,'(A)') 'Need to increase npes_pp in main program'
            STOP
          END IF
          newpe            = .FALSE.
          toget_temp(k,ii) = j
          sumget(j,ii)     = k
          k                = k + 1
          lenget(i)        = lenget(i) + 1
        END IF
      END DO DO_find_data_inner1
      lengetsum(i) = lengetsum(i-1) + lenget(i)
      IF (.NOT. newpe) THEN
        pesget(ii) = i
        getpes(i) = ii
        WRITE(details,'("No. of elements of PE ",I4," required by PE ",I4,     &
                       & ": ",I8)') i, numpe, lenget(i)
      END IF
    END DO DO_find_data_outer

    numpesget = ii
    IF (local) THEN     ! For local data
      newpe = .TRUE.    ! Could avoid using newpe, but be consistent with above.
      i = numpe
      k = 1
      DO_find_data_inner2: DO j = 1, neq_pp1
        IF (preq_pp(j,i).EQ.1) THEN
          IF (newpe) THEN
            ii = ii + 1
          END IF
          IF (ii > npes_pp) THEN
            PRINT *, 'PE: ', numpe,                                            &
                     '2 Number of PEs to get exceeds dimension of toget:',     &
                      ii, '>', npes_pp
            WRITE(details,'(A)') 'Need to increase npes_pp in main program'
            STOP
          END IF
          newpe = .FALSE.
          toget_temp(k,ii) = j
          sumget(j,ii) = k
          k = k + 1
          lenget(i) = lenget(i) + 1
        END IF
      END DO DO_find_data_inner2
      lengetsum(i) = lengetsum(i-1) + lenget(i)
      DO iii = numpe + 1, npes
        lengetsum(iii) = lengetsum(iii) + lenget(numpe)
      END DO
      IF (.NOT. newpe) THEN
        pesget(ii) = i
        getpes(i) = ii
        WRITE(details,'("No. of elements of PE ",I4," required by PE ",I4,    &
                       & ": ",I8)') i, numpe, lenget(i)
      END IF
    END IF

    len_pl_pp = lengetsum(npes)
    WRITE(details,*)'Total number of unique elements required'
    WRITE(details,*)'i.e. length of pl_pp required: ', len_pl_pp
    WRITE(details,'("PE: ",I4," Number of remote PEs required: ",I8)') &
                   numpe, numpesget
    rem_acc = lengetsum(npes) - lenget(numpe)
    loc_acc = lenget(numpe)
    IF (loc_acc > 0) THEN
      rem_loc = rem_acc/REAL(loc_acc,iwp)
    ELSE
      rem_loc = 0
    ENDIF
    WRITE(*,'("PE: ",I4," Accesses - remote, local, remote/local: ",2I6,F8.2, &
              & " From ",I6, " PEs")') numpe, rem_acc, loc_acc, rem_loc,numpesget
    CALL MPI_REDUCE(rem_loc,sum_rem_loc,1,MPI_REAL8,MPI_SUM,      &
            npes-1,MPI_COMM_WORLD,ier)
    IF(ier .NE. MPI_SUCCESS) THEN
      CALL MPERROR('Error in REDUCE - rem_loc',ier)
    END IF
    CALL MPI_REDUCE(rem_acc,sum_rem_acc,1,MPI_INTEGER,MPI_SUM,    &
            npes-1,MPI_COMM_WORLD,ier)
    IF(ier .NE. MPI_SUCCESS) THEN
      CALL MPERROR('Error in REDUCE - rem_acc',ier)
    END IF
    CALL MPI_REDUCE(numpesget,sum_numpesget,1,MPI_INTEGER,MPI_SUM,&
      npes-1, MPI_COMM_WORLD,ier)
    IF(ier .NE. MPI_SUCCESS) THEN
      CALL MPERROR('Error in REDUCE - numpesget',ier)
    END IF
    IF (numpe == npes) THEN
      WRITE(*,'("Average accesses ratio - remote/local: ", F8.2)')            &
                sum_rem_loc/REAL(npes,iwp)
      WRITE(*,'("Total remote accesses                : ", I8)') sum_rem_acc
      IF (sum_numpesget > 0) THEN
        WRITE(*,'("Average remote accesses per PE       : ", F8.2)')          &
        sum_rem_acc/REAL(sum_numpesget,iwp)
      END IF
    ENDIF

    DEALLOCATE(preq_pp)

!------------------------------------------------------------------------------
! 4. Calculate ggl_pp
!    First find locations within preq_pp by summing (0 if element not required, 
!    1 if element required).
!------------------------------------------------------------------------------

    DO j = 1, nels_pp
      DO i = 1,ntot
        IF (gg_pp(i,j) .EQ. 0) THEN
          ggl_pp(i,j) = 0
        ELSE
          IF (gg_pp(i,j) <= threshold .OR. threshold == 0) THEN
            pe_number = (gg_pp(i,j) - 1)/neq_pp1 + 1
              k = gg_pp(i,j) - (pe_number - 1)*neq_pp1
          ELSE
            pe_number = num_neq_pp1 + (gg_pp(i,j) - threshold - 1)/           &
                       (neq_pp1 - 1) + 1
            k         = gg_pp(i,j) - threshold -                              &
                       (pe_number - num_neq_pp1 - 1)*(neq_pp1 - 1)
          ENDIF
          sumtemp1 = lengetsum(pe_number - 1)
          sumtemp2 = sumget(k,getpes(pe_number))
          ggl_pp(i,j) = sumtemp1 + sumtemp2
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
! 5. Send PE requirement request to each PE
!------------------------------------------------------------------------------

    DO ii = 1, numpesget
      i = pesget(ii)
      CALL MPI_ISEND (toget_temp(1,ii),lenget(i),MPI_INTEGER,i-1,numpe,  &
            MPI_COMM_WORLD,vrequest(ii),ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (A4) isend',ier)
      END IF
      CALL MPI_TEST(vrequest(ii),lflag,vstatus(1,ii),ier)
      WRITE(details,'("PE: ",I4," request sent to PE: ",I5)') numpe, i
    END DO

!------------------------------------------------------------------------------
! 6. Receive PE request. Now receive corresponding data from other PEs
!------------------------------------------------------------------------------

    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    IF (ier .NE. MPI_SUCCESS) THEN
      CALL MPERROR('Error in (A5) barrier before receives',ier)
    END IF

    toput_temp = 0

!------------------------------------------------------------------------------
! 7. Test for existence of message
!------------------------------------------------------------------------------

    count     = 0
    numpesput = 0
    lenput    = 0
    putpes    = 0
    ii        = 0
 
    DO                       ! Do until count limit reached
      CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,  &
          lflag,vstatus(1,numpesput+1),ier)
      IF (ier .NE. MPI_SUCCESS) THEN
        CALL MPERROR('Error in (A5) iprobe',ier)
      END IF
      IF (lflag) THEN
        count = 0
        ii = ii + 1
        numpesput = numpesput + 1
        IF (numpesput > npes_pp) THEN
          print *, 'PE: ', numpe, &
                'Number of PEs to put to exceeds dimension of toput:',        &
                 numpesput, '>', npes_pp
          STOP
        END IF
            

        CALL MPI_GET_COUNT (vstatus(1,numpesput),MPI_INTEGER,recbufsize,ier)
        IF (ier .NE. MPI_SUCCESS) THEN
          CALL MPERROR('Error in (A5) get_count',ier)
        END IF
        pe_number = vstatus(MPI_tag,numpesput)
        lenput(pe_number) = recbufsize 
        pesput(ii) = pe_number
        putpes(pe_number) = ii
        CALL MPI_RECV (toput_temp(1,ii),lenput(pe_number),MPI_INTEGER,        &
        MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus(1,numpesput),ier)
        IF (ier .NE. MPI_SUCCESS) THEN
          CALL MPERROR('Error in (A5) receive',ier)
        END IF
      ELSE
        count = count + 1
        IF (count .EQ. 100) THEN  
        ! See issue 20 http://code.google.com/p/parafem/issues/list
          EXIT
        END IF
      END IF

    END DO

!------------------------------------------------------------------------------
! 8. Make sure all ISENDs have completed and free up internal book-keeping 
!    of requests
!------------------------------------------------------------------------------

    CALL MPI_WAITALL(numpesget,vrequest,vstatus,ier)
    IF (ier /= MPI_SUCCESS) THEN
      CALL MPERROR ('Error in MPI_WAITALL', ier)
    END IF

!------------------------------------------------------------------------------
! 9. Set numpesgetput to max of numpesget and numpesput
!------------------------------------------------------------------------------

    numpesgetput = MAX(numpesget, numpesput)

!------------------------------------------------------------------------------
! 10. Print more information about required communication.
!------------------------------------------------------------------------------

    DO i = 1,numpesput
      WRITE(details,'("Number of elements of PE ",I4," required by PE ",      &
                     & I4,": ",I8)') numpe, pesput(i), lenput(pesput(i))
    END DO
    WRITE(details,'("PE: ", I4," Number of PEs to send data to: ", I4)')      &
        numpe, numpesput
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    IF (ier .NE. MPI_SUCCESS) THEN
      CALL MPERROR('Error in (A5) barrier',ier)
    END IF

!------------------------------------------------------------------------------
! 11. Now numpesget, numpesput, toget and toput are known, reallocate toget 
!     and toput. Note that numpesget is number of remote PEs required, but if 
!     as would normally be expected, some of the data is local (on numpe), 
!     toget requires an extra element for the local data. Numpesput is number 
!     to be sent to (excluding numpe).
!------------------------------------------------------------------------------

    IF (local) THEN
      ALLOCATE ( toget(neq_pp1, numpesget + 1), toput(neq_pp1,numpesput) )      
      toget = 0 
      toput = 0
      toget = toget_temp(:, 1:numpesget + 1)
    ELSE
      ALLOCATE ( toget(neq_pp1, numpesget), toput(neq_pp1,numpesput) )      
      toget = 0
      toput = 0
      toget = toget_temp(:, 1:numpesget)
    END IF
    toput = toput_temp(:, 1:numpesput)
    DEALLOCATE ( toget_temp, toput_temp )

!------------------------------------------------------------------------------
! 12. The following barrier is probably necessary, but not checked thoroughly
!------------------------------------------------------------------------------

    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    IF (ier .NE. MPI_SUCCESS) THEN
      CALL MPERROR('Error in last barrier in make_ggl',ier)
    END IF

  END SUBROUTINE MAKE_GGL

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE SCATTER_NODES(npes,nn,nels_pp,g_num_pp,nod,numvar,nodes_pp,      &
                    node_start,node_end,element_contribution,nodal_value,flag)

    !/****f* gather_scatter/scatter_nodes
    !*  NAME
    !*    SUBROUTINE: scatter_nodes
    !*  SYNOPSIS
    !*    Usage:      CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,        &
    !*                            numvar,nodes_pp,node_start,node_end,        &
    !*                            element_contribution,nodal_value,flag)
    !*  FUNCTION
    !*    Assembly element contributions of a nodal value
    !*    
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    flag                  : Integer
    !*                          : Indicator (1 average sum, 0 total sum)
    !*
    !*    nels_pp               : Integer
    !*                          : Number of elements per processor
    !*
    !*    nn                    : Integer
    !*                          : Total number of nodes
    !*
    !*    nod                   : Integer
    !*                          : Number of nodes per element
    !*
    !*    node_end              : Integer
    !*                          : Last node stored in the process
    !*
    !*    node_start            : Integer
    !*                          : First node stored in the process
    !*
    !*    nodes_pp              : Integer
    !*                          : Number of nodes stored in the process
    !*
    !*    npes                  : Integer
    !*                          : Number of processes
    !*
    !*    numvar                : Integer
    !*                          : Number of components of the variable
    !*                            (1-scalar, 3-vector, 6-tensor)
    !*
    !*    g_num_pp(nod,nels_pp) : Integer
    !*                          : Elements connectivity
    !*
    !*    element_contribution(ntot,nels_pp) : Real
    !*                                       : Elements contribution to the 
    !*                                         nodal variable
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    nodal_value(nodes_pp*nodof)        : Real
    !*                                       : Value after assembling elements
    !*                                         contribution
    !*  AUTHOR
    !*    F. Calvo
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE

    INTEGER,   INTENT(IN)  :: flag, nels_pp, nn, nod, node_end, node_start, &
                              nodes_pp, npes, numvar, g_num_pp(:,:)
    REAL(iwp), INTENT(IN)  :: element_contribution(:,:)
    REAL(iwp), INTENT(OUT) :: nodal_value(:)
    INTEGER                :: i, j, k, iel, startNode, endNode, sizee, idx1
    INTEGER                :: idx2, bufsize, ier
    INTEGER, ALLOCATABLE   :: itemp(:), itempdest(:)
    REAL(iwp), ALLOCATABLE :: temp(:)

!------------------------------------------------------------------------------
! 1. Allocate internal arrays
!------------------------------------------------------------------------------
    
    sizee = (nodes_pp+1)*numvar
    ALLOCATE( itemp(sizee), itempdest(sizee), temp(sizee) )
    itemp = 0 ; itempdest = 0; temp = 0.0_iwp

!------------------------------------------------------------------------------
! 2. Every process (first loop) sends the range of nodes assigned to that
!    process. The other processes check whether the elements contain nodes
!    in that range for the contribution. The contributions are summed up
!    and the process collects the results
!------------------------------------------------------------------------------
    
    bufsize = 1

    DO i = 1,npes
      startNode = node_start
      endNode   = node_end
      CALL MPI_BCAST(startNode,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST(endNode,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)

      sizee = nodes_pp*numvar
      CALL MPI_BCAST(sizee,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)

      temp = 0.0_iwp
      itemp = 0

      DO iel = 1,nels_pp
        DO j = 1,nod
          IF (g_num_pp(j,iel)>=startNode .and. g_num_pp(j,iel)<=endNode) THEN
            idx1 = (g_num_pp(j,iel)-startNode)*numvar
            idx2 = (j-1)*numvar
            DO k = 1,numvar
              temp(idx1+k)  = temp(idx1+k) + element_contribution(idx2+k,iel)
              itemp(idx1+k) = itemp(idx1+k) + 1
            END DO
          END IF
        END DO
      END DO

      CALL MPI_REDUCE(temp, nodal_value, sizee, MPI_REAL8, MPI_SUM, i-1, &
                      MPI_COMM_WORLD, ier)
      CALL MPI_REDUCE(itemp, itempdest, sizee, MPI_INTEGER, MPI_SUM, i-1, &
                      MPI_COMM_WORLD, ier)
    END DO

!------------------------------------------------------------------------------
! 3. Compute average if required
!------------------------------------------------------------------------------

    IF (flag==1) THEN
      DO i = 1,nodes_pp*numvar
        IF (itempdest(i) == 0) THEN
          write(*,*)"Error in average: itemp = 0"
        END IF
        nodal_value(i) = nodal_value(i)/itempdest(i)
      END DO
    END IF

!------------------------------------------------------------------------------
! 4. Deallocate internal arrays
!------------------------------------------------------------------------------

    DEALLOCATE(itemp,itempdest,temp)
      
  END SUBROUTINE SCATTER_NODES
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------ 
!------------------------------------------------------------------------------

  SUBROUTINE SCATTER_NODES_NS(npes,nn,nels_pp,g_num_pp,nod,numvar,nodes_pp,   &
                    node_start,node_end,element_contribution,nodal_value,flag)

    !/****f* gather_scatter/scatter_nodes_ns
    !*  NAME
    !*    SUBROUTINE: scatter_nodes_ns
    !*  SYNOPSIS
    !*    Usage:      CALL scatter_nodes_ns(npes,nn,nels_pp,g_num_pp,nod,     &
    !*                            numvar,nodes_pp,node_start,node_end,        &
    !*                            element_contribution,nodal_value,flag)
    !*  FUNCTION
    !*    Assemble element contributions of a nodal value for program p126
    !*    in Smith and Griffiths "Programming the Finite Element Method" (Ed 4)
    !*    Deals with coupled 20 node velocity elements and 8 node pressure 
    !*    elements
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    flag                  : Integer
    !*                          : Indicator (1 average sum, 0 total sum)
    !*
    !*    nels_pp               : Integer
    !*                          : Number of elements per processor
    !*
    !*    nn                    : Integer
    !*                          : Total number of nodes
    !*
    !*    nod                   : Integer
    !*                          : Number of nodes per element
    !*
    !*    node_end              : Integer
    !*                          : Last node stored in the process
    !*
    !*    node_start            : Integer
    !*                          : First node stored in the process
    !*
    !*    nodes_pp              : Integer
    !*                          : Number of nodes stored in the process
    !*
    !*    npes                  : Integer
    !*                          : Number of processes
    !*
    !*    numvar                : Integer
    !*                          : Number of components of the variable
    !*                            (1-scalar, 3-vector, 6-tensor)
    !*
    !*    g_num_pp(nod,nels_pp) : Integer
    !*                          : Elements connectivity
    !*
    !*    element_contribution(ntot,nels_pp) : Real
    !*                                       : Elements contribution to the 
    !*                                         nodal variable
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    nodal_value(nodes_pp*nodof)        : Real
    !*                                       : Value after assembling elements
    !*                                         contribution
    !*  AUTHOR
    !*    L. Margetts
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE

    INTEGER,   INTENT(IN)  :: flag, nels_pp, nn, nod, node_end, node_start, &
                              nodes_pp, npes, numvar, g_num_pp(:,:)
    REAL(iwp), INTENT(IN)  :: element_contribution(:,:)
    REAL(iwp), INTENT(OUT) :: nodal_value(:)
    INTEGER                :: i, j, k, iel, startNode, endNode, sizee, idx1
    INTEGER                :: idx2, bufsize, ier
    REAL(iwp)              :: zero = 0.0_iwp
    INTEGER, ALLOCATABLE   :: itemp(:), itempdest(:)
    REAL(iwp), ALLOCATABLE :: temp(:)

!------------------------------------------------------------------------------
! 1. Allocate internal arrays
!------------------------------------------------------------------------------
    
    sizee = (nodes_pp+1)*numvar
    ALLOCATE( itemp(sizee), itempdest(sizee), temp(sizee) )
    itemp = 0 ; itempdest = 0 ; temp = 0.0_iwp

!------------------------------------------------------------------------------
! 2. Every process (first loop) sends the range of nodes assigned to that
!    process. The other processes check whether the elements contain nodes
!    in that range for the contribution. The contributions are summed up
!    and the process collects the results
!------------------------------------------------------------------------------
    
    bufsize = 1

    DO i = 1,npes
      startNode = node_start
      endNode   = node_end
      CALL MPI_BCAST(startNode,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST(endNode,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)

      sizee = nodes_pp*numvar
      CALL MPI_BCAST(sizee,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)

      temp = 0.0
      itemp = 0

      DO iel = 1,nels_pp
        DO j = 1,nod
          IF (g_num_pp(j,iel)>=startNode .and. g_num_pp(j,iel)<=endNode) THEN
            idx1 = (g_num_pp(j,iel)-startNode)*numvar
!           idx2 = (j-1)*numvar

!           DO k = 1,numvar
!             temp(idx1+k)  = temp(idx1+k) + element_contribution(idx2+k,iel)
!             itemp(idx1+k) = itemp(idx1+k) + 1
!           END DO
            
            idx2 = j
            temp(idx1+1)  = temp(idx1+1) + element_contribution(idx2,iel)
            itemp(idx1+1) = itemp(idx1+1) + 1
 
            IF(j==1) THEN
              idx2 = 20+1
              temp(idx1+2)  = temp(idx1+2) + element_contribution(idx2,iel)
              itemp(idx1+2) = itemp(idx1+2) + 1
            ELSE IF(j==3) THEN
              idx2 = 20+2
              temp(idx1+2)  = temp(idx1+2) + element_contribution(idx2,iel)
              itemp(idx1+2) = itemp(idx1+2) + 1
            ELSE IF(j==5) THEN
              idx2 = 20+3
              temp(idx1+2)  = temp(idx1+2) + element_contribution(idx2,iel)
              itemp(idx1+2) = itemp(idx1+2) + 1
            ELSE IF(j==7) THEN
              idx2 = 20+4
              temp(idx1+2)  = temp(idx1+2) + element_contribution(idx2,iel)
              itemp(idx1+2) = itemp(idx1+2) + 1
            ELSE IF(j==13) THEN
              idx2 = 20+5
              temp(idx1+2)  = temp(idx1+2) + element_contribution(idx2,iel)
              itemp(idx1+2) = itemp(idx1+2) + 1
            ELSE IF(j==15) THEN
              idx2 = 20+6
              temp(idx1+2)  = temp(idx1+2) + element_contribution(idx2,iel)
              itemp(idx1+2) = itemp(idx1+2) + 1
            ELSE IF(j==17) THEN
              idx2 = 20+7
              temp(idx1+2)  = temp(idx1+2) + element_contribution(idx2,iel)
              itemp(idx1+2) = itemp(idx1+2) + 1
            ELSE IF(j==19) THEN
              idx2 = 20+8
              temp(idx1+2)  = temp(idx1+2) + element_contribution(idx2,iel)
              itemp(idx1+2) = itemp(idx1+2) + 1
            ELSE
              temp(idx1+2)  = temp(idx1+2) + zero
              itemp(idx1+2) = itemp(idx1+2) + 1
            END IF

            idx2 = j+28
            temp(idx1+3)  = temp(idx1+3) + element_contribution(idx2,iel)
            itemp(idx1+3) = itemp(idx1+3) + 1

            idx2 = j+48
            temp(idx1+4)  = temp(idx1+4) + element_contribution(idx2,iel)
            itemp(idx1+4) = itemp(idx1+4) + 1
 

          END IF
        END DO
      END DO

      CALL MPI_REDUCE(temp, nodal_value, sizee, MPI_REAL8, MPI_SUM, i-1, &
                      MPI_COMM_WORLD, ier)
      CALL MPI_REDUCE(itemp, itempdest, sizee, MPI_INTEGER, MPI_SUM, i-1, &
                      MPI_COMM_WORLD, ier)
    END DO

!------------------------------------------------------------------------------
! 3. Compute average if required
!------------------------------------------------------------------------------

    IF (flag==1) THEN
      DO i = 1,nodes_pp*numvar
        IF (itempdest(i) == 0) THEN
          write(*,*)"Error in average: itemp = 0"
        END IF
        nodal_value(i) = nodal_value(i)/itempdest(i)
      END DO
    END IF

!------------------------------------------------------------------------------
! 4. Deallocate internal arrays
!------------------------------------------------------------------------------

    DEALLOCATE(itemp,itempdest,temp)
      
  END SUBROUTINE SCATTER_NODES_NS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
END MODULE GATHER_SCATTER
