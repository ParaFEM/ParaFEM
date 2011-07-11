PROGRAM mat2numpe

  !/****h* /mat2numpe
  !*  NAME
  !*    PROGRAM: mat2numpe
  !*  SYNOPSIS
  !*    Usage:   ./mat2numpe <job_name>
  !*  FUNCTION
  !*    Reads a ParaFEM file <job_name.d>, swaps the material ID for the 
  !*    processor ID, and outputs a new ParaFEM file <job_name-npes.d>
  !*
  !*    This enables the user to visualize the processor partitions in the 
  !*    ParaFEM-Viewer.
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
  USE gather_scatter
  USE timing

  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER     :: ndim  = 3
  INTEGER               :: nod   = 4
  INTEGER               :: partitioner = 1
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

  INTEGER,ALLOCATABLE   :: g_num_d(:,:)
  INTEGER,ALLOCATABLE   :: nels_pp_store(:)
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
    WRITE(*,'(A)') " Usage: mat2numpe <job_name> "
    WRITE(*,*)
    WRITE(*,'(A)') "        mat2numpe expects argument <job_name> and 2 input"
    WRITE(*,'(A)') "        files <job_name>.dat and <job_name>.d"
    WRITE(*,*)
    STOP 
  END IF

  CALL GETARG(1, job_name)

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
  READ(10,*) element, nels, nn, npes
  CLOSE(10)

  ALLOCATE(g_num_d(nod+5,nels))
  ALLOCATE(g_coord(ndim,nn))
  ALLOCATE(nels_pp_store(npes))

  g_num_d       = 0
  g_coord       = zero
  nels_pp_store = 0

!------------------------------------------------------------------------------
! 4. Read ".d" geometry file
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ") -1) // ".d"
  OPEN(11,FILE=fname,STATUS='OLD',ACTION='READ')

  READ(11,*)
  READ(11,*)

  DO i = 1,nn
    READ(11,*) id, g_coord(:,i)    
  END DO

  READ(11,*)
  DO iel = 1,nels 
    READ(11,*) g_num_d(:,iel)
  END DO

!------------------------------------------------------------------------------
! 4. Compute nels_pp and iel_start for each processor and store in array
!------------------------------------------------------------------------------

  numpe = 0

  DO i = 1, npes 
    numpe = numpe + 1
    CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)
    nels_pp_store(i) = nels_pp
  END DO

!------------------------------------------------------------------------------
! 5. Write <job_name-npes.d> file
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ") -1) // "-npes.d"
  OPEN(20,FILE=fname,STATUS='REPLACE',ACTION='WRITE') 
  
  WRITE(20,'(A)') "*THREE_DIMENSIONAL"
  WRITE(20,'(A)') "*NODES"

  DO i = 1, nn
    WRITE(20,'(3E20.12)') g_coord(:,i)
  END DO

  WRITE(20,'(A)') "*ELEMENTS"

  iel = 0

  DO i = 1, npes
    nels_pp = nels_pp_store(i)
    DO j = 1, nels_pp 
      iel = iel + 1
      WRITE(20,'(9I8)') g_num_d(1:8,iel), i
    END DO
  END DO

 
END PROGRAM mat2numpe
