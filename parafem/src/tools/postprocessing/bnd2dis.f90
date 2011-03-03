PROGRAM bnd2dis

  !/****h* /bnd2dis
  !*  NAME
  !*    PROGRAM: bnd2dis
  !*  SYNOPSIS
  !*    Usage:   ./bnd2dis <job_name>
  !*  FUNCTION
  !*    Reads a ParaFEM boundary condition file <job_name.bnd> and outputs
  !*    a dummy displacements file <job_name_dummy.dis> where restrained
  !*    freedoms are given a value of 1.0.
  !*
  !*    This enables the user to visualize restrained freedoms in the 
  !*    ParaFEM-Viewer.
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  This utility is only needed as an interim measure until the 
  !*  ParaFEM-Viewer is updated to view <job_name>.bnd files
  !*
  !*/

  USE precision
  USE global_variables
  USE timing

  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER,PARAMETER     :: ndim = 3,nodof = 3
  INTEGER               :: i,j,k,count,nn,nr,skip
  INTEGER               :: argc,iargc
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp,one = 1.0_iwp
  CHARACTER(LEN=50)     :: job_name
  CHARACTER(LEN=50)     :: fname

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  INTEGER,ALLOCATABLE   :: rest(:,:)
  REAL(iwp),ALLOCATABLE :: disp(:),timest(:)

!------------------------------------------------------------------------------
! 3. Read job_name from the command line.
!    Read control data and mesh data
!------------------------------------------------------------------------------

  ALLOCATE(timest(1))
  timest    = zero
  timest(1) = elap_time()

  argc = iargc()

  IF(argc /= 1) THEN
    WRITE(*,'(/A)') " Usage: bnd2dis <job_name> "
    WRITE(*,*)
    WRITE(*,'(A)')  "        bnd2dis expects argument <job_name> and 2 input"
    WRITE(*,'(A/)') "        files <job_name>.dat and <job_name>.bnd"
    STOP 
  END IF

  CALL GETARG(1, job_name)

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
  READ(10,*) 
  READ(10,*)
  READ(10,*) skip, nn, nr
  CLOSE(10)

  PRINT *, "nn = ",nn
  PRINT *, "nr = ",nr

  ALLOCATE(rest(nodof+1,nr))
  ALLOCATE(disp(ndim))

  rest  = 0
  disp  = zero

!------------------------------------------------------------------------------
! 4. Read ".bnd" boundary conditions file
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ") -1) // ".bnd"
  OPEN(11,FILE=fname,STATUS='OLD',ACTION='READ')

  DO i = 1,nr
    READ(11,*) rest(:,i)   
    PRINT *, i, rest(:,i) 
  END DO

  PRINT *, "Read ",nr," restrained nodes."

!------------------------------------------------------------------------------
! 5. Write <job_name_dummy.bnd> file
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ") -1) // "-dummy.dis"
  OPEN(20,FILE=fname,STATUS='REPLACE',ACTION='WRITE') 

  WRITE(20,'(A)') "*DISPLACEMENT"
  WRITE(20,'(A)') "      1"

  count = 0

  DO i = 1,nn

    disp  = zero

    DO j = 1,nr
      IF(rest(1,j) == i) THEN
        DO k = 1,nodof
          IF(rest(k+1,j) == 0) THEN
            disp(k) = one
          END IF
        END DO
        count = count + 1
        CYCLE ! found
      END IF
    END DO

    WRITE(20,'(I10,3E20.12)') i, disp(:)

  END DO

  CLOSE(20)

  PRINT *, "bnd2dis completed in:",elap_time()-timest(1),"s."
 
END PROGRAM bnd2dis
