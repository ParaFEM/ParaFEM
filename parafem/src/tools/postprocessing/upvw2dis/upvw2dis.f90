PROGRAM upvw2dis

  !/****h* /upvw2dis
  !*  NAME
  !*    PROGRAM: upvw2dis
  !*  SYNOPSIS
  !*    Usage:   ./upvw2dis <file>
  !*  FUNCTION
  !*    Splits the mixed velocity and pressure data in the .upvw file output by 
  !*    p126 into separate fields. This is required to view the results in the
  !*    ParaFEM-Viewer
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
  INTEGER               :: nodof = 4
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
  CHARACTER(LEN=50)     :: program_name

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: upvw(:,:)

!------------------------------------------------------------------------------
! 3. Read job_name from the command line.
!    Read control data and mesh data
!------------------------------------------------------------------------------

  argc = iargc()

  IF(argc /= 1) THEN
    WRITE(*,*)
    WRITE(*,'(A)') " Usage: upvw2dis <job_name> "
    WRITE(*,*)
    WRITE(*,'(A)') "        upvw2dis expects argument <job_name> and 2 input"
    WRITE(*,'(A)') "        files <job_name>.dat and <job_name>.upvw"
    WRITE(*,*)
    WRITE(*,'(A)') "********************* PROGRAM ABORTED *********************"
    WRITE(*,*)
    STOP 
  END IF

  CALL GETARG(1, job_name)

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
  READ(10,*) program_name
  READ(10,*) element,meshgen,nels,nn
  CLOSE(10)

!------------------------------------------------------------------------------
! 4. Read ".upvw" file
!------------------------------------------------------------------------------
 
  ALLOCATE(upvw(nodof,nn))
  upvw   = 0

  fname = job_name(1:INDEX(job_name, " ") -1) // ".upvw"
  OPEN(11,FILE=fname,STATUS='OLD',ACTION='READ')

  READ(11,*)
  READ(11,*)

  DO i = 1,nn
    READ(11,*) id, upvw(:,i)    
  END DO

!------------------------------------------------------------------------------
! 5. Write ".dis" and ".ttr" files
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dis"
  OPEN(20,FILE=fname,STATUS='REPLACE',ACTION='WRITE') 
  
  WRITE(20,'(A)') "*DISPLACEMENT"
  WRITE(20,'(A)') " 1" 

  DO i = 1, nn
    WRITE(20,'(I8,3E20.12)') i, upvw(1,i), upvw(3,i), upvw(4,i)
  END DO
  CLOSE(20)

  fname = job_name(1:INDEX(job_name, " ") -1) // ".ttr"
  OPEN(30,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  WRITE(30,'(A)') "*TEMPERATURE"
  WRITE(30,'(A)') " 1" 

  DO i = 1, nn
    WRITE(30,'(I8, E20.12)') i, upvw(2,i)
  END DO
  CLOSE(30)
  
END PROGRAM upvw2dis
