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
! 5. Write ".dis" file
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dis"
  OPEN(20,FILE=fname,STATUS='REPLACE',ACTION='WRITE') 
  
  WRITE(20,'(A)') "*DISPLACEMENT"
  WRITE(20,'(A)') " 1" 

  DO i = 1, nn
    WRITE(20,'(I8,3E20.12)') i, upvw(1,i), upvw(3,i), upvw(4,i)
  END DO

  WRITE(20,'(A)') "*DISPLACEMENT"
  WRITE(20,'(A)') " 2" 

  DO i = 1, nn
    WRITE(20,'(I8, E20.12,2A)') i, upvw(2,i),"  0.0000E+00", "  0.0000E+00"
  END DO

  
END PROGRAM upvw2dis
