PROGRAM rfem

!/****h* tools/preprocessing/rfem
!*  NAME
!*    PROGRAM rfem
!*  FUNCTION
!*    Creates box of finite elements and randomly generated field of elastic
!*    properties
!*  AUTHOR
!*    Lee Margetts
!*  COPYRIGHT
!*    (c) University of Manchester 2012
!****
!*/

  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  CHARACTER(LEN=15)      :: field_name

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------   

  INTEGER, ALLOCATABLE   :: 
  REAL(iwp), ALLOCATABLE :: 

!------------------------------------------------------------------------------
! 3. Read job_name from the command line
!------------------------------------------------------------------------------

  argc = iargc()
  IF (argc /= 1) THEN
    PRINT*
    PRINT*, "Usage:  rfem <job_name>"
    PRINT*
    PRINT*, "        program expects <job_name>.?? and outputs"
    PRINT*, "        <job_name>.d" 
    PRINT*
    STOP
  END IF
  CALL GETARG(1, job_name)

!------------------------------------------------------------------------------
! 4. Read control data
!------------------------------------------------------------------------------

  IF (INDEX(job_name,".mg") /= 0) THEN
    job_name = job_name(1:INDEX(job_name,".mg")-1)
  END IF
  fname = job_name(1:INDEX(job_name," ")-1) // ".mg"
  OPEN (10, file=fname, status='old', action='read')

  READ(10,*) program_name

  PRINT *
  PRINT *, "Generating input deck for program ", program_name
  PRINT *

!------------------------------------------------------------------------------
! 5. Select random field generator
!------------------------------------------------------------------------------

  SELECT CASE(field_type)

  CASE('sim3de')

    
  CASE DEFAULT
  
    PRINT *
    PRINT *, "Wrong value given in variable field_name"
    PRINT *
      
  END SELECT
  

END PROGRAM rfem