PROGRAM ndttrget
IMPLICIT NONE

INTEGER              :: nstep,nres,i,j,k
REAL                 :: ttr
CHARACTER(LEN=50)    :: job_name,arg1,arg2,fname

!------------------------------------------------------------------------------
! 1. Read arguments
!------------------------------------------------------------------------------

CALL GETARG(1,job_name)
CALL GETARG(2,arg1)
CALL GETARG(3,arg2)

read(arg1,*) nstep !Convert string to integer
read(arg2,*) nres

PRINT *, "job_name = ", job_name
PRINT *, "nstep = ", nstep
PRINT *, "nres = ", nres

!------------------------------------------------------------------------------
! 2. Open files, read and write ttr
!------------------------------------------------------------------------------

fname=job_name(1:INDEX(job_name, " ")-1)//".ndttr"
PRINT *, "open: ", fname
OPEN(10,file=fname,status='replace',action='write')

DO i=1,nstep
  WRITE(arg1,'(i10)')i !Convert integer to string
  arg1=ADJUSTL(arg1)   !Remove trailing spaces
  IF (i<10) fname=job_name(1:INDEX(job_name, " ")-1)//".ensi.NDTTR-00000"//arg1
  IF (i<100 .AND. i>9) fname=job_name(1:INDEX(job_name, " ")-1)//".ensi.NDTTR-0000"//arg1
  IF (i<1000 .AND. i>99) fname=job_name(1:INDEX(job_name, " ")-1)//".ensi.NDTTR-000"//arg1
  IF (i<10000 .AND. i>999) fname=job_name(1:INDEX(job_name, " ")-1)//".ensi.NDTTR-00"//arg1
  IF (i<100000 .AND. i>9999) fname=job_name(1:INDEX(job_name, " ")-1)//".ensi.NDTTR-0"//arg1
  IF (i<1000000 .AND. i>99999) fname=job_name(1:INDEX(job_name, " ")-1)//".ensi.NDTTR-"//arg1
  PRINT *, "read: ", fname
  OPEN(11,file=fname,status='old',action='read')
  
  DO j=1,4+nres-1      !Skip header and ttr values before nres
     READ (11,*)
  END DO
  READ (11,*)ttr
  WRITE(10,'(E16.8)')ttr
  CLOSE(11)
END DO
  
CLOSE(10)

END PROGRAM ndttrget
