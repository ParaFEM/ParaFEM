PROGRAM ttrb2ttr
IMPLICIT NONE

INTEGER              :: i,j,k,l
INTEGER              :: nn
INTEGER              :: nstep
INTEGER              :: npes
INTEGER              :: bufsize,step,argc,iargc
INTEGER, ALLOCATABLE :: rl(:),ps(:)
REAL, ALLOCATABLE    :: ttr(:)
CHARACTER(LEN=50)    :: fname,job_name,label

argc = iargc()
CALL GETARG(1,job_name)
fname=job_name(1:INDEX(job_name, " ")-1)//".ttrb"
OPEN(10,file=fname,status='old',action='read',access='sequential',          &
     form='unformatted')
fname=job_name(1:INDEX(job_name, " ")-1)//".ttr"
OPEN(11,file=fname,status='replace',action='write')
fname=job_name(1:INDEX(job_name, " ")-1)//".npp"
OPEN(12,file=fname,status='old',action='read')

!------------------------------------------------------------------------------
! 0. Read record lengths
!------------------------------------------------------------------------------

READ(12,*) nn
READ(12,*) nstep
READ(12,*) npes

ALLOCATE(rl(npes))

DO i=1,npes
  READ(12,*) rl(i)
END DO

bufsize = maxval(rl)

!------------------------------------------------------------------------------
! 1. Read binary
!------------------------------------------------------------------------------

ALLOCATE(ttr(bufsize))
ttr=0.0

DO i=1,nstep+1   ! +1 for step 0
  READ(10) label
  WRITE(11,'(A)') label
  READ(10) step
  WRITE(11,*) step
  l=0
  DO j=1,npes
    ttr=0.0
    READ(10) ttr(1:rl(j))
    DO k=1,rl(j)
     l=l+1
     WRITE(11,'(I8,(1P,E12.4))') l, ttr(k)
    END DO
  END DO
END DO

END PROGRAM ttrb2ttr
