PROGRAM ndttrget
!------------------------------------------------------------------------------
!      Program ndttrget    Tool to produce read temperature values from 
!                          .ensi.NDTTR files (output from program xx12.f90)
!                          for a specific node and writes to .ndttr
!                          Under development, currently only for xx12 using 
!                          4-node tetrahedra
!
!      USAGE: ndttrget <job_name> <nstep> <nres> <format>
!------------------------------------------------------------------------------
  
  USE, INTRINSIC :: ISO_C_BINDING ! to output C binary file
  
  IMPLICIT NONE
  
  INTEGER,PARAMETER    :: iwp=SELECTED_REAL_KIND(15)
  INTEGER              :: nstep,nres,i,j,k
  INTEGER(KIND=C_INT)  :: int_in
  REAL(KIND=C_FLOAT)   :: ttr_in
  REAL(iwp)            :: ttr
  CHARACTER(LEN=50)    :: job_name,arg,fname,stepnum,format
  CHARACTER(LEN=80)    :: cbuffer

!------------------------------------------------------------------------------
! 1. Read arguments
!------------------------------------------------------------------------------
  
  CALL GETARG(1,job_name)
  CALL GETARG(2,arg)
  read(arg,*) nstep !Convert string to integer
  CALL GETARG(3,arg)
  read(arg,*) nres !Convert string to integer
  CALL GETARG(4,format)
  
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
    
    WRITE(stepnum,'(I0.6)') i
    
    SELECT CASE (format)
      CASE('bin')
        fname=job_name(1:INDEX(job_name, " ")-1) // ".bin.ensi.NDTTR-"//stepnum
        OPEN(11,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',     &
                           ACCESS='STREAM')
    
        PRINT *, "read: ", fname
        !Skip header
        READ(11)   cbuffer
        READ(11)   cbuffer
        READ(11)   int_in
        READ(11)   cbuffer
        !Skip ttr values before nres
        DO j=1,nres
           READ(11) ttr_in
        END DO
        
        ttr=ttr_in
        WRITE(10,'(E16.8)') ttr
        CLOSE(11)

      CASE('ascii')
        fname=job_name(1:INDEX(job_name, " ")-1) // ".ensi.NDTTR-"//stepnum
        OPEN(11,file=fname,status='old',action='read')
        PRINT *, "read: ", fname
        
        DO j=1,4+nres-1      !Skip header and ttr values before nres
          READ (11,*)
        END DO
        READ (11,*) ttr
        
        WRITE(10,'(E16.8)') ttr
        CLOSE(11)
        
      CASE default
        PRINT*, "Invalid format flag"
        PRINT*, "Program aborting"
    END SELECT
    
  END DO
    
  CLOSE(10)

END PROGRAM ndttrget
