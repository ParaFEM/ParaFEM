PROGRAM ensibin2ascii
!------------------------------------------------------------------------------
!      Program ensibin2ascii  Tool to convert binary ensight gold formatted  
!                             variable files into ascii. Under development,
!                             currently for 4-node tetrahedra
!
! Usage: ensibin2ascii <job_name> <file extension> <number of variables> <ndim>
!
!------------------------------------------------------------------------------
  
  USE, INTRINSIC :: ISO_C_BINDING ! to output C binary file
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------
  
  INTEGER,PARAMETER    :: iwp=SELECTED_REAL_KIND(15)
  CHARACTER(LEN=80)    :: fname,job_name,ext,arg,cbuffer
  INTEGER              :: nvar,ndim,i,j
  INTEGER(KIND=C_INT)  :: int_in
  REAL(KIND=C_FLOAT)   :: var_in
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------
  

  
!------------------------------------------------------------------------------
! 3. Read file_name from the command line.
!    Read control data
!------------------------------------------------------------------------------
  
  CALL GETARG(1,job_name)
  CALL GETARG(2,ext)
  CALL GETARG(3,arg)
  read(arg,*) nvar !str2int
  CALL GETARG(4,arg)
  read(arg,*) ndim !str2int
    
!------------------------------------------------------------------------------
! 4. Open files for reading/writing
!------------------------------------------------------------------------------
  
  PRINT *, "Starting conversion of binary ENSI variable file into ascii"

  fname = job_name(1:INDEX(job_name, " ") -1) // ".bin.ensi."//ext
  OPEN(10,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',           &
                     ACCESS='STREAM')
  
  fname = job_name(1:INDEX(job_name, " ") -1) // ".ensi."//ext
  OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
  READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
  READ(10) int_in  ; WRITE(11,'(I0)')  int_in
  READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
  
  DO i=1,nvar
    DO j=1,ndim
        READ(10) var_in  ; WRITE(11,'(E16.9)',ADVANCE='no')  var_in
        IF(j == ndim) WRITE(11,*)''
    END DO
  END DO
  
  CLOSE(10)
  CLOSE(11)
  PRINT *, "Conversion complete"
  
!  CALL shutdown()
  
END PROGRAM ensibin2ascii
