PROGRAM endian_convert     
!------------------------------------------------------------------------------
!      Program endian_convert
!
!------------------------------------------------------------------------------
  
! Usage: endian_convert <binary/ascii> <variable/geometry> <job_name>
!                       <file extension> <number of variables> <ndim>
  
  USE, INTRINSIC :: ISO_C_BINDING ! to output C binary file
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------
  
  INTEGER,PARAMETER    :: iwp=SELECTED_REAL_KIND(15)
  CHARACTER(LEN=80)    :: fname,bin_asc,var_geo,job_name,ext,arg,cbuffer
  INTEGER              :: nvar,ndim,nels,i,j
  INTEGER(KIND=C_INT)  :: int_in
  REAL(KIND=C_FLOAT)   :: var_in
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------ 
  
!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions
!------------------------------------------------------------------------------
  
  CALL GETARG(1,bin_asc)
  CALL GETARG(2,var_geo)
  CALL GETARG(3,job_name)
  CALL GETARG(4,ext)
  CALL GETARG(5,arg)
  read(arg,*) nvar !str2int
  CALL GETARG(6,arg)
  read(arg,*) ndim !str2int
  
!------------------------------------------------------------------------------
! 4. Open files for reading/writing
!------------------------------------------------------------------------------
  
  PRINT *, "Starting endian conversion of file"
  
  fname = job_name(1:INDEX(job_name, " ") -1) // ".bin.ensi."//ext
! GWL: SWAP is not a valid CONVERT type in ifort. So we use 'NATIVE' instead
! GWL: and rely on the user setting the env var FORT_CONVERT10=LITTLE_ENDIAN or
! GWL: BIG_ENDIAN at runtime to override the 'NATIVE' setting.
!  OPEN(10,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',           &
!                     ACCESS='STREAM',CONVERT='SWAP')
  OPEN(10,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',           &
                     ACCESS='STREAM',CONVERT='NATIVE')
  
  SELECT CASE (var_geo)
    CASE('var')
      SELECT CASE (bin_asc)
        CASE('binary')
          fname = job_name(1:INDEX(job_name, " ") -1) // ".bin.ensi.endianswap."//ext
          OPEN(11,FILE=fname,STATUS='REPLACE',FORM='UNFORMATTED',ACTION='WRITE',  &
                             ACCESS='STREAM')
          
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) int_in  ; WRITE(11) int(int_in,kind=c_int)
          READ(10) cbuffer ; WRITE(11) cbuffer
          
          DO i=1,nvar
            DO j=1,ndim
                READ(10) var_in
                WRITE(11) real(var_in,kind=c_float)
!                IF(j == ndim) WRITE(11,*)''
            END DO
          END DO
          PRINT *, "Conversion complete"
          
        CASE('ascii')
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
          PRINT *, "Conversion complete"
          
        CASE default
          PRINT*, "Invalid binary/ascii flag"
          PRINT*, "Program aborting"
      END SELECT
    
    CASE('geo')
      SELECT CASE (bin_asc)
        CASE('binary')
          fname = job_name(1:INDEX(job_name, " ") -1) // ".bin.ensi.endianswap."//ext
          OPEN(11,FILE=fname,STATUS='REPLACE',FORM='UNFORMATTED',ACTION='WRITE',  &
                             ACCESS='STREAM')
          
          !-Read/Write .geo header
          PRINT *, "Reading and writing .geo file header"
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) int_in  ; WRITE(11) int(int_in,kind=c_int)
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) int_in  ; WRITE(11) int(int_in,kind=c_int)
          nvar=int_in
          
          !-Read/Write nodal coordinates 
          PRINT *, "Reading and writing nodal coordinates"
          DO j=1,ndim
            DO i=1,nvar
              READ(10) var_in ; WRITE(11) real(var_in,kind=c_float)
            END DO
          END DO
          
          !-Read/Write elements
          PRINT *, "Reading and writing elements"
          READ(10) cbuffer ; WRITE(11) cbuffer
          READ(10) int_in  ; WRITE(11) int(int_in,kind=c_int)
          nels=int_in
          DO j=1,nels
            DO i=1,4
              READ(10) int_in  ; WRITE(11) int(int_in,kind=c_int)
            END DO
          END DO
          PRINT *, "Conversion complete"
          
        CASE('ascii')
          fname = job_name(1:INDEX(job_name, " ") -1) // ".ensi."//ext
          OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
          
          !-Read/Write .geo header
          PRINT *, "Reading and writing .geo file header"
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) int_in  ; WRITE(11,'(I0)')  int_in
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) int_in  ; WRITE(11,'(I0)')  int_in
          
          !-Read/Write nodal coordinates 
          PRINT *, "Reading and writing nodal coordinates"
          DO j=1,ndim
            DO i=1,nvar
!              READ(10) var_in  ; WRITE(11,'(E19.13)')  var_in
              READ(10) var_in  ; WRITE(11,'(ES19.12E3)')  var_in
            END DO
          END DO
          
          !-Read/Write elements
          PRINT *, "Reading and writing elements"
          READ(10) cbuffer ; WRITE(11,'(A)')   cbuffer
          READ(10) int_in  ; WRITE(11,'(I0)')  int_in
          nels=int_in
          DO j=1,nels
            DO i=1,4
              READ(10) int_in  ; WRITE(11,'(I0,A)',ADVANCE='no')  int_in,' '
              IF(i == 4) WRITE(11,*)''
            END DO
          END DO
          PRINT *, "Conversion complete"
          
        CASE default
          PRINT*, "Invalid binary/ascii flag"
          PRINT*, "Program aborting"
      END SELECT
    END SELECT
  
  CLOSE(10)
  CLOSE(11)
  
  !CALL shutdown()
  
END PROGRAM endian_convert
