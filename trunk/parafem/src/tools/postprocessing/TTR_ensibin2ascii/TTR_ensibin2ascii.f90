PROGRAM TTR_ensibin2ascii
!------------------------------------------------------------------------------
!      Program TTR_ensibin2ascii  Tool to convert binary ensight gold formatted  
!                                 scalar result files into ascii. Under 
!                                 development, currently for 4-node tetrahedra
!------------------------------------------------------------------------------
  
  USE, INTRINSIC :: ISO_C_BINDING ! to output C binary file
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------
  
  INTEGER,PARAMETER   :: iwp=SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER  :: ndim=3,nodof=1,nprops=5
  INTEGER             :: meshgen,partitioner,np_types,nels,nn,nr,nip,nod
  INTEGER             :: loaded_nodes,fixed_freedoms,nstep,npri
  INTEGER             :: limit,el_print,i_o
  INTEGER             :: i,j,k,prnwidth
  INTEGER(KIND=C_INT) :: int_in
  REAL(KIND=C_FLOAT)  :: ttr_in
  REAL(iwp),PARAMETER :: zero = 0.0_iwp
  REAL(iwp)           :: val0,dtim,theta,tol,etype,real_time,ttr
  CHARACTER(LEN=50)   :: fname,job_name,stepnum
  CHARACTER(LEN=15)   :: element
  CHARACTER(LEN=80)   :: cbuffer
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------
  
  REAL(iwp),ALLOCATABLE :: temp_real(:),prop(:,:),val_f(:),val(:,:)
  INTEGER,ALLOCATABLE   :: temp_int(:),node(:),sense(:)
  
!------------------------------------------------------------------------------
! 3. Read job_name from the command line.
!    Read control data
!------------------------------------------------------------------------------
  
  CALL GETARG(1,job_name)
  
  PRINT *, "Starting conversion of ENSI result files"

  fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
  OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
  READ(10,*) element,meshgen,partitioner,np_types,                            &
             nels,nn,nr,nip,nod,loaded_nodes,fixed_freedoms,                  &
             val0,                                                            &
             dtim,nstep,npri,theta,                                           &
             tol,limit,el_print,i_o
  CLOSE(10)
  
!----------------------------------------------------------------------------
! 4. Open files
!----------------------------------------------------------------------------
  
  k=nstep/npri+1
  DO j=1,k
    WRITE(stepnum,'(I0.6)') j
    
    fname   = job_name(1:INDEX(job_name, " ")-1) // ".bin.ensi.NDTTR-"//stepnum
    OPEN(12,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',         &
                       ACCESS='STREAM')
    
    fname   = job_name(1:INDEX(job_name, " ")-1)//".ensi.NDTTR-"//stepnum
    OPEN(13, file=fname, status='replace', action='write')
    
!----------------------------------------------------------------------------
! 5. Read/Write variable files
!----------------------------------------------------------------------------
    
    READ(12)   cbuffer
    READ(12)   cbuffer
    READ(12)   int_in
    READ(12)   cbuffer
    
    WRITE(13,'(A)')   "Alya Ensight Gold --- Scalar per-node variable file"
    WRITE(13,'(A)')   "part"
    WRITE(13, '(I12)') 1
    WRITE(13,'(A)')   "coordinates"
    
    PRINT *, "Reading/Writing results for file",j
    DO i = 1,nn
      READ(12) ttr_in
      ttr = ttr_in
      WRITE(13, '(E15.9)') ttr
    END DO
    
    CLOSE(12);CLOSE(13)
  END DO
  
  PRINT *, "Conversion complete"
  
  CALL shutdown()
  
END PROGRAM TTR_ensibin2ascii
