PROGRAM rfembc

!/****h* tools/preprocessing/rfembc
!*  NAME
!*    PROGRAM rfembc
!*  FUNCTION
!*    Loads an unstructured mesh and finds zmin and zmax nodes. Generates a
!*    .bnd file for zmin and a .fix file for zmax
!*  AUTHORS
!*    Louise M. Lever, Lee Margetts
!*  COPYRIGHT
!*    (c) University of Manchester 2012
!****
!*/

  USE PRECISION
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER                :: argc,iargc
  INTEGER                :: ndim, i
  CHARACTER(LEN=50)      :: job_name,res_name,in_dat_name,out_dat_name,d_name
  CHARACTER(LEN=50)      :: bnd_name,fix_name
  CHARACTER(LEN=15)      :: fixed_value_arg
  CHARACTER(LEN=15)      :: element
  INTEGER                :: nels,nn,nr,nod,nip,loaded_nodes
  INTEGER                :: limit,mesh,fixed_freedoms,partition,np_types=1
  REAL(iwp)              :: e,v,tol,mises=70.0
  REAL(iwp)              :: fixed_value
  INTEGER                :: idummy
  REAL(iwp)              :: rdummy
  INTEGER                :: zord=3
  REAL(iwp)              :: zmin_extent,zmax_extent
  INTEGER                :: num_bnd=0,num_fix=0

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------   

  REAL(iwp), ALLOCATABLE :: m_coord(:,:)

!------------------------------------------------------------------------------
! 3. Read job_name from the command line
!------------------------------------------------------------------------------

  argc = iargc()
  IF (argc /= 2) THEN
    PRINT*
    PRINT*, "Usage:  rfembc <job_name> <fixed_value>"
    PRINT*
    PRINT*, "        program expects <job_name>.d and outputs"
    PRINT*, "        <job_name>-rfem.dat,"
    PRINT*, "        <job_name>-rfem.bnd and <job_name>-rfem.fix"
    PRINT*
    STOP
  END IF
  CALL GETARG(1, job_name)
  CALL GETARG(2, fixed_value_arg)
  READ(fixed_value_arg,*) fixed_value

!------------------------------------------------------------------------------
! 4. Open result file
!------------------------------------------------------------------------------
     
  res_name = job_name(1:INDEX(job_name," ")-1) // ".res"
  OPEN(11,FILE=res_name,STATUS='REPLACE',ACTION='WRITE')
  
!------------------------------------------------------------------------------
! 5. Read model dat, coords and elems; determine extents
!------------------------------------------------------------------------------
     
  IF (INDEX(job_name,".dat") /= 0) THEN
     job_name = job_name(1:INDEX(job_name,".dat")-1)
  END IF
  in_dat_name = job_name(1:INDEX(job_name," ")-1) // ".dat"
  out_dat_name = job_name(1:INDEX(job_name," ")-1) // "-rfem.dat"
  d_name = job_name(1:INDEX(job_name," ")-1) // ".d"

  OPEN (14, file=in_dat_name, status='old', action='read')
  READ(14,*) element,mesh,partition,nels,nn,nr,nip,nod,loaded_nodes,        &
             fixed_freedoms,e,v,tol,limit
  CLOSE(14)
  
  ndim = 3
  ALLOCATE(m_coord(ndim,nn))
  
  OPEN (15, file=d_name, status='old', action='read')
  READ(15,*)   !headers
  READ(15,*)   !headers

  READ(15,*) idummy,m_coord(:,1)
  zmin_extent = m_coord(zord,1)
  zmax_extent = m_coord(zord,1)
  DO i = 2,nn
     READ(15,*) idummy,m_coord(:,i)
     IF ( m_coord(zord,i) < zmin_extent ) THEN
        zmin_extent = m_coord(zord,i)
     END IF
     IF ( m_coord(zord,i) > zmax_extent ) THEN
        zmax_extent = m_coord(zord,i)
     END IF
  END DO

  CLOSE(15)
  
  WRITE(11, '(A,I12)') "Number of model nodes = ", nn
  WRITE(11, '(A,I12)') "Number of model elems = ", nels
  WRITE(11, '(A,I12)') "Number of model nods = ", nod
  
  WRITE(11, '(A,2F)') "ZMIN/ZMAX of model: ", zmin_extent, zmax_extent
  
!------------------------------------------------------------------------------
! 6. Write .bnd and .fix files where nodes match zmin/zmax respectively
!------------------------------------------------------------------------------

  bnd_name = job_name(1:INDEX(job_name," ")-1) // "-rfem.bnd"
  fix_name = job_name(1:INDEX(job_name," ")-1) // "-rfem.fix"

  OPEN(16,FILE=bnd_name,STATUS='REPLACE',ACTION='WRITE')
  OPEN(17,FILE=fix_name,STATUS='REPLACE',ACTION='WRITE')

  DO i = 1,nn
     IF( m_coord(zord,i) == zmin_extent ) THEN
        WRITE(16,'(I12,A)') i, " 0 0 0"
        num_bnd = num_bnd + 1
     END IF
     IF( m_coord(zord,i) == zmax_extent ) THEN
        WRITE(17,'(I12,A,E14.6)') i, " 3", fixed_value
        num_fix = num_fix + 1
     END IF
  END DO
  
  CLOSE(16)
  CLOSE(17)

  WRITE(11, '(A,I12)') "Number of bound nodes = ", num_bnd
  WRITE(11, '(A,I12)') "Number of fixed elems = ", num_fix

  DEALLOCATE(m_coord)

!------------------------------------------------------------------------------
! 7. Rewrite .dat file with updated np_types(1), nr(num_bnd),
!    fixed_freedoms(num_fix) and mises(70.0)
!------------------------------------------------------------------------------

  OPEN (14, file=out_dat_name, status='REPLACE', action='WRITE')
  WRITE(14,*) element
  WRITE(14,*) mesh
  WRITE(14,*) partition
  WRITE(14,*) np_types
  WRITE(14, '(6I)') nels, nn, num_bnd, nip, nod, loaded_nodes, num_fix
  WRITE(14, '(E14.6,I, E14.6)') tol, limit, mises
  CLOSE(14)
  
!------------------------------------------------------------------------------
! 8. Close result file
!------------------------------------------------------------------------------
     
  CLOSE(11)

END PROGRAM rfembc
