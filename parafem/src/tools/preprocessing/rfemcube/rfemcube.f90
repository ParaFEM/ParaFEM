PROGRAM rfemcube

!/****h* tools/preprocessing/rfemcube
!*  NAME
!*    PROGRAM rfemcube
!*  FUNCTION
!*    Creates box of finite elements
!*  AUTHOR
!*    Louise M. Lever
!*  COPYRIGHT
!*    (c) University of Manchester 2012
!****
!*/

  USE PRECISION
  USE GEOMETRY
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER                :: argc,iargc
  INTEGER                :: iel,i,j,k
  INTEGER                :: nels,nxe,nye,nze,ndim,nod,nn,nip
  INTEGER                :: output
  REAL(iwp)              :: aa,bb,cc
  REAL(iwp)              :: eavg,egeo,ehrm
  REAL(iwp)              :: emn,esd
  REAL(iwp)              :: thx,thy,thz
  REAL(iwp)              :: v
  REAL(iwp),PARAMETER    :: zero=0.0_iwp
  REAL(iwp)              :: tol=1.000000e-06,mises=70.0
  INTEGER                :: limit=2000
  CHARACTER(LEN=6)       :: varfnc
  CHARACTER(LEN=15)      :: rfield
  CHARACTER(LEN=50)      :: in_rfem_name, out_model_name,fname

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------   
  INTEGER, ALLOCATABLE   :: num(:),g_num(:,:) 
  REAL(iwp), ALLOCATABLE :: coord(:,:),g_coord(:,:)

!------------------------------------------------------------------------------
! 3. Read in_rfem_name and out_model_name from the command line
!------------------------------------------------------------------------------

  argc = iargc()
  IF( (argc /= 2) ) THEN
     PRINT*
     PRINT*, "Usage:  rfemcube <in_rfem_name> <out_model_name>"
     PRINT*
     PRINT*, "        program expects as input:"
     PRINT*, "          <in_rfem_name>.rf"
     PRINT*
     PRINT*, "        and outputs:"
     PRINT*, "          <out_model_name>.d" 
     PRINT*, "          <out_model_name>.dat" 
     PRINT*, "          <in_rfem_name>.res" 
     PRINT*
     STOP
  END IF
  CALL GETARG(1, in_rfem_name)
  CALL GETARG(2, out_model_name)
  
!------------------------------------------------------------------------------
! 4. Read control data
!------------------------------------------------------------------------------

  IF (INDEX(in_rfem_name,".rf") /= 0) THEN
     in_rfem_name = in_rfem_name(1:INDEX(in_rfem_name,".rf")-1)
  END IF

  fname = in_rfem_name(1:INDEX(in_rfem_name," ")-1) // ".rf"
  OPEN (10, file=fname, status='old', action='read')
  
  READ(10,*) rfield
  
!------------------------------------------------------------------------------
! 5. Read based on selected random field generator
!------------------------------------------------------------------------------

  SELECT CASE(rfield)
     
  CASE('sim3de')
     
     READ(10,*) output
     READ(10,*) nels,nxe,nze
     READ(10,*) aa,bb,cc
     READ(10,*) thx,thy,thz
     READ(10,*) emn,esd
     READ(10,*) varfnc
     READ(10,*) v
     
     nye    = nels/nxe/nze
     nn     = (nxe+1)*(nye+1)*(nze+1)
     ndim   = 3
     nip    = 8
     nod    = 8
     
     fname = in_rfem_name(1:INDEX(in_rfem_name," ")-1) // ".res"
     OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

!------------------------------------------------------------------------------
! 8. Generate finite element mesh for a regular cuboid
!------------------------------------------------------------------------------
    
     ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),num(nod),g_num(nod,nels))
     
     coord = zero ; g_coord = zero
     num   = 0    ; g_num   = 0
     
     DO iel = 1, nels
        CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
        g_coord(:, g_num(:,iel)) = TRANSPOSE(coord)
     END DO
     
     fname = out_model_name(1:INDEX(out_model_name," ")-1) // ".d"
     OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
     
     WRITE(13,'(A)') "*THREE_DIMENSIONAL"
     WRITE(13,'(A)') "*NODES"
     
     DO i = 1,nn
        WRITE(13,'(I8,3E14.6)') i, g_coord(:,i)
     END DO
     
     WRITE(13,'(A)') "*ELEMENTS"
     
     DO iel = 1, nels
        WRITE(13,'(I8,A,8I8,I8)') iel, " 3 8 1 ",g_num(1,iel),g_num(4,iel),   &
             g_num(8,iel),g_num(5,iel),g_num(2,iel),   &
             g_num(3,iel),g_num(7,iel),g_num(6,iel),   &
             iel  ! each element has its own material  
     END DO
     
     CLOSE(13)
     
     DEALLOCATE(coord,g_coord,num,g_num)
     
!------------------------------------------------------------------------------
! 10. Write .dat file for RFEMBC
!------------------------------------------------------------------------------

     fname = out_model_name(1:INDEX(out_model_name," ")-1) // ".dat"
     OPEN (14, file=fname, status='REPLACE', action='WRITE')
     WRITE(14,*) "hexahedron"
     WRITE(14,*) "1"
     WRITE(14,*) "1"
     WRITE(14,*) nels
     WRITE(14, '(7I8)') nels, nn, 0, nip, nod, 0, 0
     WRITE(14, '(E14.6,I12, E14.6)') tol, limit, mises
     CLOSE(14)

!------------------------------------------------------------------------------
! 10. Cleanup
!------------------------------------------------------------------------------

  CASE DEFAULT
     
     PRINT *
     PRINT *, "  *** Program aborted ***"
     PRINT *
     PRINT *, "  Wrong value given in variable: FIELD_TYPE"
     PRINT *, "  The only accepted value is:    SIM3DE"
     PRINT *
     
  END SELECT
  
  CLOSE(10)
  CLOSE(11)
  
END PROGRAM rfemcube
