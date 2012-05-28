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

  USE PRECISION
  USE GEOMETRY

  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER                :: argc,iargc
  INTEGER                :: iel,i
  INTEGER                :: nels,nxe,nye,nze,ndim,nod,nn
  INTEGER                :: iefld,istat,kseed,output
  REAL(iwp)              :: aa,bb,cc
  REAL(iwp)              :: eavg,egeo,ehrm
  REAL(iwp)              :: emn,esd
  REAL(iwp)              :: thx,thy,thz
  REAL(iwp),PARAMETER    :: zero=0.0_iwp
  CHARACTER(LEN=6)       :: varfnc
  CHARACTER(LEN=15)      :: rfield,job,sub1,sub2
  CHARACTER(LEN=50)      :: job_name,fname,program_name
  LOGICAL                :: debug=.false.
  LOGICAL                :: shofld=.false.
  LOGICAL                :: lunif=.false.
  LOGICAL                :: dcheck=.false.

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------   
  INTEGER, ALLOCATABLE   :: num(:),g_num(:,:) 
  INTEGER, ALLOCATABLE   :: ieplt(:)
  REAL(iwp), ALLOCATABLE :: coord(:,:),g_coord(:,:)
  REAL(iwp), ALLOCATABLE :: efld(:)   ! elastic modulus random field

!------------------------------------------------------------------------------
! 3. Read job_name from the command line
!------------------------------------------------------------------------------

  argc = iargc()
  IF (argc /= 1) THEN
    PRINT*
    PRINT*, "Usage:  rfem <job_name>"
    PRINT*
    PRINT*, "        program expects <job_name>.dat and outputs"
    PRINT*, "        <job_name>.d" 
    PRINT*
    STOP
  END IF
  CALL GETARG(1, job_name)

!------------------------------------------------------------------------------
! 4. Read control data
!------------------------------------------------------------------------------

  IF (INDEX(job_name,".dat") /= 0) THEN
    job_name = job_name(1:INDEX(job_name,".dat")-1)
  END IF
  fname = job_name(1:INDEX(job_name," ")-1) // ".dat"
  OPEN (10, file=fname, status='old', action='read')

  READ(10,*) rfield

!------------------------------------------------------------------------------
! 5. Select random field generator
!------------------------------------------------------------------------------

  SELECT CASE(rfield)

  CASE('sim3de')

    READ(10,*) output
    READ(10,*) nels,nxe,nze
    READ(10,*) aa,bb,cc
    READ(10,*) thx,thy,thz
    READ(10,*) emn,esd
    READ(10,*) varfnc

    nye  = nels/nxe/nze
    nn   = (nxe+1)*(nye+1)*(nze+1)
    ndim = 3
    nod  = 8

    IF(output /= 1 .AND. output /=2) THEN
      PRINT *, output, " = Incorrect value for output. Accepted values are &
                           &'1' or '2'"
      STOP
    END IF

!------------------------------------------------------------------------------
! 6. Generate spatially random field for Young's modulus in a regular cuboid
!------------------------------------------------------------------------------

   ALLOCATE(efld(nels))
   ALLOCATE(ieplt(3))

   lunif  = .false.  ! not used
   shofld = .false.  ! not used 
   ieplt  = 0        ! not used
   job    = ""       ! not used
   sub1   = ""       ! not used
   sub2   = ""       ! not used 
   debug  = .false.  ! not used
   dcheck = .false.  ! not used

   istat  = 11       ! output debugging info to .res file
   efld   = zero
   kseed  = 0        ! important variable

   fname = job_name(1:INDEX(job_name," ")-1) // ".res"
   OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

   fname = job_name(1:INDEX(job_name," ")-1) // ".efld"
   OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

   CALL sim3de(efld,nxe,nye,nze,aa,bb,cc,thx,thy,thz,emn,esd,lunif,kseed, &
               shofld,ieplt,iefld,job,sub1,sub2,varfnc,debug,istat,dcheck,   &
               nxe,nye,nze,eavg,egeo,ehrm)

   DO iel = 1,nels
     WRITE(12,'(I10,E12.4)') iel, efld(iel)
   END DO

   CLOSE(11)
   CLOSE(12)

!------------------------------------------------------------------------------
! 7. Generate finite element mesh for a regular cuboid
!------------------------------------------------------------------------------


  IF(output == 2 ) THEN

    ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),num(nod),g_num(nod,nels))

    coord = zero ; g_coord = zero
    num   = 0    ; g_num   = 0

    DO iel = 1, nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:, g_num(:,iel)) = TRANSPOSE(coord)
    END DO

    fname = job_name(1:INDEX(job_name," ")-1) // ".d"
    OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

    WRITE(13,'(A)') "*THREE_DIMENSIONAL"
    WRITE(13,'(A)') "*NODES"
  
    DO i = 1,nn
      WRITE(13,'(I12,3E14.6)') i, g_coord(:,i)
    END DO

    WRITE(13,'(A)') "*ELEMENTS"

    DO iel = 1, nels
      WRITE(13,'(I12,A,8I12,A)') iel, " 3 8 1 ",g_num(1,iel),g_num(4,iel),   &
                                   g_num(8,iel),g_num(5,iel),g_num(2,iel),   &
                                   g_num(3,iel),g_num(7,iel),g_num(6,iel),   &
                                   " 1"   
    END DO

    CLOSE(13)

  END IF

  CASE DEFAULT
  
    PRINT *
    PRINT *, "  *** Program aborted ***"
    PRINT *
    PRINT *, "  Wrong value given in variable: FIELD_TYPE"
    PRINT *, "  The only accepted value is:    SIM3DE"
    PRINT *
      
  END SELECT
  

END PROGRAM rfem
