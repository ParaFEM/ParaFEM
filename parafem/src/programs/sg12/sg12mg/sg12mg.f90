PROGRAM sg12mg

!/****h* Programs/Smith&Griffiths/sg12mg
!*  NAME
!*    PROGRAM sg12mg
!*  FUNCTION
!*    Creates simple meshes for the programs in Smith and Griffiths Edition 4.
!*  AUTHOR
!*    Lee Margetts
!*  COPYRIGHT
!*    (c) University of Manchester 2007-2010
!****
!*/

  USE precision
  USE geometry
  USE loading
  
  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER                :: nr
  INTEGER                :: nn
  INTEGER                :: nels
  INTEGER                :: nle
  INTEGER                :: nxe
  INTEGER                :: nye
  INTEGER                :: nze 
  INTEGER                :: nod 
  INTEGER                :: ndim
  INTEGER                :: nodof
  INTEGER                :: nip
  INTEGER                :: nmodes
  INTEGER                :: lalfa
  INTEGER                :: leig
  INTEGER                :: lx
  INTEGER                :: lz
  INTEGER                :: i
  INTEGER                :: iel
  INTEGER                :: argc
  INTEGER                :: cjits
  INTEGER                :: ell
  INTEGER                :: fixed_nodes
  INTEGER                :: iargc
  INTEGER                :: limit
  INTEGER                :: loaded_freedoms
  INTEGER                :: fixed_freedoms
  INTEGER                :: nev
  INTEGER                :: ncv
  INTEGER                :: maxitr
  INTEGER                :: meshgen
  REAL(iwp)              :: aa
  REAL(iwp)              :: bb 
  REAL(iwp)              :: cc  
  REAL(iwp)              :: cjtol  
  REAL(iwp)              :: rho   
  REAL(iwp)              :: e 
  REAL(iwp)              :: v  
  REAL(iwp)              :: visc
  REAL(iwp)              :: el  
  REAL(iwp)              :: er  
  REAL(iwp)              :: acc  
  REAL(iwp)              :: kappa  
  REAL(iwp)              :: penalty  
  REAL(iwp)              :: tol 
  REAL(iwp)              :: x0 
  CHARACTER(LEN=15)      :: element
  CHARACTER(LEN=50)      :: program_name
  CHARACTER(LEN=50)      :: job_name
  CHARACTER(LEN=50)      :: fname
  CHARACTER(LEN=1)       :: bmat
  CHARACTER(LEN=2)       :: which

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------   

  INTEGER, ALLOCATABLE   :: g_num(:,:)
  INTEGER, ALLOCATABLE   :: rest(:,:)
  INTEGER, ALLOCATABLE   :: nf(:,:)
  INTEGER, ALLOCATABLE   :: no(:)
  REAL(iwp), ALLOCATABLE :: g_coord(:,:)
  REAL(iwp), ALLOCATABLE :: coord(:,:)
  REAL(iwp), ALLOCATABLE :: val(:)

!------------------------------------------------------------------------------
! 3. Read job_name from the command line
!------------------------------------------------------------------------------

  argc = iargc()
  IF (argc /= 1) THEN
    PRINT*
    PRINT*, "Usage:  sg12mg <job_name>"
    PRINT*
    PRINT*, "        program expects <job_name>.mg and outputs any or all of"
    PRINT*, "        <job_name>.dat <job_name>.d" 
    PRINT*, "        <job_name>.bnd <job_name>.lds"
    PRINT*
    STOP
  END IF
  CALL GETARG(1, job_name)

!------------------------------------------------------------------------------
! 4. Read control data
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name," ")-1) // ".mg"
  OPEN (10, file=fname, status='old', action='read')

  READ(10,*) program_name
  
!------------------------------------------------------------------------------
! 6. Select program using the SELECT CASE CONSTRUCT
!------------------------------------------------------------------------------

  SELECT CASE(program_name)

!------------------------------------------------------------------------------
! 7. Program p121
!------------------------------------------------------------------------------

  CASE('p121')

  READ(10,*) nels, nxe, nze, nip
  READ(10,*) aa, bb, cc, e, v
  READ(10,*) tol, limit

  nye   = nels/nxe/nze
  nr    = ((2*nxe+1)*(nze+1)+(nxe+1)*nze)*2 +                                &
          ((2*nye-1)*nze+(nye-1)*nze)*2     +                                &
          (2*nye-1)*(nxe+1)                 +                                &
          (nye-1)*nxe
  nle   = nxe/5
  ndim  = 3
  nodof = 3
  nod   = 20
  nn    = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1) +                      &
          (nxe+1)*(nze+1)*nye

  loaded_freedoms = 3*nle*nle + 4*nle + 1

!------------------------------------------------------------------------------
! 7.2 Allocate dynamic arrays
!------------------------------------------------------------------------------

  ALLOCATE(coord(nod,ndim))
  ALLOCATE(g_coord(ndim,nn))
  ALLOCATE(g_num(nod,nels))
  ALLOCATE(rest(nr,nodof+1))
  ALLOCATE(val(loaded_freedoms))
  ALLOCATE(no(loaded_freedoms))
  
  coord    = 0.0_iwp
  g_coord  = 0.0_iwp
  g_num    = 0
  rest     = 0
  val      = 0.0_iwp
  no       = 0

!------------------------------------------------------------------------------
! 7.3 Find nodal coordinates and element steering array
!     Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

  DO iel = 1, nels
    CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
    g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
  END DO

  fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
  OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  WRITE(11,'(A)') "*THREE_DIMENSIONAL"
  WRITE(11,'(A)') "*NODES"

  DO i = 1, nn
    WRITE(11,'(I12,3F12.4)') i, g_coord(:,i)
  END DO

  WRITE(11,'(A)') "*ELEMENTS"
  DO iel = 1, nels
    WRITE(11,'(I12,A,20I12,A)') iel, " 3 20 1 ", g_num(3,iel),g_num(5,iel),   &
                                 g_num(7,iel),g_num(1,iel),g_num(15,iel),     &
                                 g_num(17,iel),g_num(19,iel),g_num(13,iel),   &
                                 g_num(4,iel),g_num(6,iel),g_num(8,iel),      &
                                 g_num(2,iel),g_num(16,iel),g_num(18,iel),    &
                                 g_num(20,iel),g_num(14,iel),g_num(10,iel),   &
                                 g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1 "
  END DO

  CLOSE(11)

!------------------------------------------------------------------------------
! 7.4 Boundary conditions
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
  OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  CALL cube_bc20(rest,nxe,nye,nze)
 
  DO i = 1, nr
    WRITE(12,'(I8,3I6)') rest(i,:) 
  END DO

  CLOSE(12)

!------------------------------------------------------------------------------
! 7.5 Loading conditions
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
  OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  CALL load_p121(no,val,nle,nxe,nze)
  val = -val * aa * bb / 12._iwp
 
  DO i = 1, loaded_freedoms
    WRITE(13,'(I10,2A,3E12.4)') no(i), "  0.0000E+00  ", "0.0000E+00", val(i) 
  END DO

  CLOSE(13)

!------------------------------------------------------------------------------
! 7.6 New control data
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
  OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  WRITE(14,'(A)') "hexahedron"
  WRITE(14,'(A)') "2"              ! Abaqus node numbering scheme
  WRITE(14,'(6I9)') nels, nn, nr, nip, nod, loaded_freedoms
  WRITE(14,'(3E12.4,I8)') e, v, tol, limit

  CLOSE(14)

!------------------------------------------------------------------------------
! 8. Program p122
!------------------------------------------------------------------------------

  CASE('p122')

  PRINT*
  PRINT*, "Program p122 not supported"
  PRINT*

  STOP

!------------------------------------------------------------------------------
! 9. Program p123
!------------------------------------------------------------------------------

  CASE('p123')

  PRINT*
  PRINT*, "Program p123 not supported"
  PRINT*

  STOP

!------------------------------------------------------------------------------
! 10. Program p124
!------------------------------------------------------------------------------

  CASE('p124')

  PRINT*
  PRINT*, "Program p124 not supported"
  PRINT*

  STOP

!------------------------------------------------------------------------------
! 11. Program p125
!------------------------------------------------------------------------------

  CASE('p125')

  PRINT*
  PRINT*, "Program p125 not supported"
  PRINT*

  STOP

!------------------------------------------------------------------------------
! 12. Program p126
!------------------------------------------------------------------------------

  CASE('p126')

  READ(10,*) nels, nxe, nze, nip
  READ(10,*) aa, bb, cc
  READ(10,*) visc, rho, tol, limit
  READ(10,*) cjtol, cjits, penalty
  READ(10,*) x0, ell, kappa
  
  nye            = nels/nxe/nze
  fixed_freedoms = 3*nxe*nye+2*nxe+2*nye+1
  nr             = 3*nxe*nye*nze+4*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+2
  ndim           = 3
  nodof          = 3
  nod            = 20
  nn             = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1) +              &
                   (nxe+1)*(nze+1)*nye

!------------------------------------------------------------------------------
! 12.1 Allocate dynamic arrays
!------------------------------------------------------------------------------

  ALLOCATE(coord(nod,ndim))
  ALLOCATE(g_coord(ndim,nn))
  ALLOCATE(g_num(nod,nels))
  ALLOCATE(rest(nr,nodof+1))
  ALLOCATE(val(fixed_freedoms))
  ALLOCATE(no(fixed_freedoms))
  
  coord    = 0.0_iwp
  g_coord  = 0.0_iwp
  g_num    = 0
  rest     = 0
  val      = 0.0_iwp
  no       = 0

!------------------------------------------------------------------------------
! 12.2 Find nodal coordinates and element steering array
!     Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

  DO iel = 1, nels
    CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
    g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
  END DO

  fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
  OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  WRITE(11,'(A)') "*THREE_DIMENSIONAL"
  WRITE(11,'(A)') "*NODES"

  DO i = 1, nn
    WRITE(11,'(I12,3F12.4)') i, g_coord(:,i)
  END DO

  WRITE(11,'(A)') "*ELEMENTS"
  DO iel = 1, nels
    WRITE(11,'(I12,A,20I12,A)') iel, " 3 20 1 ", g_num(3,iel),g_num(5,iel),   &
                                 g_num(7,iel),g_num(1,iel),g_num(15,iel),     &
                                 g_num(17,iel),g_num(19,iel),g_num(13,iel),   &
                                 g_num(4,iel),g_num(6,iel),g_num(8,iel),      &
                                 g_num(2,iel),g_num(16,iel),g_num(18,iel),    &
                                 g_num(20,iel),g_num(14,iel),g_num(10,iel),   &
                                 g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1 "
  END DO

  CLOSE(11)

!------------------------------------------------------------------------------
! 12.3 Boundary conditions
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
  OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  CALL ns_cube_bc20(rest,nxe,nye,nze)
 
  DO i = 1, nr
    WRITE(12,'(I8,3I6)') rest(i,:) 
  END DO

  CLOSE(12)

!------------------------------------------------------------------------------
! 12.4 Loading conditions for the lid
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".lid" 
  OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  CALL loading_p126(no,nxe,nye,nze)
  val = 1.0
 
  DO i = 1, fixed_freedoms
    WRITE(13,'(I10,E12.4)') no(i), val(i) 
  END DO

  CLOSE(13)

!------------------------------------------------------------------------------
! 12.5 New control data
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
  OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  meshgen = 2 ! current default

  WRITE(14,'(A)') "'p126'"
  WRITE(14,'(I4)') meshgen
  WRITE(14,'(6I9)') nels, nn, nr, nip, fixed_freedoms
  WRITE(14,'(3E12.4,I8)') visc, rho, tol, limit
  WRITE(14,'(E12.4,I8,E12.4)') cjtol, cjits, penalty
  WRITE(14,'(E12.4,I8,E12.4)') x0, ell, kappa
  
  CLOSE(14)


!------------------------------------------------------------------------------
! 13. Program p127
!------------------------------------------------------------------------------

  CASE('p127')

  PRINT*
  PRINT*, "Program p127 not supported"
  PRINT*

  STOP

!------------------------------------------------------------------------------
! 14. Program p128
!------------------------------------------------------------------------------
 
  CASE('p128')
 
  PRINT*
  PRINT*, "Program p128 not supported"
  PRINT*
  PRINT*, "      Use program p128ar instead "
  PRINT*
  PRINT*, "      The <my_job>.mg file should contain the following data:"
  PRINT*
  PRINT*, "      p128ar"
  PRINT*, "      element"
  PRINT*, "      nels, nxe, nze, nod, nip"
  PRINT*, "      aa, bb, cc"
  PRINT*, "      rho, e, v"
  PRINT*, "      nev, ncv, bmat, which"
  PRINT*, "      tol, maxitr"
  PRINT*

  STOP

!------------------------------------------------------------------------------
! 15. Program p128ar
!------------------------------------------------------------------------------
 
  CASE('p128ar')
  READ(10,*) element
  READ(10,*) meshgen
  READ(10,*) nels, nxe, nze, nod, nip
  READ(10,*) aa, bb, cc
  READ(10,*) rho, e, v
  READ(10,*) nev, ncv, bmat, which
  READ(10,*) tol, maxitr

  nye   = nels/nxe/nze
  nr    = (nxe+1) * (nze+1)
  nn    = (nxe+1) * (nze+1) * (nye+1)
  ndim  = 3
  nodof = 3

!------------------------------------------------------------------------------
! 15.1 Allocate dynamic arrays
!------------------------------------------------------------------------------

  ALLOCATE(g_coord(ndim,nn))
  ALLOCATE(coord(nod,ndim))
  ALLOCATE(g_num(nod,nels))
  ALLOCATE(rest(nr,nodof+1))
  
  g_coord = 0.0_iwp
  coord   = 0.0_iwp
  g_num   = 0
  rest    = 0

  DO i=1,nr
    rest(i,1) = i
  END DO

  DO iel=1,nels
    CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
    g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)    
  END DO
 
  fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
  OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  WRITE(11,'(A)') "*THREE_DIMENSIONAL"
  WRITE(11,'(A)') "*NODES"

  DO i = 1, nn
    WRITE(11,'(I6,3F12.4)') i, g_coord(:,i)
  END DO

  WRITE(11,'(A)') "*ELEMENTS"
  DO iel = 1, nels
    WRITE(11,'(I6,A,8I10,A)') iel, " 3 8 1 ", g_num(1,iel),g_num(4,iel),     & 
                              g_num(8,iel),g_num(5,iel),g_num(2,iel),        &
                              g_num(3,iel),g_num(7,iel),g_num(6,iel), " 1 "
  END DO

  CLOSE(11)  

  fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
  OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
 
  DO i = 1, nr
    WRITE(12,'(4I6)') rest(i,:) 
  END DO

  CLOSE(12)
 
!------------------------------------------------------------------------------
! 15.6 Write new control data file
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
  OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  meshgen = 2 ! Default

  WRITE(14,'(A)')         program_name
  WRITE(14,'(A)')         "'hexahedron'"
  WRITE(14,'(I8)')         meshgen
  WRITE(14,'(6I9)')       nels, nn, nr, nod, nip
  WRITE(14,'(3E12.4,I8)') rho, e, v
  WRITE(14,'(2I8,5A)')    nev, ncv, " '", bmat, "' '", which,"'"
  WRITE(14,'(E12.4,I8)')  tol, maxitr
  CLOSE(14)

!------------------------------------------------------------------------------
! 15.7 Deallocate arrays
!------------------------------------------------------------------------------

  DEALLOCATE(coord)
  DEALLOCATE(g_coord)
  DEALLOCATE(g_num)
  DEALLOCATE(rest)

!------------------------------------------------------------------------------
! 16. Program p129
!------------------------------------------------------------------------------  
  CASE('p129')
  
  DO i=1,nr
    rest(i,1) = i
  END DO

  DO iel=1,nels
    CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
    g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)    
  END DO
  
  fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
  OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  WRITE(11,'(A)') "*THREE_DIMENSIONAL"
  WRITE(11,'(A)') "*NODES"

  DO i = 1, nn
    WRITE(11,'(I6,3F12.4)') i, g_coord(:,i)
  END DO

  WRITE(11,'(A)') "*ELEMENTS"
  DO iel = 1, nels
    WRITE(11,'(I6,A,20I10,A)') iel, " 3 20 1 ", g_num(3,iel),g_num(5,iel),    &
                                 g_num(7,iel),g_num(1,iel),g_num(15,iel),    &
                                 g_num(17,iel),g_num(19,iel),g_num(13,iel),    &
                                 g_num(4,iel),g_num(6,iel),g_num(8,iel),    &
                                 g_num(2,iel),g_num(16,iel),g_num(18,iel),     &
                                 g_num(20,iel),g_num(14,iel),g_num(10,iel),    &
                                 g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1 "
  END DO

  CLOSE(11)  

  fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
  OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
 
  DO i = 1, nr
    WRITE(12,'(4I6)') rest(i,:) 
  END DO

  CLOSE(12)

!------------------------------------------------------------------------------
! 16. Program p1210
!------------------------------------------------------------------------------

  CASE('p1210')

  PRINT*
  PRINT*, "Program p1210 not supported"
  PRINT*

  STOP

!------------------------------------------------------------------------------
! 17. Case default
!------------------------------------------------------------------------------

  CASE DEFAULT

  PRINT*
  PRINT*, "Mesh only generated for programs: "
  PRINT*, "  p121, p122, p123, p124, p125, p126"
  PRINT*, "  p127, p128, p129 and p1210"
  PRINT*

  END SELECT
  
END PROGRAM sg12mg
