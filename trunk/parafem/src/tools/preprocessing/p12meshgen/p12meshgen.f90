PROGRAM p12meshgen

!/****h* tools/preprocessing/p12meshgen
!*  NAME
!*    PROGRAM p12meshgen
!*  FUNCTION
!*    Creates simple meshes using files with the .mg extension.
!*  AUTHOR
!*    Lee Margetts
!*  COPYRIGHT
!*    (c) University of Manchester 2007-2014
!****
!*/

  USE precision
  USE geometry
  USE loading
  USE input
  
  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER                :: nr,nn,nels,nres,nle,nxe,nye,nze,nlen
  INTEGER                :: nod,ndim,nodof,nip 
  INTEGER                :: nmodes,lalfa,leig,lx,lz
  INTEGER                :: i,j,iel,l,m,n,iargc,argc
  INTEGER                :: plasits,cjits,limit,maxitr,incs
  INTEGER                :: ell
  INTEGER                :: fixed_nodes,fixed_freedoms,loaded_freedoms
  INTEGER                :: nev,ncv
  INTEGER                :: meshgen,partitioner
  INTEGER                :: nstep,npri,count
  INTEGER                :: np_types,nprops
  INTEGER                :: prnwidth,remainder
  REAL(iwp),PARAMETER    :: zero = 0.0_iwp, twelth = 1.0_iwp/12.0_iwp
  REAL(iwp)              :: aa,bb,cc
  REAL(iwp)              :: kx,ky,kz,cp
  REAL(iwp)              :: cjtol  
  REAL(iwp)              :: rho,visc,e,v   
  REAL(iwp)              :: el,er,acc  
  REAL(iwp)              :: kappa,penalty,tol,x0  
  REAL(iwp)              :: alpha1,beta1,theta,omega,dtim
  REAL(iwp)              :: zbox,zele,z1,z3,ttl
  REAL(iwp)              :: val0
  REAL(iwp)              :: sbary
  REAL(iwp)              :: phi,c,psi,cons,plastol
  CHARACTER(LEN=15)      :: element
  CHARACTER(LEN=50)      :: program_name,argv,problem_type
  CHARACTER(LEN=50)      :: iotype
  CHARACTER(LEN=1)       :: bmat
  CHARACTER(LEN=2)       :: which
  LOGICAL                :: solid=.true.
!-------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------- 

  INTEGER, ALLOCATABLE   :: g_num(:,:),rest(:,:),nf(:,:),no(:),no_f(:)
  INTEGER, ALLOCATABLE   :: num(:),etype(:)
  REAL(iwp), ALLOCATABLE :: g_coord(:,:),coord(:,:),val(:),val_f(:),qinc(:)
  REAL(iwp), ALLOCATABLE :: oldlds(:)

!-------------------------------------------------------------------------
! 3. Read data file base name "argv" from the command line
!-------------------------------------------------------------------------

  argc = iargc()
  IF (argc /= 1) THEN
    PRINT*
    PRINT*,"Usage: p12meshgen <base_name>"
    PRINT*
    PRINT*,"       program expects <base_name>.mg and outputs any or all of"
    PRINT*,"       <base_name>.dat <base_name>.d" 
    PRINT*,"       <base_name>.bnd <base_name>.lds"
    PRINT*
    STOP
  END IF
  CALL getname(argv,nlen)

!------------------------------------------------------------------------------
! 4. Read control data
!------------------------------------------------------------------------------

! IF (INDEX(argv,".mg") /= 0) THEN
!   argv = argv(1:INDEX(argv,".mg")-1)
! END IF
! fname = argv(1:INDEX(argv," ")-1) // ".mg"

  OPEN (10, FILE=argv(1:nlen)//'.mg', status='old', action='read')

  READ(10,*) program_name

  PRINT *
  PRINT *, "Generating input deck for program ", program_name
  PRINT *

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! 5. Select program using the SELECT CASE CONSTRUCT
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

  SELECT CASE(program_name)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Program p121
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

  CASE('p121')
  
    READ(10,*) iotype,nels,nxe,nze,nod,nip,aa,bb,cc,e,v,tol,limit

    nye = nels/nxe/nze

!-------------------------------------------------------------------------
! p121.1 Select 8 node or 20 node hexahedra
!-------------------------------------------------------------------------

    SELECT CASE(nod)
  
    CASE(20)
   
      nr    = ((2*nxe+1)*(nze+1)+(nxe+1)*nze)*2 +                        &
              ((2*nye-1)*nze+(nye-1)*nze)*2     +                        &
               (2*nye-1)*(nxe+1)+(nye-1)*nxe

      nn    = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1) +              &
               (nxe+1)*(nze+1)*nye

      ndim  = 3
      nodof = 3
      nle   = nxe/5
      
      fixed_freedoms  = 0
      loaded_freedoms = 3*nle*nle + 4*nle + 1
      element         = 'hexahedron'

      ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),num(nod),&
               rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms))
            
      coord = 0.0_iwp ; g_coord = 0.0_iwp ; val = 0.0_iwp
      g_num = 0       ; rest    = 0       ; no  = 0       ; num = 0
          
!------------------------------------------------------------------------------
! p121.2  Find nodal coordinates and element steering array
!         Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

      DO iel = 1, nels
        CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
        g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
      END DO
          
!------------------------------------------------------------------------------
! p121.3  Boundary conditions
!------------------------------------------------------------------------------
            
      CALL cube_bc20(rest,nxe,nye,nze)
          
!------------------------------------------------------------------------------
! p121.4  Loading conditions
!------------------------------------------------------------------------------
            
      CALL load_p121(nle,nod,nxe,nze, no,val)
      val = -val * aa * bb * (25._iwp / 12._iwp)
              
!------------------------------------------------------------------------------
! p121.5  Create input deck for 8 node hexahedra
!------------------------------------------------------------------------------
            
      CASE(8)
          
        nr    = ((nxe+1)*(nze+1))*2 + ((nye-1)*(nze+1))*2 + ((nxe-1)*(nze-1)) 
        ndim  = 3
        nodof = 3
        nn    = (nxe+1)*(nye+1)*(nze+1)
           
        nle             = nxe/5
        loaded_freedoms = (nle+1)*(nle+1)
        fixed_freedoms  = 0
        element         = 'hexahedron'

        ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),           &
                 rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),  &
                 num(nod))
            
        coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
        g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0

!------------------------------------------------------------------------------
! p121.6  Find nodal coordinates and element steering array
!         write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

        DO iel = 1, nels
          CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
          g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
        END DO

!------------------------------------------------------------------------------
! p121.7  Boundary conditions
!------------------------------------------------------------------------------
         
        CALL cube_bc8(rest,nxe,nye,nze)
          
!------------------------------------------------------------------------------
! p121.8  Loading conditions
!------------------------------------------------------------------------------
         
        CALL load_p121(nle,nod,nxe,nze, no,val)
        val = val * aa * bb 

!------------------------------------------------------------------------------
! p121.9  Default case and error message
!------------------------------------------------------------------------------
              
        CASE DEFAULT
        
          PRINT *
          PRINT *, "Wrong value given in variable NOD"
          PRINT *, "  Accepted values are 8 and 20"
          PRINT *, "  Here NOD = ", nod
          PRINT *
             
          iotype = 'aborted'
             
        END SELECT

!------------------------------------------------------------------------------
! p121.10  Output model
!------------------------------------------------------------------------------
         
        SELECT CASE(iotype)

        CASE('parafem')

!------------------------------------------------------------------------------
! p121.11  Output geometry file ".d" 
!------------------------------------------------------------------------------
         
          OPEN(11,FILE=argv(1:nlen)//'.d',STATUS='REPLACE',ACTION='WRITE')
            
          WRITE(11,'(A)') "*THREE_DIMENSIONAL"
          WRITE(11,'(A)') "*NODES"
          
          DO i = 1, nn
            WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
          END DO
          
          WRITE(11,'(A)') "*ELEMENTS"
            
          IF(nod==20) THEN
            DO iel = 1, nels
               WRITE(11,'(I12,A,20I12,A)')iel," 3 20 1 ",g_num(3,iel),        &
                  g_num(5,iel),g_num(7,iel),g_num(1,iel),g_num(15,iel),       &
                  g_num(17,iel),g_num(19,iel),g_num(13,iel),g_num(4,iel),     &
                  g_num(6,iel),g_num(8,iel),g_num(2,iel),g_num(16,iel),       &
                  g_num(18,iel),g_num(20,iel),g_num(14,iel),g_num(10,iel),    &
                  g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1"
            END DO
          ELSE IF(nod==8) THEN
            DO iel = 1, nels
              WRITE(11,'(I12,A,8I12,A)') iel, " 3 8 1 ", g_num(1,iel),        &
                  g_num(2,iel),g_num(3,iel),g_num(4,iel),g_num(5,iel),        &
                  g_num(6,iel),g_num(7,iel),g_num(8,iel)," 1"
            END DO
          ELSE
              PRINT *, "Wrong number of nodes in element"
          END IF
          
          CLOSE(11)

!------------------------------------------------------------------------------
! p121.11  Output ".bnd" file
!------------------------------------------------------------------------------

          OPEN(12,FILE=argv(1:nlen)//'.bnd',STATUS='REPLACE',ACTION='WRITE')
          
          DO i = 1, nr
            WRITE(12,'(I8,3I6)') rest(i,:) 
          END DO
          
          CLOSE(12)

!------------------------------------------------------------------------------
! p121.12  Output ".lds" file
!------------------------------------------------------------------------------

          OPEN(13,FILE=argv(1:nlen)//'.lds',STATUS='REPLACE',ACTION='WRITE')
          
          DO i = 1, loaded_freedoms
            WRITE(13,'(I12,2A,3E16.8)') no(i),"  0.00000000E+00  ",           &
                                        "0.00000000E+00",val(i) 
          END DO
            
          CLOSE(13)

!------------------------------------------------------------------------------
! p121.13  Output new ".dat" file
!------------------------------------------------------------------------------

          OPEN(14,FILE=argv(1:nlen)//'.dat',STATUS='REPLACE',ACTION='WRITE')
          
          WRITE(14,'(A)') "'hexahedron'"
          IF(nod==8) THEN
            WRITE(14,'(A)') "1"            ! SGM node numbering scheme
          ELSE
            WRITE(14,'(A)') "2"            ! Abaqus node numbering scheme
          END IF
          WRITE(14,'(A)') "1"              ! Internal mesh partitioning
          WRITE(14,'(3I12,2I5,2I9)')nels,nn,nr,nip,nod,loaded_freedoms
          WRITE(14,'(3E12.4,I8,A)') e,v,tol,limit," 1" 
          
          CLOSE(14)
             
!------------------------------------------------------------------------------
! p121.14  Output for visualization in ParaView (Ensight Gold)
!------------------------------------------------------------------------------

          CASE('paraview')

            ! modify mesh_ensi for optional arguments

            ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 
            etype=0; nf=0

            nstep=1; npri=1; dtim=1.0; solid=.true. 

            CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,          &
                           oldlds(1:),nstep,npri,dtim,solid)

          CASE DEFAULT

            PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

          END SELECT 

!------------------------------------------------------------------------------ 
!------------------------------------------------------------------------------
! Program p122
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p122')
  
    problem_type='ed4' ! hardwired
    READ(10,*) iotype, nels, nxe, nze, nod, nip
    READ(10,*) aa, bb, cc, incs
    READ(10,*) phi,c,psi,e,v
    READ(10,*) plasits,cjits,plastol,cjtol

    ALLOCATE(qinc(incs)); qinc=0.0_iwp

    READ(10,*) qinc(:)
    
    nye   = nels/nxe/nze

!------------------------------------------------------------------------------
! p122.1  Select 8 node or 20 node hexahedra
!------------------------------------------------------------------------------

    SELECT CASE(nod)
   
!------------------------------------------------------------------------------
! p122.2  Create input deck for 20 node hexahedra
!------------------------------------------------------------------------------
  
    CASE(20)
    
      nr    = ((2*nxe+1)*(nze+1)+(nxe+1)*nze)*2 +                           &
              ((2*nye-1)*nze+(nye-1)*nze)*2+(2*nye-1)*(nxe+1)+(nye-1)*nxe
      ndim  = 3
      nodof = 3
      nn    = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1)+(nxe+1)*(nze+1)*nye
      
      IF(problem_type == 'ed4') THEN
        nle             = nxe/5
        loaded_freedoms = 3*nle*nle + 4*nle + 1
        fixed_freedoms  = 0
      ELSE IF(problem_type == 'boussinesq') THEN
        loaded_freedoms = 1
        fixed_freedoms  = 0
      ELSE
        PRINT *, "Problem type: ",problem_type," not recognised."
      END IF
  
!------------------------------------------------------------------------------
! p122.3  Allocate dynamic arrays
!------------------------------------------------------------------------------
  
      ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),              &
               rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),     &
               num(nod))
    
      coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
      g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0
  
!------------------------------------------------------------------------------
! p122.4  Find nodal coordinates and element steering array
!------------------------------------------------------------------------------

     DO iel = 1, nels
       CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
       g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
     END DO
  
!------------------------------------------------------------------------------
! p122.5  Boundary conditions
!------------------------------------------------------------------------------
    
     CALL cube_bc20(rest,nxe,nye,nze)
  
!------------------------------------------------------------------------------
! p122.6  Loading conditions
!------------------------------------------------------------------------------

    CALL load_p121(nle,nod,nxe,nze, no,val)
    val = -val * aa * bb / 12._iwp
  
!------------------------------------------------------------------------------
! p122.7  Create input deck for 8 node hexahedra
!------------------------------------------------------------------------------
    
   CASE(8)
  
     nr    = ((nxe+1)*(nze+1))*2 + ((nye-1)*(nze+1))*2 + ((nxe-1)*(nze-1)) 
     ndim  = 3
     nodof = 3
     nn    = (nxe+1)*(nye+1)*(nze+1)
    
     IF(problem_type == 'ed4') THEN
       nle             = nxe/5
       loaded_freedoms = (nle+1)*(nle+1)
       fixed_freedoms  = 0
     ELSE IF(problem_type == 'boussinesq') THEN
       loaded_freedoms = 1
       fixed_freedoms  = 0
     ELSE
       PRINT *, "Problem type: ",problem_type," not recognised."
     END IF

!------------------------------------------------------------------------------
! p122.8  Allocate dynamic arrays
!------------------------------------------------------------------------------

     ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),               &
              rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),      &
              num(nod))
    
     coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
     g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0

!------------------------------------------------------------------------------
! p122.9  Find nodal coordinates and element steering array
!------------------------------------------------------------------------------

     DO iel = 1, nels
       CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
       g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
     END DO

!------------------------------------------------------------------------------
! p122.10  Boundary conditions
!------------------------------------------------------------------------------

     CALL cube_bc8(rest,nxe,nye,nze)
  
!------------------------------------------------------------------------------
! p122.11  Loading conditions
!------------------------------------------------------------------------------

     CALL load_p121(nle,nod,nxe,nze, no,val)
     val = val * aa * bb 

!------------------------------------------------------------------------------
! p122.12  Default case and error message
!------------------------------------------------------------------------------
  
     CASE DEFAULT
  
     PRINT *
     PRINT *, "Wrong value given in variable NOD"
     PRINT *, "  Accepted values are 8 and 20"
     PRINT *, "  Here NOD = ", nod
     PRINT *
      
   END SELECT

!------------------------------------------------------------------------------
! p122.13  Output model
!------------------------------------------------------------------------------

   SELECT CASE(iotype)

   CASE('parafem')
   
!------------------------------------------------------------------------------
! p122.14  Output geometry file ".d" 
!------------------------------------------------------------------------------
         
     OPEN(11,FILE=argv(1:nlen)//'.d',STATUS='REPLACE',ACTION='WRITE')
       
     WRITE(11,'(A)') "*THREE_DIMENSIONAL"
     WRITE(11,'(A)') "*NODES"
     
     DO i = 1, nn
       WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
     END DO
     
     WRITE(11,'(A)') "*ELEMENTS"
       
     IF(nod==20) THEN
       DO iel = 1, nels
          WRITE(11,'(I12,A,20I12,A)')iel," 3 20 1 ",g_num(3,iel),        &
             g_num(5,iel),g_num(7,iel),g_num(1,iel),g_num(15,iel),       &
             g_num(17,iel),g_num(19,iel),g_num(13,iel),g_num(4,iel),     &
             g_num(6,iel),g_num(8,iel),g_num(2,iel),g_num(16,iel),       &
             g_num(18,iel),g_num(20,iel),g_num(14,iel),g_num(10,iel),    &
             g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1"
       END DO
     ELSE IF(nod==8) THEN
       DO iel = 1, nels
         WRITE(11,'(I12,A,8I12,A)') iel, " 3 8 1 ", g_num(1,iel),        &
             g_num(2,iel),g_num(3,iel),g_num(4,iel),g_num(5,iel),        &
             g_num(6,iel),g_num(7,iel),g_num(8,iel)," 1"
       END DO
     ELSE
         PRINT *, "Wrong number of nodes in element"
     END IF
     
     CLOSE(11)

!------------------------------------------------------------------------------
! p122.15  Output ".bnd" file
!------------------------------------------------------------------------------

     OPEN(12,FILE=argv(1:nlen)//'.bnd',STATUS='REPLACE',ACTION='WRITE')
     
     DO i = 1, nr
       WRITE(12,'(I8,3I6)') rest(i,:) 
     END DO
     
     CLOSE(12)

!------------------------------------------------------------------------------
! p122.16  Output ".lds" file
!------------------------------------------------------------------------------

     OPEN(13,FILE=argv(1:nlen)//'.lds',STATUS='REPLACE',ACTION='WRITE')
     
     DO i = 1, loaded_freedoms
       WRITE(13,'(I12,2A,3E16.8)') no(i),"  0.00000000E+00  ",           &
                                   "0.00000000E+00",val(i) 
     END DO
       
     CLOSE(13)

!------------------------------------------------------------------------------
! p122.17  New control data
!------------------------------------------------------------------------------

     OPEN(14,FILE=argv(1:nlen)//'.dat',STATUS='REPLACE',ACTION='WRITE')
  
     WRITE(14,'(A)') "'hexahedron'"
     IF(nod==8) THEN
       WRITE(14,'(A)') "1"            ! Abaqus node numbering scheme
     ELSE
       WRITE(14,'(A)') "2"            ! Abaqus node numbering scheme
     END IF
     WRITE(14,'(A)') "1"              ! Internal mesh partitioning
     WRITE(14,'(5I9,A,I9)') nels, nn, nr, nip, nod, "  0  ",loaded_freedoms
     WRITE(14,'(6E12.4,I8,A)') phi, c, psi, e, v
     WRITE(14,'(3I6,2E12.4)') incs, plasits, cjits, plastol, cjtol
     DO i=1,incs
       WRITE(14,'(E12.4)') qinc(i)
     END DO
     
     CLOSE(14)
             
!------------------------------------------------------------------------------
! p122.18  Output for visualization in ParaView (Ensight Gold)
!------------------------------------------------------------------------------
  
   CASE('paraview')

    ! modify mesh_ensi for optional arguments

      ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 
      etype=0; nf=0

      nstep=incs; npri=1; dtim=1.0; solid=.true. 

      CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,                &
                     oldlds(1:),nstep,npri,dtim,solid)


   CASE DEFAULT

     PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

   END SELECT 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE DEFAULT

    PRINT*
    PRINT*, "Mesh only generated for programs: "
    PRINT*, "  p121, p122, p123, p124, p125, p126"
    PRINT*, "  p127, p128, p129 and p1210"
    PRINT*

  END SELECT
          
END PROGRAM p12meshgen
