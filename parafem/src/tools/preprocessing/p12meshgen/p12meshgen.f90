PROGRAM p12meshgen

!/****h* tools/preprocessing/p12meshgen
!*  NAME
!*    PROGRAM p12meshgen
!*  FUNCTION
!*    Creates simple meshes using files with the .mg extension.
!*  AUTHOR
!*    Lee Margetts
!*  COPYRIGHT
!*    (c) University of Manchester 2021
!****
!*/

  USE precision
  USE geometry
  USE loading
  USE input
  USE new_library
  
  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER                :: nr,nn,nels,nres,nle,nxe,nye,nze,nlen
  INTEGER                :: nod,ndim,nodof,nip 
  INTEGER                :: nmodes,lalfa,leig,lx,lz
  INTEGER                :: i,j,k,iel,l,m,n,iargc,argc
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
  LOGICAL                :: found=.false.

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  INTEGER, ALLOCATABLE   :: g_num(:,:),rest(:,:),nf(:,:),no(:),no_f(:)
  INTEGER, ALLOCATABLE   :: num(:),etype(:)
  REAL(iwp), ALLOCATABLE :: g_coord(:,:),coord(:,:),val(:),val_f(:),qinc(:)
  REAL(iwp), ALLOCATABLE :: oldlds(:)

!------------------------------------------------------------------------------
! 3. Read data file base name "argv" from the command line
!------------------------------------------------------------------------------

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
  CALL getname_mg(argv,nlen)

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

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 5. Select program using the SELECT CASE CONSTRUCT
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SELECT CASE(program_name)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program p121
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p121')
  
    READ(10,*) iotype,nels,nxe,nze,nod,nip,aa,bb,cc,e,v,tol,limit

    nye = nels/nxe/nze

!------------------------------------------------------------------------------
! p121.1 Select 8 node or 20 node hexahedra
!------------------------------------------------------------------------------

    SELECT CASE(nod)
  
    CASE(20)
   
      nr    = ((2*nxe+1)*(nze+1)+(nxe+1)*nze)*2 +                             &
              ((2*nye-1)*nze+(nye-1)*nze)*2     +                             &
               (2*nye-1)*(nxe+1)+(nye-1)*nxe

      nn    = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1) +                   &
               (nxe+1)*(nze+1)*nye

      ndim  = 3
      nodof = 3
      nle   = nxe/5
      
      fixed_freedoms  = 0
      loaded_freedoms = 3*nle*nle + 4*nle + 1
      element         = 'hexahedron'

      ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),num(nod),     &
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

        ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),            &
                 rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),   &
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
            WRITE(12,'(I15,3I6)') rest(i,:) 
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
          WRITE(14,'(3E12.4,I8)') e,v,tol,limit 
          
          CLOSE(14)
             
!------------------------------------------------------------------------------
! p121.14  Output for visualization in ParaView (Ensight Gold)
!------------------------------------------------------------------------------

          CASE('paraview')

            ! modify mesh_ensi for optional arguments

            ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 

            etype  = 1   ! Only one material type in this mesh
            nf     = 0
            oldlds = zero

            nstep=1; npri=1; dtim=1.0; solid=.true. 

            k=0 
            DO j=1,loaded_freedoms
              k=k+1
              found=.false.
              DO i=1,nn
                IF(i==no(k)) THEN
                  l=i*3
                  oldlds(l)=val(k)
                  found=.true.
                END IF
                IF(found)CYCLE
              END DO
            END DO

            CALL rest_to_nf(rest,nf)

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

    element = 'hexahedron'

!------------------------------------------------------------------------------
! p122.1  Select 8 node or 20 node hexahedra
!------------------------------------------------------------------------------

    SELECT CASE(nod)
   
!------------------------------------------------------------------------------
! p122.2  Create input deck for 20 node hexahedra
!------------------------------------------------------------------------------
  
    CASE(20)
    
      nr    = ((2*nxe+1)*(nze+1)+(nxe+1)*nze)*2 +                             &
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

      ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),              &
               rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),     &
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
           WRITE(11,'(I12,A,20I12,A)')iel," 3 20 1 ",g_num(3,iel),            &
              g_num(5,iel),g_num(7,iel),g_num(1,iel),g_num(15,iel),           &
              g_num(17,iel),g_num(19,iel),g_num(13,iel),g_num(4,iel),         &
              g_num(6,iel),g_num(8,iel),g_num(2,iel),g_num(16,iel),           &
              g_num(18,iel),g_num(20,iel),g_num(14,iel),g_num(10,iel),        &
              g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1"
        END DO
      ELSE IF(nod==8) THEN
        DO iel = 1, nels
          WRITE(11,'(I12,A,8I12,A)') iel, " 3 8 1 ", g_num(1,iel),            &
              g_num(2,iel),g_num(3,iel),g_num(4,iel),g_num(5,iel),            &
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
        WRITE(13,'(I12,2A,3E16.8)') no(i),"  0.00000000E+00  ",               &
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

      etype  = 1
      nf     = 0
      oldlds = zero

      nstep=incs; npri=1; dtim=1.0; solid=.true. 
    
      k=1 
      DO j=1,loaded_freedoms
        k=k+1
        found=.false.
        DO i=1,nn
          IF(i==no(k)) THEN
            l=i*3
            oldlds(l)=val(k)
            found=.true.
          END IF
          IF(found) EXIT
        END DO
      END DO
        
      CALL rest_to_nf(rest,nf)

      CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,                &
                     oldlds(1:),nstep,npri,dtim,solid)


    CASE DEFAULT

      PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

    END SELECT 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program p123
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p123')

    READ(10,*) iotype, nels, nxe, nze, nip
    READ(10,*) aa, bb, cc, kx, ky, kz
    READ(10,*) tol, limit
    READ(10,*) loaded_freedoms, fixed_freedoms
 
    PRINT *, "Read .mg file"

!------------------------------------------------------------------------------
! p123.1 Initialize variables
!------------------------------------------------------------------------------

    nye   = nels/nxe/nze
    ndim  = 3
    nod   = 8
    nr    = (nxe+1)*(nye+1) + (nxe+1)*nze + nye*nze
    nn    = (nxe+1)*(nye+1)*(nze+1)
    nodof = 1
    nres  = nxe*(nze-1)+1

    element='hexahedron'

!------------------------------------------------------------------------------
! p123.2 Allocate dynamic arrays
!------------------------------------------------------------------------------
  
    ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),                &
             rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),       &
             num(nod),val_f(fixed_freedoms),no_f(fixed_freedoms))
    
    coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
    g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0
    val_f    = 0.0_iwp ; no_f    = 0
  
!------------------------------------------------------------------------------
! p123.3 Find nodal coordinates and element steering array
!        Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

    DO iel = 1, nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
    END DO
    
    SELECT CASE(iotype)

    CASE('parafem')

    OPEN(11,FILE=argv(1:nlen)//".d",STATUS='REPLACE',ACTION='WRITE')
    
    WRITE(11,'(A)') "*THREE_DIMENSIONAL"
    WRITE(11,'(A)') "*NODES"
  
    DO i = 1, nn
      WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
    END DO
 
    DEALLOCATE(g_coord)
 
    WRITE(11,'(A)') "*ELEMENTS"
    
    DO iel = 1, nels
      WRITE(11,'(I12,A,8I10,A)') iel, " 3 8 1 ", g_num(1,iel),g_num(4,iel),  &
                                   g_num(8,iel),g_num(5,iel),g_num(2,iel),   &
                                   g_num(3,iel),g_num(7,iel),g_num(6,iel),   &
                                    " 1"
    END DO
    
    CLOSE(11)

    PRINT *, "Output nodal coordinates and element steering array"

!------------------------------------------------------------------------------
! p123.4 Boundary conditions
!------------------------------------------------------------------------------
  
    OPEN(12,FILE=argv(1:nlen)//".bnd",STATUS='REPLACE',ACTION='WRITE')
  
    CALL box_bc8(rest,nxe,nye,nze)
  
    DO i = 1, nr
      WRITE(12,'(I8,3I6)') rest(i,:) 
    END DO
  
    CLOSE(12)

    PRINT *, "Output boundary conditions"

!------------------------------------------------------------------------------
! p123.5 Loading conditions
!------------------------------------------------------------------------------

    IF(loaded_freedoms > 0) THEN
     
      OPEN(13,FILE=argv(1:nlen)//".lds",STATUS='REPLACE',ACTION='WRITE')
     
      no   = nres
      val  = 10.0_iwp
  
      DO i = 1, loaded_freedoms
        WRITE(13,'(I10,E16.8)') no(i),val(i)
      END DO

      CLOSE(13)

      PRINT *, "Output fixed loads"

    END IF

    IF(fixed_freedoms>0) THEN

      OPEN(14,FILE=argv(1:nlen)//".fix",STATUS='REPLACE',ACTION='WRITE')

      no_f  = nres
      val_f = 100.0_iwp

      DO i = 1, fixed_freedoms
        WRITE(14,'(I10,E16.8)') no_f(i),val_f(i)
      END DO

      CLOSE(14)

      PRINT *, "Output fixed freedoms"

    END IF

!------------------------------------------------------------------------------
! p123.6 New control data
!------------------------------------------------------------------------------

      OPEN(15,FILE=argv(1:nlen)//".dat",STATUS='REPLACE',ACTION='WRITE')
  
      WRITE(15,'(A)') "'hexahedron'"
      WRITE(15,'(A)') "2"              ! Abaqus node numbering scheme
      WRITE(15,'(A)') "1"              ! Internal mesh partitioning
      WRITE(15,'(3I12,4I6)') nels, nn, nr, nip, nod, loaded_freedoms, fixed_freedoms
      WRITE(15,'(4E12.4,2I8)') kx, ky, kz, tol, limit, nres

      CLOSE(15)

      PRINT *, "Output new control data file"
      PRINT *, "Job completed"
      PRINT *

    CASE('paraview')

      ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 
      etype=0; nf=0

      nstep=1; npri=1; dtim=1.0; solid=.true. 

      CALL rest_to_nf(rest,nf)

      CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,                &
                     oldlds(1:),nstep,npri,dtim,solid)

    CASE DEFAULT

      PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

    END SELECT

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program p124
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p124')

    READ(10,*) iotype, nels, nxe, nze, nip
    READ(10,*) aa, bb, cc, kx, ky, kz, rho, cp
    READ(10,*) dtim, nstep, theta
    READ(10,*) npri, tol, limit, val0
    READ(10,*) np_types, loaded_freedoms, fixed_freedoms

    element='hexahedron'
  
    IF(np_types>1) PRINT *, "NP_TYPES will be set to 1"
    PRINT *, "Read .mg file"

!------------------------------------------------------------------------------
! p124.1 Initialize variables
!------------------------------------------------------------------------------

    nye    = nels/nxe/nze
    ndim   = 3
    nod    = 8
    nr     = (nxe+1)*(nye+1) + (nxe+1)*nze + nye*nze
    nn     = (nxe+1)*(nye+1)*(nze+1)
    nodof  = 1
    nres   = (nxe)*(nze-1)+1
!   dtim   = 0.01_iwp
!   nstep  = 150
    theta  = 0.5
!   npri   = 10
    nprops = 5

!------------------------------------------------------------------------------
! p124.2 Allocate dynamic arrays
!------------------------------------------------------------------------------
  
    ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),              &
             rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),     &
             num(nod),val_f(fixed_freedoms),no_f(fixed_freedoms))
    
    coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
    g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0
    val_f    = 0.0_iwp ; no_f    = 0
  
!------------------------------------------------------------------------------
! p124.3 Find nodal coordinates and element steering array
!        Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

    DO iel = 1, nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
    END DO

    SELECT CASE(iotype)

    CASE('parafem')
    
      OPEN(11,FILE=argv(1:nlen)//".d",STATUS='REPLACE',ACTION='WRITE')
    
      WRITE(11,'(A)') "*THREE_DIMENSIONAL"
      WRITE(11,'(A)') "*NODES"
  
      DO i = 1, nn
        WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
      END DO
  
      WRITE(11,'(A)') "*ELEMENTS"
    
      DO iel = 1, nels
        WRITE(11,'(I12,A,8I12,A)') iel, " 3 8 1 ", g_num(1,iel),g_num(4,iel),  &
                                     g_num(8,iel),g_num(5,iel),g_num(2,iel),   &
                                     g_num(3,iel),g_num(7,iel),g_num(6,iel),   &
                                      " 1"
      END DO
    
      CLOSE(11)

      PRINT *, "Output nodal coordinates and element steering array"

!------------------------------------------------------------------------------
! p124.4 Boundary conditions
!------------------------------------------------------------------------------
  
      OPEN(12,FILE=argv(1:nlen)//".bnd",STATUS='REPLACE',ACTION='WRITE')
  
      CALL box_bc8(rest,nxe,nye,nze)
  
      DO i = 1, nr
        WRITE(12,'(I12,3I6)') rest(i,:) 
      END DO
  
      CLOSE(12)

      PRINT *, "Output boundary conditions"

!------------------------------------------------------------------------------
! p124.5 Loading conditions
!------------------------------------------------------------------------------

      IF(loaded_freedoms > 0) THEN
     
        OPEN(13,FILE=argv(1:nlen)//".lds",STATUS='REPLACE',ACTION='WRITE')
     
        no   = nres
        val  = 10.0_iwp
  
        DO i = 1, loaded_freedoms
          WRITE(13,'(I12,E16.8)') no(i),val(i)
        END DO

        CLOSE(13)

        PRINT *, "Output fixed loads"
 
      END IF

      IF(fixed_freedoms>0) THEN

        OPEN(14,FILE=argv(1:nlen)//".fix",STATUS='REPLACE',ACTION='WRITE')

        no_f  = nres
        val_f = 100.0_iwp
 
        DO i = 1, fixed_freedoms
          WRITE(14,'(I12,E16.8)') no_f(i),val_f(i)
        END DO

        CLOSE(14)

        PRINT *, "Output fixed freedoms"

      END IF

!------------------------------------------------------------------------------
! p124.6 Materials file
!------------------------------------------------------------------------------

      OPEN(15,FILE=argv(1:nlen)//".mat",STATUS='REPLACE',ACTION='WRITE')

      WRITE(15,'(A,2I5)')    "*MATERIAL", np_types, nprops
      WRITE(15,'(A)')        "<edit material_name>"
      WRITE(15,'(A,5E12.4)') " 1", kx, ky, kz, rho, cp
     
      CLOSE(15)
     
!------------------------------------------------------------------------------
! p124.7 New control data
!------------------------------------------------------------------------------

      OPEN(16,FILE=argv(1:nlen)//".dat",STATUS='REPLACE',ACTION='WRITE')
  
      WRITE(16,'(A)') "'hexahedron'"
      WRITE(16,'(A)') "2"              ! Abaqus node numbering scheme
      WRITE(16,'(A)') "1"              ! Internal mesh partitioning
      WRITE(16,'(A)') "1"              ! Number of materials
      WRITE(16,'(2I12,5I9)') nels, nn, nr, nip, nod, loaded_freedoms, fixed_freedoms
      WRITE(16,'(E12.4)') val0
      WRITE(16,'(E12.4,2I8,E12.4)') dtim, nstep, npri, theta 
      WRITE(16,'(E12.4,2I8)') tol, limit, nres

      CLOSE(15)

      PRINT *, "Output new control data file"
      PRINT *, "Some values have default values"
      PRINT *, "Job completed"
      PRINT *

    CASE('paraview')

      ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 
      etype=0; nf=0

      solid=.true. 

      CALL rest_to_nf(rest,nf)

      CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,                &
                     oldlds(1:),nstep,npri,dtim,solid)

    CASE DEFAULT

      PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

    END SELECT

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program p125
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p125')

    READ(10,*) iotype, nels, nxe, nze, nip
    READ(10,*) aa, bb, cc, kx, ky, kz
    READ(10,*) dtim, nstep
    READ(10,*) npri, val0
    
    loaded_freedoms=0; fixed_freedoms=0

    element='hexahedron'

    PRINT *, "Read .mg file"

!------------------------------------------------------------------------------
! p125.1 Initialize variables
!------------------------------------------------------------------------------

    nye   = nels/nxe/nze
    ndim  = 3
    nod   = 8
    nr    = (nxe+1)*(nye+1) + (nxe+1)*nze + nye*nze
    nn    = (nxe+1)*(nye+1)*(nze+1)
    nodof = 1
    nres  = nxe*(nze-1)+1

!------------------------------------------------------------------------------
! p125.2 Allocate dynamic arrays
!------------------------------------------------------------------------------
  
    ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels))
    
    coord    = 0.0_iwp ; g_coord = 0.0_iwp ; g_num = 0 
  
!------------------------------------------------------------------------------
! p125.3 Find nodal coordinates and element steering array
!------------------------------------------------------------------------------

    DO iel = 1, nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
    END DO

!------------------------------------------------------------------------------
! p125.4 Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

    SELECT CASE(iotype)

    CASE('parafem')

    OPEN(11,FILE=argv(1:nlen)//".d",STATUS='REPLACE',ACTION='WRITE')
    
    WRITE(11,'(A)') "*THREE_DIMENSIONAL"
    WRITE(11,'(A)') "*NODES"
  
    DO i = 1, nn
      WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
    END DO

    DEALLOCATE(g_coord)
  
    WRITE(11,'(A)') "*ELEMENTS"
    
    DO iel = 1, nels
      WRITE(11,'(I12,A,8I12,A)') iel, " 3 8 1 ", g_num(1,iel),g_num(4,iel),   &
                                   g_num(8,iel),g_num(5,iel),g_num(2,iel),    &
                                   g_num(3,iel),g_num(7,iel),g_num(6,iel),    &
                                    " 1"
    END DO

    CLOSE(11)

    PRINT *, "Output nodal coordinates and element steering array"

!------------------------------------------------------------------------------
! p125.5 Boundary conditions
!------------------------------------------------------------------------------
  
    ALLOCATE(rest(nr,nodof+1))
    rest  = 0

    OPEN(12,FILE=argv(1:nlen)//".bnd",STATUS='REPLACE',ACTION='WRITE')
  
    CALL box_bc8(rest,nxe,nye,nze)
  
    DO i = 1, nr
      WRITE(12,'(I12,3I6)') rest(i,:) 
    END DO
  
    CLOSE(12)

    PRINT *, "Output boundary conditions"

!------------------------------------------------------------------------------
! p125.6 Loading conditions
!------------------------------------------------------------------------------

    IF(loaded_freedoms > 0) THEN
     
      OPEN(13,FILE=argv(1:nlen)//".lds",STATUS='REPLACE',ACTION='WRITE')
     
      no   = nres
      val  = 10.0_iwp
  
      DO i = 1, loaded_freedoms
        WRITE(13,'(I11,E16.8)') no(i),val(i)
      END DO

      CLOSE(13)

      PRINT *, "Output fixed loads"

    END IF

    IF(fixed_freedoms>0) THEN

      OPEN(14,FILE=argv(1:nlen)//".fix",STATUS='REPLACE',ACTION='WRITE')

      no_f  = nres
      val_f = 100.0_iwp

      DO i = 1, fixed_freedoms
        WRITE(14,'(I12,E16.8)') no_f(i),val_f(i)
      END DO

      CLOSE(14)

      PRINT *, "Output fixed freedoms"

    END IF

!------------------------------------------------------------------------------
! p125.7 New control data
!------------------------------------------------------------------------------

    OPEN(15,FILE=argv(1:nlen)//".dat",STATUS='REPLACE',ACTION='WRITE')
  
    WRITE(15,'(A)') "'hexahedron'"
    WRITE(15,'(A)') "2"              ! Abaqus node numbering scheme
    WRITE(15,'(A)') "1"              ! Internal mesh partitioning
    WRITE(15,'(7I9)') nels, nn, nr, nip, nod, loaded_freedoms, fixed_freedoms
    WRITE(15,'(4E12.4,I8)') kx, ky, kz
    WRITE(15,'(E12.4,I8,I8,E12.4)') dtim, nstep, npri
    WRITE(15,'(I8,E12.4)') nres, val0
     
    CLOSE(15)
     
    PRINT *, "Output new control data file"
    PRINT *, "Some values have default values"
    PRINT *, "Job completed"
    PRINT *

    CASE('paraview')

      ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim),num(nod)) 
      etype=0; nf=0; num=0

      solid=.true. 

      CALL rest_to_nf(rest,nf)

      CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,                &
                     oldlds(1:),nstep,npri,dtim,solid)

    CASE DEFAULT

      PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

    END SELECT 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program p126
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p126')

    READ(10,*) iotype, nels, nxe, nze, nip
    READ(10,*) aa, bb, cc
    READ(10,*) visc, rho, tol, limit
    READ(10,*) cjtol, cjits, penalty
    READ(10,*) x0, ell, kappa
  
    nye            = nels/nxe/nze
    fixed_freedoms = 3*nxe*nye+2*nxe+2*nye+1
    nr             = 3*nxe*nye*nze+4*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+2
    ndim           = 3
    nodof          = 4
    nod            = 20
    nn             = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1) +            &
                     (nxe+1)*(nze+1)*nye

    element="hexahedron"

!------------------------------------------------------------------------------
! p126.1 Allocate dynamic arrays
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
! p126.2 Find nodal coordinates and element steering array
!------------------------------------------------------------------------------

    DO iel = 1, nels
      CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
    END DO

!------------------------------------------------------------------------------
! p126.3 Find boundary conditions and apply velocities to lid
!------------------------------------------------------------------------------

    CALL ns_cube_bc20(rest,nxe,nye,nze)
    CALL loading_p126(no,nxe,nye,nze)
    val = 1.0
 
!------------------------------------------------------------------------------
! p126.4 Select data format
!------------------------------------------------------------------------------

    SELECT CASE(iotype)

!------------------------------------------------------------------------------
! p126.5 Output in proprietary ParaFEM format (based on Kidger's DanFE)
!------------------------------------------------------------------------------

    CASE('parafem')
 
      PRINT *, "  Writing nodal coordinates and element steering array"
      OPEN(11,FILE=argv(1:nlen)//".d",STATUS='REPLACE',ACTION='WRITE')
  
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
! p126.6 Write boundary conditions in ParaFEM format
!------------------------------------------------------------------------------

      PRINT *, "  Writing boundary conditions"

      OPEN(12,FILE=argv(1:nlen)//".bnd",STATUS='REPLACE',ACTION='WRITE')

      DO i = 1, nr
        WRITE(12,'(I8,4I6)') rest(i,:) 
      END DO

      CLOSE(12)

!------------------------------------------------------------------------------
! p126.7 Write loading conditions for the lid in ParaFEM format
!------------------------------------------------------------------------------

      PRINT *, "  Writing loading conditions for the lid"

      OPEN(13,FILE=argv(1:nlen)//".lid",STATUS='REPLACE',ACTION='WRITE')

       DO i = 1, fixed_freedoms
        WRITE(13,'(I10,E12.4)') no(i), val(i) 
      END DO

      CLOSE(13)

!------------------------------------------------------------------------------
! p126.8 Write new control data file in ParaFEM format
!------------------------------------------------------------------------------

      PRINT *, "  Writing new control data file"
 
      OPEN(14,FILE=argv(1:nlen)//".dat",STATUS='REPLACE',ACTION='WRITE')

      meshgen     = 2 ! current default
      partitioner = 1 ! current default
      nres        = nze*(nxe+1)+3*nxe+1

      WRITE(14,'(A)') "'p126'"
      WRITE(14,'(2I4)') meshgen, partitioner
      WRITE(14,'(7I9)') nels, nn, nres, nr, fixed_freedoms, nip
      WRITE(14,'(3E12.4,I8)') visc, rho, tol, limit
      WRITE(14,'(E12.4,I8,E12.4)') cjtol, cjits, penalty
      WRITE(14,'(E12.4,I8,E12.4)') x0, ell, kappa
  
      CLOSE(14)

!------------------------------------------------------------------------------
! 12.9 Write a .dis file so that boundary conditions can be checked
!------------------------------------------------------------------------------

      PRINT *, "  Writing dummy .dis file with flags for boundary conditions"

      OPEN(15,FILE=argv(1:nlen)//".dis",STATUS='REPLACE',ACTION='WRITE')

      WRITE(15,'(A)') "*DISPLACEMENT"
      WRITE(15,'(A)') " 1"
 
      j = 1
 
      DO i = 1,nn
        IF(rest(j,1) == i) THEN
          j = j+1
          IF(rest(j,2) == 0) THEN
            WRITE(15,'(I8,3A)') i, " 1.0000E+00", "  0.0000E+00", "  0.0000E+00"
          ELSE
            WRITE(15,'(I8,3A)') i, " 0.0000E+00", "  0.0000E+00", "  0.0000E+00"
          END IF
        ELSE
          WRITE(15,'(I8,3A)') i, " 0.0000E+00", "  0.0000E+00", "  0.0000E+00"
        END IF
      END DO

      WRITE(15,'(A)') "*DISPLACEMENT"
      WRITE(15,'(A)') " 2"
  
      j = 1

      DO i = 1,nn
        IF(rest(j,1) == i) THEN
          j = j + 1
          IF(rest(j,3) == 0) THEN
            WRITE(15,'(I8,3A)') i, " 1.0000E+00", "  0.0000E+00", "  0.0000E+00"
          ELSE
            WRITE(15,'(I8,3A)') i, " 0.0000E+00", "  0.0000E+00", "  0.0000E+00"
          END IF
        ELSE
          WRITE(15,'(I8,3A)') i, " 0.0000E+00", "  0.0000E+00", "  0.0000E+00"
        END IF
      END DO

      WRITE(15,'(A)') "*DISPLACEMENT"
      WRITE(15,'(A)') " 3"
  
      j = 1

      DO i = 1,nn
        IF(rest(j,1) == i) THEN
          j = j + 1
          IF(rest(j,4) == 0) THEN
            WRITE(15,'(I8,3A)') i, " 1.0000E+00", "  0.0000E+00", "  0.0000E+00"
          ELSE
            WRITE(15,'(I8,3A)') i, " 0.0000E+00", "  0.0000E+00", "  0.0000E+00"
          END IF
        ELSE
          WRITE(15,'(I8,3A)') i, " 0.0000E+00", "  0.0000E+00", "  0.0000E+00"
        END IF
      END DO

      WRITE(15,'(A)') "*DISPLACEMENT"
      WRITE(15,'(A)') " 4"
 
      j = 1
 
      DO i = 1,nn
        IF(rest(j,1) == i) THEN
          j = j + 1
          IF(rest(j,5) == 0) THEN
            WRITE(15,'(I8,3A)') i, " 1.0000E+00", "  0.0000E+00", "  0.0000E+00"
          ELSE
            WRITE(15,'(I8,3A)') i, " 0.0000E+00", "  0.0000E+00", "  0.0000E+00"
          END IF
        ELSE
          WRITE(15,'(I8,3A)') i, " 0.0000E+00", "  0.0000E+00", "  0.0000E+00"
        END IF
      END DO

      CLOSE(15)

!------------------------------------------------------------------------------
! p126.10 Output in ASCII Ensight Gold format
!------------------------------------------------------------------------------

    CASE('paraview')
    
      ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 
      etype=0; nf=0

      nstep=1; npri= 1; dtim=1.0_iwp
      solid=.true. ! need to review this 

      CALL rest_to_nf(rest,nf)

      CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,                &
                     oldlds(1:),nstep,npri,dtim,solid)

    CASE DEFAULT

      PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

    END SELECT 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 13. Program p127
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p127')

    PRINT*
    PRINT*, "Program p127 is not supported"
    PRINT*, "Consult 'Programming the Finite Element Method', 5th Edition, &
             &for details"

    PRINT*

!   STOP
    
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program p128
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p128')
 
    READ(10,*) iotype,nels,nxe,nze,nip
    READ(10,*) aa,bb,cc
    READ(10,*) rho,e,v
    READ(10,*) nmodes
    READ(10,*) el,er
    READ(10,*) lalfa,leig,lx,lz,acc

    nye   = nels/nxe/nze
    nr    = (nxe+1) * (nze+1)
    nn    = (nxe+1) * (nze+1) * (nye+1)
    ndim  = 3
    nodof = 3
    nod   = 8

    element="hexahedron"
    
!------------------------------------------------------------------------------
! p128.1 Allocate dynamic arrays
!------------------------------------------------------------------------------

    ALLOCATE(g_coord(ndim,nn))
    ALLOCATE(coord(nod,ndim))
    ALLOCATE(g_num(nod,nels))
    ALLOCATE(rest(nr,nodof+1))
  
    g_coord = 0.0_iwp
    coord   = 0.0_iwp
    g_num   = 0
    rest    = 0

!------------------------------------------------------------------------------
! p128.2 Find nodal coordinates and element steering array
!------------------------------------------------------------------------------
      
    DO iel=1,nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)    
    END DO

!------------------------------------------------------------------------------
! p128.3 Find boundary conditions
!------------------------------------------------------------------------------

    DO i=1,nr
      rest(i,1) = i
    END DO

!------------------------------------------------------------------------------
! p128.4 Select data format
!------------------------------------------------------------------------------

    SELECT CASE(iotype)

!------------------------------------------------------------------------------
! p128.5 Output in proprietary ParaFEM format (based on Kidger's DanFE)
!------------------------------------------------------------------------------

      CASE('parafem')
 
        OPEN(11,FILE=argv(1:nlen)//".d",STATUS='REPLACE',ACTION='WRITE')
  
        WRITE(11,'(A)') "*THREE_DIMENSIONAL"
        WRITE(11,'(A)') "*NODES"

        DO i = 1, nn
          WRITE(11,'(I10,3F12.4)') i, g_coord(:,i)
        END DO

        WRITE(11,'(A)') "*ELEMENTS"
        DO iel = 1, nels
          WRITE(11,'(I10,A,8I10,A)') iel," 3 8 1 ",g_num(1,iel),g_num(4,iel), & 
                                  g_num(8,iel),g_num(5,iel),g_num(2,iel),     &
                                  g_num(3,iel),g_num(7,iel),g_num(6,iel), " 1 "
        END DO

        CLOSE(11)  

        OPEN(12,FILE=argv(1:nlen)//".bnd",STATUS='REPLACE',ACTION='WRITE')
 
        DO i = 1, nr
          WRITE(12,'(I10,3I3)') rest(i,:) 
        END DO

        CLOSE(12)
 
!------------------------------------------------------------------------------
! p128.6 Write new control data file in ParaFEM format
!------------------------------------------------------------------------------

        OPEN(14,FILE=argv(1:nlen)//".dat",STATUS='REPLACE',ACTION='WRITE')

        meshgen     = 2 ! default
        partitioner = 1 ! default

        WRITE(14,'(I2)')        meshgen
        WRITE(14,'(I2)')        partitioner
        WRITE(14,'(4I10)')      nels, nn, nr, nip
        WRITE(14,'(3E12.4)')    rho, e, v
        WRITE(14,'(I8)')        nmodes
        WRITE(14,'(2E12.4)')    el,er
        WRITE(14,'(4I6,E12.4)') lalfa,leig,lx,lz,acc

        CLOSE(14)

!------------------------------------------------------------------------------
! p128.7 Deallocate arrays
!------------------------------------------------------------------------------

        DEALLOCATE(coord)
        DEALLOCATE(g_coord)
        DEALLOCATE(g_num)
        DEALLOCATE(rest)

!------------------------------------------------------------------------------
! p128.8 Output in ASCII Ensight Gold format
!------------------------------------------------------------------------------

      CASE('paraview')

        ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 
        etype=0; nf=0

        nstep=nmodes; npri= 1; dtim=1.0_iwp
        solid=.true. ! need to review this 

        CALL rest_to_nf(rest,nf)

        CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,              &
                       oldlds(1:),nstep,npri,dtim,solid)

      CASE DEFAULT

        PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

      END SELECT     

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program p129
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p129')
 
    READ(10,*) iotype,nels,nxe,nze,nip
    READ(10,*) aa,bb,cc,rho
    READ(10,*) e,v
    READ(10,*) alpha1,beta1
    READ(10,*) nstep,npri,theta
    READ(10,*) omega,tol,limit

    nye   = nels/nxe/nze
    nr    = 3*nxe*nze+2*nxe+2*nze+1
    ndim  = 3
    nodof = 3
    nod   = 20
    nn    = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1)+(nxe+1)*(nze+1)*nye

    element = 'hexahedron' 

    loaded_freedoms = (2*nxe)+1 ! Differs from the problem in the book

    nres = 3*(nye*(nxe+1)*(nze+1)+nr*(nye-1)+(nxe+1))

!------------------------------------------------------------------------------
! p129.1 Allocate dynamic arrays
!------------------------------------------------------------------------------

    ALLOCATE(coord(nod,ndim))
    ALLOCATE(g_coord(ndim,nn))
    ALLOCATE(g_num(nod,nels))
    ALLOCATE(rest(nr,nodof+1))
    ALLOCATE(no(loaded_freedoms))
    ALLOCATE(val(loaded_freedoms))
 
    coord   = 0.0_iwp
    g_coord = 0.0_iwp
    g_num   = 0
    rest    = 0
    val     = 0.0_iwp
    no      = 0
 
!------------------------------------------------------------------------------
! p129.2 Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention
!------------------------------------------------------------------------------
 
    DO iel=1,nels
      CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)    
    END DO
  
!------------------------------------------------------------------------------
! p129.3 Boundary conditions
!------------------------------------------------------------------------------

    DO i=1,nr
      rest(i,1) = i
    END DO

!------------------------------------------------------------------------------
! p129.4 Loading conditions
!------------------------------------------------------------------------------

    count = 1

    DO i = 1,loaded_freedoms
      no(i) = nn - ((2*nxe)+1) + i
      IF(i==1.OR.i==((2*nxe)+1)) THEN
        val(count) = 25._iwp/12._iwp 
      ELSE IF(mod(i,2)==0) THEN
        val(count) = 25._iwp/3._iwp
      ELSE
        val(count) = 25._iwp/6._iwp
      END IF
      count = count + 1
    END DO

!------------------------------------------------------------------------------
! p129.5 Select data format
!------------------------------------------------------------------------------

    SELECT CASE(iotype)

!------------------------------------------------------------------------------
! p129.6 Output in proprietary ParaFEM format (based on Kidger's DanFE)
!------------------------------------------------------------------------------

      CASE('parafem')

        OPEN(11,FILE=argv(1:nlen)//'.d',STATUS='REPLACE',ACTION='WRITE')
  
        WRITE(11,'(A)') "*THREE_DIMENSIONAL"
        WRITE(11,'(A)') "*NODES"

        DO i = 1, nn
          WRITE(11,'(I12,3F12.4)') i, g_coord(:,i)
        END DO

        WRITE(11,'(A)') "*ELEMENTS"
        DO iel = 1, nels
          WRITE(11,'(I10,A,20I12,A)') iel, " 3 20 1 ", g_num(3,iel),         &
                   g_num(5,iel),g_num(7,iel),g_num(1,iel),g_num(15,iel),     &
                   g_num(17,iel),g_num(19,iel),g_num(13,iel),g_num(4,iel),   &
                   g_num(6,iel),g_num(8,iel),g_num(2,iel),g_num(16,iel),     &
                   g_num(18,iel),g_num(20,iel),g_num(14,iel),g_num(10,iel),   &
                   g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1 "
        END DO

        CLOSE(11)  

        OPEN(12,FILE=argv(1:nlen)//'.bnd',STATUS='REPLACE',ACTION='WRITE')

        DO i = 1, nr
          WRITE(12,'(4I6)') rest(i,:) 
        END DO

        CLOSE(12)

        OPEN(13,FILE=argv(1:nlen)//'.lds',STATUS='REPLACE',ACTION='WRITE')

        DO i = 1, loaded_freedoms
          WRITE(13,'(I12,2A,E16.8)')  no(i),"  0.00000000E+00  ",             &
                                      "  0.00000000E+00  ", val(i)
        END DO

        CLOSE(13)

!------------------------------------------------------------------------------
! p129.7 New control data
!------------------------------------------------------------------------------
 
        OPEN(14,FILE=argv(1:nlen)//'.dat',STATUS='REPLACE',ACTION='WRITE')

        meshgen     = 2 ! current default
        partitioner = 1 ! current default
  
        WRITE(14,'(A)')             "'hexahedron'"
        WRITE(14,'(2I4)')           meshgen,partitioner
        WRITE(14,'(3I12,2I5,3I12)') nels,nn,nr,nip,nod,loaded_freedoms,nres
        WRITE(14,'(5E14.6)')        rho,e,v,alpha1,beta1
        WRITE(14,'(2I6,3E14.6,I6)') nstep,npri,theta,omega,tol,limit

        CLOSE(14)

!------------------------------------------------------------------------------
! p129.8 Output in ASCII Ensight Gold format
!------------------------------------------------------------------------------

      CASE('paraview')

        ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 
        etype=0; nf=0; oldlds=0.0_iwp

        dtim=1.0_iwp
        solid=.true. ! need to review this 

        CALL rest_to_nf(rest,nf)

        CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,              &
                       oldlds(1:),nstep,npri,dtim,solid)

      CASE DEFAULT

        PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

    END SELECT 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program xx3
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('xx3')
  
    READ(10,*) iotype,nels,nxe,nze,nod,nip,aa,bb,cc,e,v,tol,limit

    nye = nels/nxe/nze

!------------------------------------------------------------------------------
! xx3.1 Select 8 node or 20 node hexahedra
!------------------------------------------------------------------------------

    SELECT CASE(nod)
  
    CASE(20)
   
      nr    = ((2*nxe+1)*(nze+1)+(nxe+1)*nze)*2 +                             &
              ((2*nye-1)*nze+(nye-1)*nze)*2     +                             &
               (2*nye-1)*(nxe+1)+(nye-1)*nxe

      nn    = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1) +                   &
               (nxe+1)*(nze+1)*nye

      ndim  = 3
      nodof = 3
      nle   = nxe/5
      
      fixed_freedoms  = 0
      loaded_freedoms = 3*nle*nle + 4*nle + 1
      element         = 'hexahedron'

      ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),num(nod),     &
               rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms))
            
      coord = 0.0_iwp ; g_coord = 0.0_iwp ; val = 0.0_iwp
      g_num = 0       ; rest    = 0       ; no  = 0       ; num = 0
          
!------------------------------------------------------------------------------
! xx3.2  Find nodal coordinates and element steering array
!         Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

      DO iel = 1, nels
        CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
        g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
      END DO
          
!------------------------------------------------------------------------------
! xx3.3  Boundary conditions
!------------------------------------------------------------------------------
            
      CALL cube_bc20(rest,nxe,nye,nze)
          
!------------------------------------------------------------------------------
! xx3.4  Loading conditions
!------------------------------------------------------------------------------
            
      CALL load_p121(nle,nod,nxe,nze, no,val)
      val = -val * aa * bb * (25._iwp / 12._iwp)
              
!------------------------------------------------------------------------------
! xx3.5  Create input deck for 8 node hexahedra
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

        ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),            &
                 rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),   &
                 num(nod))
            
        coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
        g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0

!------------------------------------------------------------------------------
! xx3.6  Find nodal coordinates and element steering array
!         write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

        DO iel = 1, nels
          CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
          g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
        END DO

!------------------------------------------------------------------------------
! xx3.7  Boundary conditions
!------------------------------------------------------------------------------
         
        CALL cube_bc8(rest,nxe,nye,nze)
          
!------------------------------------------------------------------------------
! xx3.8  Loading conditions
!------------------------------------------------------------------------------
         
        CALL load_p121(nle,nod,nxe,nze, no,val)
        val = val * aa * bb 

!------------------------------------------------------------------------------
! xx3.9  Default case and error message
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
! xx3.10  Output model
!------------------------------------------------------------------------------
         
        SELECT CASE(iotype)

        CASE('parafem')

!------------------------------------------------------------------------------
! xx3.11  Output geometry file ".d" 
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
! xx3.12  Output ".bnd" file
!-----------------------------------------------------------------------------

          OPEN(12,FILE=argv(1:nlen)//'.bnd',STATUS='REPLACE',ACTION='WRITE')
          
          DO i = 1, nr
            WRITE(12,'(I15,3I6)') rest(i,:) 
          END DO
          
          CLOSE(12)

!------------------------------------------------------------------------------
! xx3.13  Output ".lds" file
!------------------------------------------------------------------------------

          OPEN(13,FILE=argv(1:nlen)//'.lds',STATUS='REPLACE',ACTION='WRITE')
          
          DO i = 1, loaded_freedoms
            WRITE(13,'(I12,2A,3E16.8)') no(i),"  0.00000000E+00  ",           &
                                        "0.00000000E+00",val(i) 
          END DO
            
          CLOSE(13)

!------------------------------------------------------------------------------
! xx3.14  Output new ".dat" file
!------------------------------------------------------------------------------

          OPEN(14,FILE=argv(1:nlen)//'.dat',STATUS='REPLACE',ACTION='WRITE')
          
          WRITE(14,'(A)') "'gpu'"
          WRITE(14,'(A)') "'hexahedron'"
          IF(nod==8) THEN
            WRITE(14,'(A)') "1"            ! SGM node numbering scheme
          ELSE
            WRITE(14,'(A)') "2"            ! Abaqus node numbering scheme
          END IF
          WRITE(14,'(A)') "1"              ! Internal mesh partitioning
          WRITE(14,'(3I12,2I5,2I9)')nels,nn,nr,nip,nod,loaded_freedoms
          WRITE(14,'(3E12.4,I8)') e,v,tol,limit 
          
          CLOSE(14)
             
          PRINT *, "Program xx3 can run in cpu only mode or offload to gpu"
          PRINT *, "To select cpu only or gpu offload, edit the dat file"
          PRINT *, "Two options are permitted: 'cpu' or 'gpu'"
          PRINT *

!------------------------------------------------------------------------------
! xx3.15  Output for visualization in ParaView (Ensight Gold)
!------------------------------------------------------------------------------

          CASE('paraview')

            ! modify mesh_ensi for optional arguments

            ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 

            etype  = 1   ! Only one material type in this mesh
            nf     = 0
            oldlds = zero

            nstep=1; npri=1; dtim=1.0; solid=.true. 

            k=0 
            DO j=1,loaded_freedoms
              k=k+1
              found=.false.
              DO i=1,nn
                IF(i==no(k)) THEN
                  l=i*3
                  oldlds(l)=val(k)
                  found=.true.
                END IF
                IF(found)CYCLE
              END DO
            END DO

            CALL rest_to_nf(rest,nf)

            CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,          &
                           oldlds(1:),nstep,npri,dtim,solid)

          CASE DEFAULT

            PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

          END SELECT 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Program xx9
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('xx9')
  
    READ(10,*) iotype,nels,nxe,nze,nod,nip,aa,bb,cc,e,v,tol,limit

    nye = nels/nxe/nze

!------------------------------------------------------------------------------
! xx9.1 Select 8 node or 20 node hexahedra
!------------------------------------------------------------------------------

    SELECT CASE(nod)
  
    CASE(20)
   
      nr    = ((2*nxe+1)*(nze+1)+(nxe+1)*nze)*2 +                             &
              ((2*nye-1)*nze+(nye-1)*nze)*2     +                             &
               (2*nye-1)*(nxe+1)+(nye-1)*nxe

      nn    = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1) +                   &
               (nxe+1)*(nze+1)*nye

      ndim  = 3
      nodof = 3
      nle   = nxe/5
      
      fixed_freedoms  = 0
      loaded_freedoms = 3*nle*nle + 4*nle + 1
      element         = 'hexahedron'

      ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),num(nod),     &
               rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms))
            
      coord = 0.0_iwp ; g_coord = 0.0_iwp ; val = 0.0_iwp
      g_num = 0       ; rest    = 0       ; no  = 0       ; num = 0
          
!------------------------------------------------------------------------------
! xx9.2  Find nodal coordinates and element steering array
!         Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

      DO iel = 1, nels
        CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
        g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
      END DO
          
!------------------------------------------------------------------------------
! xx9.3  Boundary conditions
!------------------------------------------------------------------------------
            
      CALL cube_bc20(rest,nxe,nye,nze)
          
!------------------------------------------------------------------------------
! xx9.4  Loading conditions
!------------------------------------------------------------------------------
            
      CALL load_p121(nle,nod,nxe,nze, no,val)
      val = -val * aa * bb * (25._iwp / 12._iwp)
              
!------------------------------------------------------------------------------
! xx9.5  Create input deck for 8 node hexahedra
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

        ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),            &
                 rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),   &
                 num(nod))
            
        coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
        g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0

!------------------------------------------------------------------------------
! xx9.6  Find nodal coordinates and element steering array
!         write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

        DO iel = 1, nels
          CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
          g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
        END DO

!------------------------------------------------------------------------------
! xx9.7  Boundary conditions
!------------------------------------------------------------------------------
         
        CALL cube_bc8(rest,nxe,nye,nze)
          
!------------------------------------------------------------------------------
! xx9.8  Loading conditions
!------------------------------------------------------------------------------
         
        CALL load_p121(nle,nod,nxe,nze, no,val)
        val = val * aa * bb 

!------------------------------------------------------------------------------
! xx9.9  Default case and error message
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
! xx9.10  Output model
!------------------------------------------------------------------------------
         
        SELECT CASE(iotype)

        CASE('parafem')

!------------------------------------------------------------------------------
! xx9.11  Output geometry file ".d" 
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
          
          PRINT *, "Completed writing '.d'   geometry file"
          CLOSE(11)

!------------------------------------------------------------------------------
! xx9.12  Output ".bnd" file
!-----------------------------------------------------------------------------

          OPEN(12,FILE=argv(1:nlen)//'.bnd',STATUS='REPLACE',ACTION='WRITE')
          
          DO i = 1, nr
            WRITE(12,'(I15,3I6)') rest(i,:) 
          END DO
          
          PRINT *, "Completed writing '.bnd' constraints file"
          CLOSE(12)

!------------------------------------------------------------------------------
! xx9.13  Output ".lds" file
!------------------------------------------------------------------------------

          OPEN(13,FILE=argv(1:nlen)//'.lds',STATUS='REPLACE',ACTION='WRITE')
          
          DO i = 1, loaded_freedoms
            WRITE(13,'(I12,2A,3E16.8)') no(i),"  0.00000000E+00  ",           &
                                        "0.00000000E+00",val(i) 
          END DO
            
          PRINT *, "Completed writing '.lds' nodal loads file"
          CLOSE(13)

!------------------------------------------------------------------------------
! xx9.14  Output new ".dat" file
!------------------------------------------------------------------------------

          OPEN(14,FILE=argv(1:nlen)//'.dat',STATUS='REPLACE',ACTION='WRITE')
          
          WRITE(14,'(A)') "'gpu'"
          WRITE(14,'(A)') "'hexahedron'"
          IF(nod==8) THEN
            WRITE(14,'(A)') "1"            ! SGM node numbering scheme
          ELSE
            WRITE(14,'(A)') "2"            ! Abaqus node numbering scheme
          END IF
          WRITE(14,'(A)') "1"              ! Internal mesh partitioning
          WRITE(14,'(3I12,2I5,2I9)')nels,nn,nr,nip,nod,loaded_freedoms
          WRITE(14,'(3E12.4,I8)') e,v,tol,limit 
          
          PRINT *, "Completed writing '.dat' control data file"
          CLOSE(14)
             
          PRINT *
          PRINT *, "Program xx9 can run in cpu only mode or offload to gpu"
          PRINT *, "To select cpu only or gpu offload, edit the dat file"
          PRINT *, "Two options are permitted: 'cpu' or 'gpu'"
          PRINT *

!------------------------------------------------------------------------------
! xx9.15  Output for visualization in ParaView (Ensight Gold)
!------------------------------------------------------------------------------

          CASE('paraview')

            ! modify mesh_ensi for optional arguments

            ALLOCATE(etype(nels),nf(nodof,nn),oldlds(nn*ndim)) 

            etype  = 1   ! Only one material type in this mesh
            nf     = 0
            oldlds = zero

            nstep=1; npri=1; dtim=1.0; solid=.true. 

            k=0 
            DO j=1,loaded_freedoms
              k=k+1
              found=.false.
              DO i=1,nn
                IF(i==no(k)) THEN
                  l=i*3
                  oldlds(l)=val(k)
                  found=.true.
                END IF
                IF(found)CYCLE
              END DO
            END DO

            CALL rest_to_nf(rest,nf)

            CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,          &
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
    PRINT*, "  p128, p129, p1210, xx3 and xx9"
    PRINT*

  END SELECT
          
END PROGRAM p12meshgen
