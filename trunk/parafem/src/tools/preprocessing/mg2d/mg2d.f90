PROGRAM mg2d

!/****h* tools/preprocessing/mg2d
!*  NAME
!*    PROGRAM mg2d
!*  FUNCTION
!*    Creates simple meshes using files with the .mg extension.
!*  AUTHOR
!*    Lee Margetts
!*  COPYRIGHT
!*    (c) University of Manchester 2007-2012
!****
!*/

  USE precision
  USE geometry
  USE loading
  
  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER                :: nr,nn,nels,nres,nle,nxe,nye,nze
  INTEGER                :: nod,ndim,nodof,nip 
  INTEGER                :: nmodes,lalfa,leig,lx,lz
  INTEGER                :: i,j,iel,l,m,n,iargc,argc
  INTEGER                :: cjits,limit,maxitr
  INTEGER                :: ell
  INTEGER                :: fixed_nodes,fixed_freedoms,loaded_freedoms
  INTEGER                :: nev,ncv
  INTEGER                :: meshgen,partitioner
  INTEGER                :: nstep,npri,count
  INTEGER                :: np_types
  INTEGER                :: prnwidth,remainder
  REAL(iwp),PARAMETER    :: zero = 0.0_iwp, twelth = 1.0_iwp/12.0_iwp
  REAL(iwp)              :: aa,bb,cc
  REAL(iwp)              :: kx,ky,kz
  REAL(iwp)              :: cjtol  
  REAL(iwp)              :: rho,visc,e,v   
  REAL(iwp)              :: el,er,acc  
  REAL(iwp)              :: kappa,penalty,tol,x0  
  REAL(iwp)              :: alpha1,beta1,theta,omega,dtim
  REAL(iwp)              :: zbox,zele,z1,z3,ttl
  REAL(iwp)              :: val0
  REAL(iwp)              :: sbary
  CHARACTER(LEN=15)      :: element
  CHARACTER(LEN=50)      :: program_name,job_name,fname,problem_type,iotype
  CHARACTER(LEN=1)       :: bmat
  CHARACTER(LEN=2)       :: which

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------   
  INTEGER, ALLOCATABLE   :: g_num(:,:),rest(:,:),nf(:,:),no(:),no_f(:),num(:)
  INTEGER, ALLOCATABLE   :: matID(:)
  REAL(iwp), ALLOCATABLE :: g_coord(:,:),coord(:,:),val(:),val_f(:)

!------------------------------------------------------------------------------
! 3. Read job_name from the command line
!------------------------------------------------------------------------------

  argc = iargc()
  IF (argc /= 1) THEN
    PRINT*
    PRINT*, "Usage:  mg2d <job_name>"
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

  IF (INDEX(job_name,".mg") /= 0) THEN
    job_name = job_name(1:INDEX(job_name,".mg")-1)
  END IF
  fname = job_name(1:INDEX(job_name," ")-1) // ".mg"
  OPEN (10, file=fname, status='old', action='read')

  READ(10,*) program_name

  PRINT *
  PRINT *, "Generating input deck for program ", program_name
  PRINT *

!------------------------------------------------------------------------------
! 5. Select program using the SELECT CASE CONSTRUCT
!------------------------------------------------------------------------------

  SELECT CASE(program_name)

  CASE('p56')

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! x. Program p56 in 4th edition
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    READ(10,*) nxe,nye,nze

    aa              =  0.5_iwp/nxe
    bb              =  3.0_iwp/nye
    cc              = -2.0_iwp/nze
    tol             =  0.00001_iwp
    limit           =  2000
    np_types        =  2
    ndim            =  3
    nodof           =  3
    nod             =  20
    nip             =  8
    nels            =  nxe*nye*nze
    nn              =  (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1)+(nxe+1)   &
                        *(nze+1)*nye
    nr              =  ((2*nxe+1)*(nze+1))+((nxe+1)*nze)+((2*nxe+1)*nye)     &
                       +((nxe+1)*nye)+(nze*nye)+(2*nze*nye)

    loaded_freedoms = 3*nxe*nxe + 4*nxe + 1

    ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),                 &
             rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),        &
             num(nod),matID(nels))
    
    coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
    g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0
    matID    = 0

    DO iel = 1, nels
      CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
    END DO
 
!------------------------------------------------------------------------------
! x. Assign material number to element. Method based on z ordinate. 
!    Only works for this simple box mesh.
!------------------------------------------------------------------------------ 

    zbox = abs((nze*cc) * 0.5_iwp)
    zele = zero 

    DO iel = 1, nels
      z1   = abs(g_coord(3,g_num(1,iel)))
      z3   = abs(g_coord(3,g_num(3,iel)))
      zele = (z1+z3) * 0.5_iwp
      PRINT *, "iel ", iel, "zele ", zele, "zbox ", zbox
      IF(zele <= zbox) THEN
        matID(iel) = 1
      ELSE
        matID(iel) = 2
      END IF
    END DO

!------------------------------------------------------------------------------
! x. Write out mesh
!------------------------------------------------------------------------------

    fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
    
    WRITE(11,'(A)') "*THREE_DIMENSIONAL"
    WRITE(11,'(A)') "*NODES"
  
    DO i = 1, nn
      WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
    END DO
  
    WRITE(11,'(A)') "*ELEMENTS"
    
    DO iel = 1, nels
      WRITE(11,'(I8,A,21I8)')    iel," 3 20 1 ",g_num(3,iel),g_num(5,iel),   &
                                 g_num(7,iel),g_num(1,iel),g_num(15,iel),    &
                                 g_num(17,iel),g_num(19,iel),g_num(13,iel),  &
                                 g_num(4,iel),g_num(6,iel),g_num(8,iel),     &
                                 g_num(2,iel),g_num(16,iel),g_num(18,iel),   &
                                 g_num(20,iel),g_num(14,iel),g_num(10,iel),  &
                                 g_num(11,iel),g_num(12,iel),g_num(9,iel),   &
                                 matID(iel) 
    END DO
  
    CLOSE(11)

!------------------------------------------------------------------------------
! x Boundary conditions
!------------------------------------------------------------------------------
    
    fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
    OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
    CALL cube_bc20_p56(rest,nxe,nye,nze)
  
    DO i = 1, nr
      WRITE(12,'(I12,3I2)') rest(i,:) 
    END DO
  
    CLOSE(12)

!------------------------------------------------------------------------------
! x. Write out loads
!------------------------------------------------------------------------------
   
    fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
    OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

    nle = nxe  
    CALL load_p121(nle,nod,nxe,nze, no,val)
    val = -(val * twelth) * aa * bb

    ttl = zero
    DO i = 1, loaded_freedoms
      ttl = ttl + val(i)
    END DO

    PRINT *, "Total load = ", ttl

    DO i = 1, loaded_freedoms
       WRITE(13,'(I10,2A,3E16.8)') no(i),                                   &
                                  "  0.00000000E+00  ","0.00000000E+00",    &
                                   val(i) 
    END DO

    CLOSE(13)

!------------------------------------------------------------------------------
! x. Write out the dat file in the 4th edition style
!------------------------------------------------------------------------------

    fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
    OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

    WRITE(14,'(I8)') nod
    WRITE(14,'(3I5,I3,ES12.2,I5,I2)') nxe,nye,nze,nip,tol,limit,np_types 
    WRITE(14,'(A)') "100  0.3"
    WRITE(14,'(A)') "100  0.3"

    prnwidth  = 24
    remainder = MOD(nels,prnwidth)
    n         = (nels-remainder)/prnwidth

    IF(nels<=prnwidth) THEN
      DO i = 1,nels
        IF(i==nels) THEN
          WRITE(14,'(I3)') matID(i)
        ELSE
          WRITE(14,'(I3)',ADVANCE='no') matID(i)
        END IF
      END DO 
    ELSE
      DO i = 1,n
        m = ((i-1)*prnwidth) + 1
        l = ((i-1)*prnwidth) + prnwidth
        WRITE(14,'(24I3)') (matID(iel),iel=m,l)
      END DO
  
      m = (n*prnwidth)+1
      l = (n*prnwidth)+remainder

      DO i = m,l
        IF(i==l) THEN
          WRITE(14,'(I3)') matID(i)
        ELSE
          WRITE(14,'(I3)',ADVANCE='no') matID(i)
        END IF
      END DO 

    END IF

! x_coords

    prnwidth  = 8
    remainder = MOD(nxe,prnwidth)
    n         = (nxe-remainder)/prnwidth

    IF(nxe<=prnwidth) THEN
      WRITE(14,'(A)',ADVANCE='no') "  0.0000"
      DO i=1,nxe
        IF(i==nxe) THEN
          WRITE(14,'(F8.4)') (0.5_iwp/nxe)*i
        ELSE
          WRITE(14,'(F8.4)',ADVANCE='no') (0.5_iwp/nxe)*i
        END IF
      END DO
    ELSE
      m = 1
      l = prnwidth
      WRITE(14,'(A)',ADVANCE='no') "  0.0000"
      WRITE(14,'(8F8.4)') (((0.5_iwp/nxe)*iel),iel=m,l)
      prnwidth  = 9
      remainder = MOD(nxe,prnwidth)
      n         = (nxe-remainder)/prnwidth
      DO j = 2,n
        m = ((j-1)*prnwidth) 
        l = ((j-1)*prnwidth) + prnwidth - 1
        WRITE(14,'(9F8.4)') (((0.5_iwp/nxe)*iel),iel=m,l)
      END DO  

      m = (n*prnwidth)
      l = (n*prnwidth)+remainder

      DO i=m,l
        IF(i==l) THEN
          WRITE(14,'(F8.4)') (0.5_iwp/nxe)*i
        ELSE
          WRITE(14,'(F8.4)',ADVANCE='no') (0.5_iwp/nxe)*i
        END IF
      END DO

    END IF

! y-coords

    prnwidth  = 8
    remainder = MOD(nye,prnwidth)
    n         = (nye-remainder)/prnwidth

    IF(nye<=prnwidth) THEN
      WRITE(14,'(A)',ADVANCE='no') "  0.0000"
      DO i=1,nye
        IF(i==nye) THEN
          WRITE(14,'(F8.4)') (3.0_iwp/nye)*i
        ELSE
          WRITE(14,'(F8.4)',ADVANCE='no') (3.0_iwp/nye)*i
        END IF
      END DO
    ELSE
      m = 1
      l = prnwidth
      WRITE(14,'(A)',ADVANCE='no') "  0.0000"
      WRITE(14,'(8F8.4)') (((3.0_iwp/nye)*iel),iel=m,l)
      prnwidth  = 9
      remainder = MOD(nye,prnwidth)
      n         = (nye-remainder)/prnwidth
      DO j = 2,n
        m = ((j-1)*prnwidth) 
        l = ((j-1)*prnwidth) + prnwidth - 1
        WRITE(14,'(9F8.4)') (((3.0_iwp/nye)*iel),iel=m,l)
      END DO  

      m = (n*prnwidth)
      l = (n*prnwidth)+remainder

      DO i=m,l
        IF(i==l) THEN
          WRITE(14,'(F8.4)') (3.0_iwp/nye)*i
        ELSE
          WRITE(14,'(F8.4)',ADVANCE='no') (3.0_iwp/nye)*i
        END IF
      END DO

    END IF

! z-coords

    prnwidth  = 8
    remainder = MOD(nze,prnwidth)
    n         = (nze-remainder)/prnwidth

    IF(nze<=prnwidth) THEN
      WRITE(14,'(A)',ADVANCE='no') "  0.0000"
      DO i=1,nze
        IF(i==nze) THEN
          WRITE(14,'(F8.4)') (-2.0_iwp/nze)*i
        ELSE
          WRITE(14,'(F8.4)',ADVANCE='no') (-2.0_iwp/nze)*i
        END IF
      END DO
    ELSE
      m = 1
      l = prnwidth
      WRITE(14,'(A)',ADVANCE='no') "  0.0000"
      WRITE(14,'(8F8.4)') (((-2.0_iwp/nze)*iel),iel=m,l)
      prnwidth  = 9
      remainder = MOD(nze,prnwidth)
      n         = (nze-remainder)/prnwidth
      DO j = 2,n
        m = ((j-1)*prnwidth) 
        l = ((j-1)*prnwidth) + prnwidth - 1
        WRITE(14,'(9F8.4)') (((-2.0_iwp/nze)*iel),iel=m,l)
      END DO  

      m = (n*prnwidth)
      l = (n*prnwidth)+remainder

      DO i=m,l
        IF(i==l) THEN
          WRITE(14,'(F8.4)') (-2.0_iwp/nze)*i
        ELSE
          WRITE(14,'(F8.4)',ADVANCE='no') (-2.0_iwp/nze)*i
        END IF
      END DO

    END IF

! restrained nodes

    WRITE(14,'(I8)') nr

    prnwidth  = 5
    remainder = MOD(nr,prnwidth)
    n         = (nr-remainder)/prnwidth

    IF(nr<=prnwidth) THEN
      DO i = 1,nr
        IF(i==nr) THEN
          WRITE(14,'(I8,3I2)') rest(i,:)
        ELSE
          WRITE(14,'(I8,3I2)',ADVANCE='no') rest(i,:)
        END IF
      END DO
    ELSE
      DO i = 1,n
        m = ((i-1)*prnwidth) + 1
        l = ((i-1)*prnwidth) + prnwidth
        DO j = m,l
          IF(j==l) THEN
            WRITE(14,'(I8,3I2)') rest(j,:)
          ELSE
            WRITE(14,'(I8,3I2)',ADVANCE='no') rest(j,:)
          END IF
        END DO
      END DO
  
      m = (n*prnwidth)+1
      l = (n*prnwidth)+remainder

      DO i = m,l
        IF(i==l) THEN
          WRITE(14,'(I8,3I2)') rest(i,:)
        ELSE
          WRITE(14,'(I8,3I2)',ADVANCE='no') rest(i,:)
        END IF
      END DO
    END IF       

! loads

    WRITE(14,'(I8)') loaded_freedoms

    DO i = 1, loaded_freedoms
       WRITE(14,'(I10,2A,3E16.8)') no(i),                                   &
                                  "  0.00000000E+00  ","0.00000000E+00",    &
                                   val(i) 
    END DO

! fixed freedoms

    WRITE(14,'(A)') " 0"
     
    CLOSE(14)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 7. Program p121
!    
!    problem_type can have the value 'ed4' or 'boussinesq'
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p121')
    problem_type='ed4' ! hardwired
!   READ(10,*) problem_type
    READ(10,*) iotype, nels, nxe, nze, nod, nip
    READ(10,*) aa, bb, cc, e, v
    READ(10,*) tol, limit

    nye   = nels/nxe/nze

!------------------------------------------------------------------------------
! 7.1 Select 8 node or 20 node hexahedra
!------------------------------------------------------------------------------

    SELECT CASE(iotype)

    CASE('parafem')

    SELECT CASE(nod)
   
!------------------------------------------------------------------------------
! 7.2 Create input deck for 20 node hexahedra
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
! 7.21 Allocate dynamic arrays
!------------------------------------------------------------------------------
  
      ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),              &
               rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),     &
               num(nod))
    
      coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
      g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0
  
!------------------------------------------------------------------------------
! 7.22 Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention 
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
       WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
     END DO
  
     WRITE(11,'(A)') "*ELEMENTS"
    
     DO iel = 1, nels
       WRITE(11,'(I12,A,20I12,A)')iel," 3 20 1 ",g_num(3,iel),g_num(5,iel),   &
                                  g_num(7,iel),g_num(1,iel),g_num(15,iel),    &
                                  g_num(17,iel),g_num(19,iel),g_num(13,iel),  &
                                  g_num(4,iel),g_num(6,iel),g_num(8,iel),     &
                                  g_num(2,iel),g_num(16,iel),g_num(18,iel),   &
                                  g_num(20,iel),g_num(14,iel),g_num(10,iel),  &
                                  g_num(11,iel),g_num(12,iel),g_num(9,iel),   &
                                  " 1"
     END DO
  
     CLOSE(11)
  
!------------------------------------------------------------------------------
! 7.23 Boundary conditions
!------------------------------------------------------------------------------
    
     fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
     OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
     CALL cube_bc20(rest,nxe,nye,nze)
  
     DO i = 1, nr
       WRITE(12,'(I12,3I2)') rest(i,:) 
     END DO
  
     CLOSE(12)
  
!------------------------------------------------------------------------------
! 7.24 Loading conditions
!------------------------------------------------------------------------------
  
     fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
     OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
     IF(problem_type == 'ed4') THEN
       CALL load_p121(nle,nod,nxe,nze, no,val)
       val = -val * aa * bb * (25._iwp / 12._iwp)
       DO i = 1, loaded_freedoms
         WRITE(13,'(I12,2A,3E16.8)') no(i),"  0.00000000E+00  ",              &
                                    "0.00000000E+00",val(i) 
       END DO
     ELSE IF (problem_type == 'boussinesq') THEN
       no  = 1
       val = -1.0_iwp
       DO i = 1, loaded_freedoms
         WRITE(13,'(I12,2A,3E16.8)') no(i),"  0.00000000E+00  ",              &
                                    "0.00000000E+00",val(i) 
       END DO
     ELSE
       PRINT *, "Problem type: ", problem_type, " not recognised.            &&
                & No values written to .lds"
     END IF
    
     CLOSE(13)
  
!------------------------------------------------------------------------------
! 7.25 New control data
!------------------------------------------------------------------------------
  
     fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
     OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
     WRITE(14,'(A)') "'hexahedron'"
     WRITE(14,'(A)') "2"                         ! Abaqus node numbering scheme
     WRITE(14,'(A)') "1"                         ! Internal mesh partitioning
     WRITE(14,'(3I12,2I5,2I9)')nels,nn,nr,nip,nod,loaded_freedoms,fixed_freedoms
     WRITE(14,'(3E12.4,I8,A)') e,v,tol,limit," 1" ! S&G partitions
  
     CLOSE(14)

!------------------------------------------------------------------------------
! 7.3 Create input deck for 8 node hexahedra
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
! 7.31 Allocate dynamic arrays
!------------------------------------------------------------------------------

     ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),               &
              rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),      &
              num(nod))
    
     coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
     g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0

!------------------------------------------------------------------------------
! 7.32 Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

     DO iel = 1, nels
       CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
       g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
     END DO
    
     fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
     OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
    
     WRITE(11,'(A)') "*THREE_DIMENSIONAL"
     WRITE(11,'(A)') "*NODES"
  
     DO i = 1, nn
       WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
     END DO
  
     WRITE(11,'(A)') "*ELEMENTS"
    
     DO iel = 1, nels
!      WRITE(11,'(I12,A,8I12,A)') iel, " 3 8 1 ", g_num(1,iel),g_num(4,iel),  &
!                                   g_num(8,iel),g_num(5,iel),g_num(2,iel),   &
!                                   g_num(3,iel),g_num(7,iel),g_num(6,iel),   &
!                                   " 1"
       WRITE(11,'(I12,A,8I12,A)') iel, " 3 8 1 ", g_num(1,iel),g_num(2,iel),  &
                                    g_num(3,iel),g_num(4,iel),g_num(5,iel),   &
                                    g_num(6,iel),g_num(7,iel),g_num(8,iel),   &
                                    " 1"
     END DO
    
     CLOSE(11)

!------------------------------------------------------------------------------
! 7.33 Boundary conditions
!------------------------------------------------------------------------------
  
     fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
     OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
     CALL cube_bc8(rest,nxe,nye,nze)
  
     DO i = 1, nr
       WRITE(12,'(I8,3I6)') rest(i,:) 
     END DO
  
     CLOSE(12)

!------------------------------------------------------------------------------
! 7.34 Loading conditions
!------------------------------------------------------------------------------

     fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
     OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
     IF(problem_type == 'ed4') THEN
       CALL load_p121(nle,nod,nxe,nze, no,val)
       val = val * aa * bb 
       DO i = 1, loaded_freedoms
         WRITE(13,'(I10,2A,3E16.8)') no(i),                                   &
                                    "  0.00000000E+00  ","0.00000000E+00",    &
                                     val(i) 
       END DO
     ELSE IF (problem_type == 'boussinesq') THEN
       no  = 1
       val = -1.0_iwp
       DO i = 1, loaded_freedoms
         WRITE(13,'(I10,2A,3E16.8)') no(i),                                   &
                                    "  0.00000000E+00  ","0.00000000E+00",    &
                                     val(i) 
       END DO
     ELSE
       PRINT *, "Problem type: ", problem_type, " not recognised.            &&
                & No values written to .lds"
     END IF
    
     CLOSE(13)

!------------------------------------------------------------------------------
! 7.35 New control data
!------------------------------------------------------------------------------

     fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
     OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
     WRITE(14,'(A)') "'hexahedron'"
     IF(nod==8) THEN
       WRITE(14,'(A)') "1"            ! Abaqus node numbering scheme
     ELSE
       WRITE(14,'(A)') "2"            ! Abaqus node numbering scheme
     END IF
     WRITE(14,'(A)') "1"              ! Internal mesh partitioning
     WRITE(14,'(7I9)') nels, nn, nr, nip, nod, loaded_freedoms, fixed_freedoms
     WRITE(14,'(3E12.4,I8,A)') e, v, tol, limit, " 1" ! S&G partitions

     CLOSE(14)

!------------------------------------------------------------------------------
! 7.4 Default case and error message
!------------------------------------------------------------------------------
  
    
   CASE DEFAULT
  
     PRINT *
     PRINT *, "Wrong value given in variable NOD"
     PRINT *, "  Accepted values are 8 and 20"
     PRINT *, "  Here NOD = ", nod
     PRINT *
      
   END SELECT
  
   CASE('paraview')

     PRINT *, "  Output for ParaView not yet implemented"; PRINT *, ""

   CASE DEFAULT

     PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

   END SELECT 

!------------------------------------------------------------------------------ 
!------------------------------------------------------------------------------
! 8. Program p122
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p122')

  PRINT*
  PRINT*, "Program p122 not supported"
  PRINT*

  STOP

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 9. Program p123
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p123')

    READ(10,*) nels, nxe, nze, nip
    READ(10,*) aa, bb, cc, kx, ky, kz
    READ(10,*) tol, limit
    READ(10,*) loaded_freedoms, fixed_freedoms
 
    PRINT *, "Read .mg file"

!------------------------------------------------------------------------------
! 9.1 Initialize variables
!------------------------------------------------------------------------------

    nye   = nels/nxe/nze
    ndim  = 3
    nod   = 8
    nr    = (nxe+1)*(nye+1) + (nxe+1)*nze + nye*nze
    nn    = (nxe+1)*(nye+1)*(nze+1)
    nodof = 1
    nres  = nxe*(nze-1)+1

!------------------------------------------------------------------------------
! 9.2 Allocate dynamic arrays
!------------------------------------------------------------------------------
  
    ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),              &
             rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),     &
             num(nod),val_f(fixed_freedoms),no_f(fixed_freedoms))
    
    coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
    g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0
    val_f    = 0.0_iwp ; no_f    = 0
  
!------------------------------------------------------------------------------
! 9.3 Find nodal coordinates and element steering array
!     Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

    DO iel = 1, nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
    END DO
    
    fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
    
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
! 9.4 Boundary conditions
!------------------------------------------------------------------------------
  
    fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
    OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
    CALL box_bc8(rest,nxe,nye,nze)
  
    DO i = 1, nr
      WRITE(12,'(I8,3I6)') rest(i,:) 
    END DO
  
    CLOSE(12)

    PRINT *, "Output boundary conditions"

!------------------------------------------------------------------------------
! 9.5 Loading conditions
!------------------------------------------------------------------------------

    IF(loaded_freedoms > 0) THEN
     
      fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
      OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
     
      no   = nres
      val  = 10.0_iwp
  
      DO i = 1, loaded_freedoms
        WRITE(13,'(I10,E16.8)') no(i),val(i)
      END DO

      CLOSE(13)

      PRINT *, "Output fixed loads"

    END IF

    IF(fixed_freedoms>0) THEN

      fname = job_name(1:INDEX(job_name, " ")-1) // ".fix" 
      OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

      no_f  = nres
      val_f = 100.0_iwp

      DO i = 1, fixed_freedoms
        WRITE(14,'(I10,E16.8)') no_f(i),val_f(i)
      END DO

      CLOSE(14)

      PRINT *, "Output fixed freedoms"

    END IF

!------------------------------------------------------------------------------
! 9.6 New control data
!------------------------------------------------------------------------------

     fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
     OPEN(15,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
     WRITE(15,'(A)') "'hexahedron'"
     WRITE(15,'(A)') "2"              ! Abaqus node numbering scheme
     WRITE(15,'(A)') "1"              ! Internal mesh partitioning
     WRITE(15,'(7I9)') nels, nn, nr, nip, nod, loaded_freedoms, fixed_freedoms
     WRITE(15,'(4E12.4,I8)') kx, ky, kz, tol, limit

     CLOSE(15)

     PRINT *, "Output new control data file"
     PRINT *, "Job completed"
     PRINT *

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 10. Program p124
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  CASE('p124')

    READ(10,*) nels, nxe, nze, nip
    READ(10,*) aa, bb, cc, kx, ky, kz
    READ(10,*) tol, limit
    READ(10,*) loaded_freedoms, fixed_freedoms
 
    PRINT *, "Read .mg file"

!------------------------------------------------------------------------------
! 10.1 Initialize variables
!------------------------------------------------------------------------------

    nye   = nels/nxe/nze
    ndim  = 3
    nod   = 8
    nr    = (nxe+1)*(nye+1) + (nxe+1)*nze + nye*nze
    nn    = (nxe+1)*(nye+1)*(nze+1)
    nodof = 1
    nres  = nxe*(nze-1)+1
    dtim  = 0.01_iwp
    nstep = 150
    theta = 0.5
    npri  = 10

!------------------------------------------------------------------------------
! 10.2 Allocate dynamic arrays
!------------------------------------------------------------------------------
  
    ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),              &
             rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),     &
             num(nod),val_f(fixed_freedoms),no_f(fixed_freedoms))
    
    coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
    g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0
    val_f    = 0.0_iwp ; no_f    = 0
  
!------------------------------------------------------------------------------
! 10.3 Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

    DO iel = 1, nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
      g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
    END DO
    
    fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
    
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
! 10.4 Boundary conditions
!------------------------------------------------------------------------------
  
    fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
    OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
    CALL box_bc8(rest,nxe,nye,nze)
  
    DO i = 1, nr
      WRITE(12,'(I8,3I6)') rest(i,:) 
    END DO
  
    CLOSE(12)

    PRINT *, "Output boundary conditions"

!------------------------------------------------------------------------------
! 10.5 Loading conditions
!------------------------------------------------------------------------------

    IF(loaded_freedoms > 0) THEN
     
      fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
      OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
     
      no   = nres
      val  = 10.0_iwp
  
      DO i = 1, loaded_freedoms
        WRITE(13,'(I10,E16.8)') no(i),val(i)
      END DO

      CLOSE(13)

      PRINT *, "Output fixed loads"

    END IF

    IF(fixed_freedoms>0) THEN

      fname = job_name(1:INDEX(job_name, " ")-1) // ".fix" 
      OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

      no_f  = nres
      val_f = 100.0_iwp

      DO i = 1, fixed_freedoms
        WRITE(14,'(I10,E16.8)') no_f(i),val_f(i)
      END DO

      CLOSE(14)

      PRINT *, "Output fixed freedoms"

    END IF

!------------------------------------------------------------------------------
! 10.6 New control data
!------------------------------------------------------------------------------

     fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
     OPEN(15,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
     WRITE(15,'(A)') "'hexahedron'"
     WRITE(15,'(A)') "2"              ! Abaqus node numbering scheme
     WRITE(15,'(A)') "1"              ! Internal mesh partitioning
     WRITE(15,'(7I9)') nels, nn, nr, nip, nod, loaded_freedoms, fixed_freedoms
     WRITE(15,'(4E12.4,I8)') kx, ky, kz, tol, limit
     WRITE(15,'(E12.4,I8,I8,E12.4)') dtim, nstep, npri, theta

     CLOSE(15)

     PRINT *, "Output new control data file"
     PRINT *, "Some values have default values"
     PRINT *, "Job completed"
     PRINT *

!------------------------------------------------------------------------------
! 11. Program p125
!------------------------------------------------------------------------------

  CASE('p125')

    READ(10,*) iotype, nels, nxe, nze, nip
    READ(10,*) aa, bb, cc, kx, ky, kz
    READ(10,*) dtim, nstep
    READ(10,*) npri, val0
    
    loaded_freedoms=0; fixed_freedoms=0
 
    PRINT *, "Read .mg file"

    SELECT CASE(iotype)

    CASE('parafem')

!------------------------------------------------------------------------------
! 10.1 Initialize variables
!------------------------------------------------------------------------------

    nye   = nels/nxe/nze
    ndim  = 3
    nod   = 8
    nr    = (nxe+1)*(nye+1) + (nxe+1)*nze + nye*nze
    nn    = (nxe+1)*(nye+1)*(nze+1)
    nodof = 1
    nres  = nxe*(nze-1)+1

!------------------------------------------------------------------------------
! 10.2 Allocate dynamic arrays
!------------------------------------------------------------------------------
  
    ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),                              &
!            rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),     &
!            num(nod),val_f(fixed_freedoms),no_f(fixed_freedoms))
             num(nod))
    
 !  coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
 !  g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0
 !  val_f    = 0.0_iwp ; no_f    = 0
   coord    = 0.0_iwp ; g_coord = 0.0_iwp ; num = 0 
  
!------------------------------------------------------------------------------
! 10.3 Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

    DO iel = 1, nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,num)
      g_coord(:,num) = TRANSPOSE(coord)
    END DO
    
    fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
    
    WRITE(11,'(A)') "*THREE_DIMENSIONAL"
    WRITE(11,'(A)') "*NODES"
  
    DO i = 1, nn
      WRITE(11,'(I12,3E14.6)') i, g_coord(:,i)
    END DO

    DEALLOCATE(g_coord)
  
    WRITE(11,'(A)') "*ELEMENTS"
    
    DO iel = 1, nels
      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,num)
      WRITE(11,'(I12,A,8I12,A)') iel, " 3 8 1 ", num(1), num(4),              &
                                 num(8), num(5), num(2),                      &
                                 num(3), num(7), num(6)," 1"
    END DO
    
    CLOSE(11)

    PRINT *, "Output nodal coordinates and element steering array"

!------------------------------------------------------------------------------
! 10.4 Boundary conditions
!------------------------------------------------------------------------------
  
    ALLOCATE(rest(nr,nodof+1))
    rest  = 0

    fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
    OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
    CALL box_bc8(rest,nxe,nye,nze)
  
    DO i = 1, nr
      WRITE(12,'(I12,3I6)') rest(i,:) 
    END DO
  
    CLOSE(12)

    PRINT *, "Output boundary conditions"

!------------------------------------------------------------------------------
! 10.5 Loading conditions
!------------------------------------------------------------------------------

    IF(loaded_freedoms > 0) THEN
     
      fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
      OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
     
      no   = nres
      val  = 10.0_iwp
  
      DO i = 1, loaded_freedoms
        WRITE(13,'(I11,E16.8)') no(i),val(i)
      END DO

      CLOSE(13)

      PRINT *, "Output fixed loads"

    END IF

    IF(fixed_freedoms>0) THEN

      fname = job_name(1:INDEX(job_name, " ")-1) // ".fix" 
      OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

      no_f  = nres
      val_f = 100.0_iwp

      DO i = 1, fixed_freedoms
        WRITE(14,'(I12,E16.8)') no_f(i),val_f(i)
      END DO

      CLOSE(14)

      PRINT *, "Output fixed freedoms"

    END IF

!------------------------------------------------------------------------------
! 10.6 New control data
!------------------------------------------------------------------------------

    fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
    OPEN(15,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
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

      PRINT *, "  Output for ParaView not yet implemented"; PRINT *, ""

    CASE DEFAULT

      PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

    END SELECT 

!------------------------------------------------------------------------------
! 12. Program p126
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


  SELECT CASE(iotype)

  CASE('parafem')
 
  PRINT *, "  Writing nodal coordinates and element steering array"
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

  PRINT *, "  Writing boundary conditions"

  fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
  OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  CALL ns_cube_bc20(rest,nxe,nye,nze)
 
  DO i = 1, nr
    WRITE(12,'(I8,4I6)') rest(i,:) 
  END DO

  CLOSE(12)

!------------------------------------------------------------------------------
! 12.4 Loading conditions for the lid
!------------------------------------------------------------------------------

  PRINT *, "  Writing loading conditions for the lid"

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

  PRINT *, "  Writing new control data file"
 
  fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
  OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

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
! 12.6 Write a .dis file so that boundary conditions can be checked
!------------------------------------------------------------------------------

  PRINT *, "  Writing dummy .dis file with flags for boundary conditions"

  fname = job_name(1:INDEX(job_name, " ")-1) // ".dis" 
  OPEN(15,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

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

  CASE('paraview')

    PRINT *, "  Output for ParaView not yet implemented"; PRINT *, ""

  CASE DEFAULT

    PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

  END SELECT 

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
 
! PRINT*
! PRINT*, "      The <my_job>.mg file should contain the following data:"
! PRINT*
! PRINT*, "      p128ar"
! PRINT*, "      element"
! PRINT*, "      nels, nxe, nze, nod, nip"
! PRINT*, "      aa, bb, cc"
! PRINT*, "      rho, e, v"
! PRINT*, "      nev, ncv, bmat, which"
! PRINT*, "      tol, maxitr"
! PRINT*

  CASE('p128')
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
   
  loaded_freedoms = (2*nxe)+1 ! Differs from the problem in the book

  nres = 3*(nye*(nxe+1)*(nze+1)+nr*(nye-1)+(nxe+1))

!------------------------------------------------------------------------------
! 16.1 Allocate dynamic arrays
!------------------------------------------------------------------------------

  SELECT CASE(iotype)

  CASE('parafem')

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
! 16.2 Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention
!------------------------------------------------------------------------------
 
  PRINT*, "  Writing nodel coordinates and element steering array"

  DO iel=1,nels
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
    WRITE(11,'(I10,A,20I12,A)') iel, " 3 20 1 ", g_num(3,iel),g_num(5,iel),   &
                                 g_num(7,iel),g_num(1,iel),g_num(15,iel),     &
                                 g_num(17,iel),g_num(19,iel),g_num(13,iel),   &
                                 g_num(4,iel),g_num(6,iel),g_num(8,iel),      &
                                 g_num(2,iel),g_num(16,iel),g_num(18,iel),    &
                                 g_num(20,iel),g_num(14,iel),g_num(10,iel),   &
                                 g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1 "
  END DO

  CLOSE(11)  

!------------------------------------------------------------------------------
! 16.3 Boundary conditions
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
  OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  rest = 0
 
  DO i=1,nr
    rest(i,1) = i
  END DO

  DO i = 1, nr
    WRITE(12,'(4I6)') rest(i,:) 
  END DO

  CLOSE(12)

!------------------------------------------------------------------------------
! 16.4 Loading conditions
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

  fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
  OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  DO i = 1, loaded_freedoms
    WRITE(13,'(I12,2A,E16.8)')  no(i),"  0.00000000E+00  ",                  &
                                      "  0.00000000E+00  ", val(i)
  END DO

  CLOSE(13)

!------------------------------------------------------------------------------
! 16.5 New control data
!------------------------------------------------------------------------------
 
  PRINT *, "  Writing new control data file"

  fname = job_name(1:INDEX(job_name, " ")-1) // ".dat"
  OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  meshgen     = 2 ! current default
  partitioner = 1 ! current default
  
  WRITE(14,'(A)')             "'hexahedron'"
  WRITE(14,'(2I4)')           meshgen,partitioner
  WRITE(14,'(3I12,2I5,3I12)') nels,nn,nr,nip,nod,loaded_freedoms,nres
  WRITE(14,'(5E14.6)')        rho,e,v,alpha1,beta1
  WRITE(14,'(2I6,3E14.6,I6)') nstep,npri,theta,omega,tol,limit

  CLOSE(14)

  CASE('paraview')

    PRINT *, "  Output for ParaView not yet implemented"; PRINT *, ""

  CASE DEFAULT

    PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

  END SELECT 

!------------------------------------------------------------------------------
! 16. Program p1210
!------------------------------------------------------------------------------

  CASE('p1210')

  READ(10,*) iotype,nels,nxe,nze,nip
  READ(10,*) aa,bb,cc
  READ(10,*) rho,e,v,sbary
  READ(10,*) dtim,nstep,npri

  nye   = nels/nxe/nze
  nr    = 3*nxe*nze+2*nxe+2*nze+1
  ndim  = 3
  nodof = 3
  nod   = 20
  nn    = (((2*nxe+1)*(nze+1))+((nxe+1)*nze))*(nye+1)+(nxe+1)*(nze+1)*nye
   
! loaded_freedoms = (2*nxe)+1 ! Differs from the problem in the book
  loaded_freedoms = 1 ! Same as problem in the book

  nres = 3*(nye*(nxe+1)*(nze+1)+nr*(nye-1)+(nxe+1))

!------------------------------------------------------------------------------
! 16.1 Allocate dynamic arrays
!------------------------------------------------------------------------------

  SELECT CASE(iotype)

  CASE('parafem')

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
! 16.2 Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention
!------------------------------------------------------------------------------
 
  PRINT*, "  Writing nodel coordinates and element steering array"

  DO iel=1,nels
    CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
    g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)    
  END DO
  
  fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
  OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  WRITE(11,'(A)') "*THREE_DIMENSIONAL"
  WRITE(11,'(A)') "*NODES"

  DO i = 1, nn
    WRITE(11,'(I10,3F12.4)') i, g_coord(:,i)
  END DO

  WRITE(11,'(A)') "*ELEMENTS"
  DO iel = 1, nels
    WRITE(11,'(I10,A,20I10,A)')  iel, " 3 20 1 ", g_num(3,iel),g_num(5,iel),  &
                                 g_num(7,iel),g_num(1,iel),g_num(15,iel),     &
                                 g_num(17,iel),g_num(19,iel),g_num(13,iel),   &
                                 g_num(4,iel),g_num(6,iel),g_num(8,iel),      &
                                 g_num(2,iel),g_num(16,iel),g_num(18,iel),    &
                                 g_num(20,iel),g_num(14,iel),g_num(10,iel),   &
                                 g_num(11,iel),g_num(12,iel),g_num(9,iel)," 1 "
  END DO

  CLOSE(11)  

!------------------------------------------------------------------------------
! 16.3 Boundary conditions
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
  OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  rest = 0
 
  DO i=1,nr
    rest(i,1) = i
  END DO

  DO i = 1, nr
    WRITE(12,'(4I10)') rest(i,:) 
  END DO

  CLOSE(12)

!------------------------------------------------------------------------------
! 16.4 Loading conditions
!------------------------------------------------------------------------------

  count = 1

! DO i = 1,loaded_freedoms
!   no(i) = nn - ((2*nxe)+1) + i
!   IF(i==1.OR.i==((2*nxe)+1)) THEN
!     val(count) = -0.1_iwp/12._iwp 
!   ELSE IF(mod(i,2)==0) THEN
!     val(count) = 0.1_iwp/3._iwp
!   ELSE
!     val(count) = -0.1_iwp/6._iwp
!   END IF
!   count = count + 1
! END DO

  no(1)  = nn-nxe
  val    = 0.0_iwp
  val(1) = -0.001_iwp

  fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
  OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  DO i = 1, loaded_freedoms
    WRITE(13,'(I12,2A,E16.8)')  no(i),"  0.00000000E+00  ",                  &
                                      "  0.00000000E+00  ", val(i)
  END DO

  CLOSE(13)

!------------------------------------------------------------------------------
! 16.5 New control data
!------------------------------------------------------------------------------
 
  PRINT *, "  Writing new control data file"

  fname = job_name(1:INDEX(job_name, " ")-1) // ".dat"
  OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

  meshgen     = 2 ! current default
  partitioner = 1 ! current default
  
  WRITE(14,'(A)')               "'hexahedron'"
  WRITE(14,'(2I4)')             meshgen,partitioner
  WRITE(14,'(I10,I3,5I10)')     nels,nip,nn,nr,nod,loaded_freedoms,nres
  WRITE(14,'(4E14.6)')          rho,e,v,sbary
  WRITE(14,'(E14.6,2I8,A)')     dtim,nstep,npri, " -1.0000"

  CLOSE(14)

  CASE('paraview')

    PRINT *, "  Output for ParaView not yet implemented"; PRINT *, ""

  CASE DEFAULT

    PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

  END SELECT 

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
  
END PROGRAM mg2d
