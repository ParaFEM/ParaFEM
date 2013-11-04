PROGRAM p12meshgen

!/****h* tools/preprocessing/p12meshgen
!*  NAME
!*    PROGRAM p12meshgen
!*  FUNCTION
!*    Creates simple meshes using files with the .mg extension.
!*  AUTHOR
!*    Lee Margetts
!*  COPYRIGHT
!*    (c) University of Manchester 2007-2013
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
  CHARACTER(LEN=50)      :: program_name,job_name,fname,problem_type
  CHARACTER(LEN=50)      :: iotype
  CHARACTER(LEN=1)       :: bmat
  CHARACTER(LEN=2)       :: which

!-------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------- 

  INTEGER, ALLOCATABLE   :: g_num(:,:),rest(:,:),nf(:,:),no(:),no_f(:)
  INTEGER, ALLOCATABLE   :: num(:),matID(:)
  REAL(iwp), ALLOCATABLE :: g_coord(:,:),coord(:,:),val(:),val_f(:),qinc(:)

!-------------------------------------------------------------------------
! 3. Read job_name from the command line
!-------------------------------------------------------------------------

  argc = iargc()
  IF (argc /= 1) THEN
    PRINT*
    PRINT*,"Usage: p12meshgen <job_name>"
    PRINT*
    PRINT*,"       program expects <job_name>.mg and outputs any or all of"
    PRINT*,"       <job_name>.dat <job_name>.d" 
    PRINT*,"       <job_name>.bnd <job_name>.lds"
    PRINT*
    STOP
  END IF
  CALL GETARG(1, job_name)

!-------------------------------------------------------------------------
! 4. Read control data
!-------------------------------------------------------------------------

  IF (INDEX(job_name,".mg") /= 0) THEN
    job_name = job_name(1:INDEX(job_name,".mg")-1)
  END IF
  fname = job_name(1:INDEX(job_name," ")-1) // ".mg"
  OPEN (10, file=fname, status='old', action='read')

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
! 6.1 Program p121
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

  CASE('p121')
  
    READ(10,*) iotype,nels,nxe,nze,nod,nip,aa,bb,cc,e,v,tol,limit

    nye = nels/nxe/nze

!-------------------------------------------------------------------------
! 6.2 Select 8 node or 20 node hexahedra
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

      ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),num(nod),&
               rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms))
            
      coord = 0.0_iwp ; g_coord = 0.0_iwp ; val = 0.0_iwp
      g_num = 0       ; rest    = 0       ; no  = 0       ; num = 0
          
!------------------------------------------------------------------------------
! 6.3  Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

      DO iel = 1, nels
        CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
        g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
      END DO
          
!------------------------------------------------------------------------------
! 6.4  Boundary conditions
!------------------------------------------------------------------------------
            
      CALL cube_bc20(rest,nxe,nye,nze)
          
!------------------------------------------------------------------------------
! 6.5  Loading conditions
!------------------------------------------------------------------------------
            
      CALL load_p121(nle,nod,nxe,nze, no,val)
      val = -val * aa * bb * (25._iwp / 12._iwp)
              
!------------------------------------------------------------------------------
! 6.6 Create input deck for 8 node hexahedra
!------------------------------------------------------------------------------
            
      CASE(8)
          
        nr    = ((nxe+1)*(nze+1))*2 + ((nye-1)*(nze+1))*2 + ((nxe-1)*(nze-1)) 
        ndim  = 3
        nodof = 3
        nn    = (nxe+1)*(nye+1)*(nze+1)
           
        nle             = nxe/5
        loaded_freedoms = (nle+1)*(nle+1)
        fixed_freedoms  = 0

        ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),g_num(nod,nels),           &
                 rest(nr,nodof+1),val(loaded_freedoms),no(loaded_freedoms),  &
                 num(nod))
            
        coord    = 0.0_iwp ; g_coord = 0.0_iwp ;   val = 0.0_iwp
        g_num    = 0       ; rest    = 0       ;   no  = 0       ; num = 0

!------------------------------------------------------------------------------
! 6.7  Find nodal coordinates and element steering array
!      Write to file using Abaqus node numbering convention 
!------------------------------------------------------------------------------

        DO iel = 1, nels
          CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
          g_coord(:,g_num(:,iel)) = TRANSPOSE(coord)
        END DO

!------------------------------------------------------------------------------
! 6.8  Boundary conditions
!------------------------------------------------------------------------------
         
        CALL cube_bc8(rest,nxe,nye,nze)
          
!------------------------------------------------------------------------------
! 6.9  Loading conditions
!------------------------------------------------------------------------------
         
        CALL load_p121(nle,nod,nxe,nze, no,val)
        val = val * aa * bb 

!------------------------------------------------------------------------------
! 6.10 Default case and error message
!------------------------------------------------------------------------------
              
        CASE DEFAULT
        
          PRINT *
          PRINT *, "Wrong value given in variable NOD"
          PRINT *, "  Accepted values are 8 and 20"
          PRINT *, "  Here NOD = ", nod
          PRINT *
             
          iotype = 'aborted'
             
        END SELECT

!-------------------------------------------------------------------------
! 6.11 Output model
!-------------------------------------------------------------------------
         
        SELECT CASE(iotype)

        CASE('parafem')

!-------------------------------------------------------------------------
! 6.12 Output ".d" file
!-------------------------------------------------------------------------
         
          fname = job_name(1:INDEX(job_name, " ")-1) // ".d" 
          OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
            
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

!-------------------------------------------------------------------------
! 6.13 Output ".bnd" file
!-------------------------------------------------------------------------

          fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd" 
          OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
          
          DO i = 1, nr
            WRITE(12,'(I8,3I6)') rest(i,:) 
          END DO
          
          CLOSE(12)

!-------------------------------------------------------------------------
! 6.14 Output ".lds" file
!-------------------------------------------------------------------------

          fname = job_name(1:INDEX(job_name, " ")-1) // ".lds" 
          OPEN(13,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
          
          DO i = 1, loaded_freedoms
            WRITE(13,'(I12,2A,3E16.8)') no(i),"  0.00000000E+00  ",           &
                                        "0.00000000E+00",val(i) 
          END DO
            
          CLOSE(13)

!-------------------------------------------------------------------------
! 6.15 Output new ".dat" file
!-------------------------------------------------------------------------

          fname = job_name(1:INDEX(job_name, " ")-1) // ".dat" 
          OPEN(14,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
          
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
             
!-------------------------------------------------------------------------
! 6.16 ParaView format (Ensight Gold)
!-------------------------------------------------------------------------

          CASE('paraview')

            PRINT *, "  Output for ParaView not yet implemented"; PRINT *, ""

          CASE DEFAULT

            PRINT *, "  Option ", iotype, " not recognised."; PRINT *, ""

          END SELECT 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

          CASE DEFAULT

            PRINT*
            PRINT*, "Mesh only generated for programs: "
            PRINT*, "  p121, p122, p123, p124, p125, p126"
            PRINT*, "  p127, p128, p129 and p1210"
            PRINT*

          END SELECT
          
END PROGRAM p12meshgen
