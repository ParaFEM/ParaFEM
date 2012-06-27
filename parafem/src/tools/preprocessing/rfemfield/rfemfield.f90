PROGRAM rfemfield

!/****h* tools/preprocessing/rfem
!*  NAME
!*    PROGRAM rfemfield
!*  FUNCTION
!*    Creates box of finite elements and randomly generated field of elastic
!*    properties. Applies field to loaded model and rewrites model and
!*    material files.
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
  INTEGER                :: iel,i,j,k
  INTEGER                :: nels,nxe,nye,nze,ndim,nod,nn
  INTEGER                :: iefld,istat,kseed,output
  INTEGER                :: np_types=2,nprops
  REAL(iwp)              :: aa,bb,cc
  REAL(iwp)              :: eavg,egeo,ehrm
  REAL(iwp)              :: emn,esd
  REAL(iwp)              :: thx,thy,thz
  REAL(iwp),PARAMETER    :: zero=0.0_iwp
  REAL(iwp)              :: v
  CHARACTER(LEN=6)       :: varfnc
  CHARACTER(LEN=15)      :: rfield,job,sub1,sub2
  CHARACTER(LEN=100)     :: model_job_name,rfem_job_name,fname,program_name
  CHARACTER(LEN=100)     :: dat_name,d_name
  CHARACTER(LEN=100)     :: instance_id_arg,instance_id
  LOGICAL                :: debug=.false.
  LOGICAL                :: shofld=.false.
  LOGICAL                :: lunif=.false.
  LOGICAL                :: dcheck=.false.

  LOGICAL                :: do_model=.false.
  CHARACTER(LEN=15)      :: m_element
  INTEGER                :: m_nels,m_nn,m_nr,m_nod,m_nip,m_loaded_nodes
  INTEGER                :: m_limit,m_mesh,m_fixed_freedoms,m_partition
  INTEGER                :: m_np_types 
  REAL(iwp)              :: m_tol,m_mises
  INTEGER                :: ord, ind, idummy
  REAL(iwp)              :: min_extents(3),max_extents(3),size_extents(3)
  REAL(iwp)              :: centroid(3)

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------   
  INTEGER, ALLOCATABLE   :: num(:),g_num(:,:) 
  INTEGER, ALLOCATABLE   :: ieplt(:)
  REAL(iwp), ALLOCATABLE :: coord(:,:),g_coord(:,:)
  REAL(iwp), ALLOCATABLE :: efld(:,:,:)   ! elastic modulus random field

  INTEGER, ALLOCATABLE   :: m_num(:,:) 
  REAL(iwp), ALLOCATABLE :: m_coord(:,:)

!------------------------------------------------------------------------------
! 3. Read model_job_name and rfem_job_name from the command line
!------------------------------------------------------------------------------

  argc = iargc()
  IF( (argc /= 3) .and. (argc /= 1) ) THEN
     PRINT*
     PRINT*, "Usage:  rfemfield <rfem_job_name> <model_job_name> <instance-id>"
     PRINT*, "   Or:  rfemfield <rfem_job_name>"
     PRINT*
     PRINT*, "        program expects as input:"
     PRINT*, "          <rfem_job_name>.rf"
     PRINT*, "          <model_job_name>.dat"
     PRINT*, "          <model_job_name>.d"
     PRINT*
     PRINT*, "        and outputs:"
     PRINT*, "          <model_job_name>-<instance-id>.dat" 
     PRINT*, "          <model_job_name>-<instance-id>.d" 
     PRINT*, "          <model_job_name>-<instance-id>.mat"
     PRINT*
     PRINT*, "        and optionally outputs:"
     PRINT*, "          <rfem_job_name>.d" 
     PRINT*, "          <rfem_job_name>.mat" 
     PRINT*
     STOP
  END IF
  CALL GETARG(1, rfem_job_name)
  IF( argc == 3 ) THEN
     CALL GETARG(2, model_job_name)
     do_model = .true.
     CALL GETARG(3, instance_id_arg)
     READ(instance_id_arg,*) instance_id
     instance_id = instance_id_arg(1:LEN_TRIM(instance_id_arg))
  END IF
  
!------------------------------------------------------------------------------
! 4. Read control data
!------------------------------------------------------------------------------

  IF (INDEX(rfem_job_name,".rf") /= 0) THEN
     rfem_job_name = rfem_job_name(1:INDEX(rfem_job_name,".rf")-1)
  END IF
  fname = rfem_job_name(1:INDEX(rfem_job_name," ")-1) // ".rf"
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
     READ(10,*) v
     
     nye    = nels/nxe/nze
     nn     = (nxe+1)*(nye+1)*(nze+1)
     nprops = nels
     ndim   = 3
     nod    = 8
     
     IF(output /= 1 .AND. output /=2) THEN
        PRINT *, output, " = Incorrect value for output. Accepted values are &
             &'1' or '2'"
        STOP
     END IF
     
     fname = rfem_job_name(1:INDEX(rfem_job_name," ")-1) // ".res"
     OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

!------------------------------------------------------------------------------
! 6. Read model dat, coords and elems; determine extents
!------------------------------------------------------------------------------

     IF( do_model ) THEN
     
        IF (INDEX(model_job_name,".dat") /= 0) THEN
           model_job_name = model_job_name(1:INDEX(model_job_name,".dat")-1)
        END IF
        dat_name = model_job_name(1:INDEX(model_job_name," ")-1) // ".dat"
        d_name = model_job_name(1:INDEX(model_job_name," ")-1) // ".d"
        
        OPEN (14, file=dat_name, status='old', action='read')
        READ(14,*) m_element,m_mesh,m_partition,m_np_types,      &
             m_nels,m_nn,m_nr,m_nip,m_nod,m_loaded_nodes,        &
             m_fixed_freedoms,m_tol,m_limit,m_mises
        CLOSE(14)
        
        ndim   = 3
        ALLOCATE(m_coord(ndim,m_nn))
        ALLOCATE(m_num(m_nod,m_nels))
        
        OPEN (15, file=d_name, status='old', action='read')
        READ(15,*)   !headers
        READ(15,*)   !headers
        READ(15,*) idummy,m_coord(:,1)
        DO ord = 1,3
           min_extents(ord) = m_coord(ord,1)
           max_extents(ord) = m_coord(ord,1)
        END DO
        DO i = 2,m_nn
           READ(15,*) idummy,m_coord(:,i)
           DO ord = 1,3
              IF ( m_coord(ord,i) < min_extents(ord) ) THEN
                 min_extents(ord) = m_coord(ord,i)
              END IF
              IF ( m_coord(ord,i) > max_extents(ord) ) THEN
                 max_extents(ord) = m_coord(ord,i)
              END IF
           END DO
        END DO
        READ(15,*)   !headers
        DO iel = 1,m_nels
           READ(15,*) idummy,idummy,idummy,idummy,m_num(:,iel),idummy
        END DO
        CLOSE(15)
        
        WRITE(11, '(A,I8)') "Number of model nodes = ", m_nn
        WRITE(11, '(A,I8)') "Number of model elems = ", m_nels
        WRITE(11, '(A,I8)') "Number of model nods = ", m_nod
        
        size_extents = max_extents - min_extents
        WRITE(11, '(A)') "Extents of model: (ord,min,max)"
        DO ord = 1,3
           WRITE(11,'(I1,2F14.6)') ord, min_extents(ord), max_extents(ord)
        END DO
        WRITE(11, '(A)') "Size of Extents of model:"
        WRITE(11,'(3F14.6)') size_extents
     END IF
     
!------------------------------------------------------------------------------
! 7. Generate spatially random field for Young's modulus in a regular cuboid
!------------------------------------------------------------------------------
    
     ALLOCATE(efld(nxe,nye,nze))
     ALLOCATE(ieplt(3))
     
     lunif  = .false.  ! not used
     shofld = .false.  ! not used 
     ieplt  = 0        ! not used
     job    = ""       ! not used
     sub1   = ""       ! not used
     sub2   = ""       ! not used 
     debug  = .true.   ! not used
     dcheck = .false.  ! not used
     
     istat  = 11       ! output debugging info to .res file
     efld   = zero
     kseed  = 0        ! important variable
     
     fname = rfem_job_name(1:INDEX(rfem_job_name," ")-1) // ".mat"
     OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
     
     CALL sim3de(efld,nxe,nye,nze,aa,bb,cc,thx,thy,thz,emn,esd,lunif,kseed, &
          shofld,ieplt,iefld,job,sub1,sub2,varfnc,debug,istat,dcheck,   &
          nxe,nye,nze,eavg,egeo,ehrm)
     
     WRITE(11,'(A,E12.4)') "eavg =", eavg
     WRITE(11,'(A,E12.4)') "egeo =", egeo
     WRITE(11,'(A,E12.4)') "Harmonic average =", ehrm
     
     iel = 0
     
     PRINT *, nprops
     PRINT *, np_types
     PRINT *, "MATERIAL"
     
     WRITE(12,'(A,2I8)') "MATERIAL", nprops, np_types
     WRITE(12,'(3A)') "E"," ", "v"
     
     DO i = 1, nze
        DO j = 1, nye
           DO k = 1, nxe
              iel = iel + 1
              WRITE(12,'(I8,2E12.4)') iel, efld(k,j,i), v
           END DO
        END DO
     END DO
     
     CLOSE(11)
     CLOSE(12)

     DEALLOCATE(ieplt)
     
!------------------------------------------------------------------------------
! 8. Generate finite element mesh for a regular cuboid
!------------------------------------------------------------------------------
    
     IF(output == 2 ) THEN
        
        ALLOCATE(coord(nod,ndim),g_coord(ndim,nn),num(nod),g_num(nod,nels))
        
        coord = zero ; g_coord = zero
        num   = 0    ; g_num   = 0
        
        DO iel = 1, nels
           CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,g_num(:,iel))
           g_coord(:, g_num(:,iel)) = TRANSPOSE(coord)
        END DO
        
        fname = rfem_job_name(1:INDEX(rfem_job_name," ")-1) // ".d"
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

     END IF
     
!------------------------------------------------------------------------------
! 9. Rewrite model dataset with new material IDs and material file
!------------------------------------------------------------------------------

     IF( do_model ) THEN

        fname = model_job_name(1:INDEX(model_job_name," ")-1) // "-" // &
             instance_id_arg(1:LEN_TRIM(instance_id_arg))// ".d"
        OPEN(16,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
        
        fname = model_job_name(1:INDEX(model_job_name," ")-1) // "-" // &
             instance_id_arg(1:LEN_TRIM(instance_id_arg)) // ".mat"
        OPEN(17,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
        
        WRITE(16,'(A)') "*THREE_DIMENSIONAL"
        WRITE(16,'(A)') "*NODES"
        
        WRITE(17,'(A,2I8)') '*MATERIAL', nprops, np_types
        WRITE(17,'(3A)') "E"," ", "v"
     
        DO i = 1,m_nn
           WRITE(16,'(I8,3E14.6)') i, m_coord(:,i)
        END DO
        
        WRITE(16,'(A)') "*ELEMENTS"

        ! TODO: m_nod is fixed here; needs to be dynamic
        DO iel = 1, m_nels
           IF( m_nod == 8 ) THEN
              WRITE(16,'(I8,A,8I12,I12)') iel, " 3 8 1 ",m_num(:,iel), iel
           ELSE
              WRITE(16,'(I8,A,20I12,I12)') iel, " 3 20 1 ",m_num(:,iel), iel
           END IF

           ! get k,j,i based on centroid of element
           ! sum and average x,y,z
           centroid = m_coord(:,m_num(1,iel))
           DO ind = 2, m_nod
              centroid = centroid + m_coord(:,m_num(ind,iel))
           END DO
           centroid = centroid / m_nod
           ! normalize to rfield grid
           centroid = centroid - min_extents
           centroid = centroid / size_extents
           centroid(1) = centroid(1) * (nxe)
           centroid(2) = centroid(2) * (nye)
           centroid(3) = centroid(3) * (nze)
           ! get i,j,k indices
           centroid = floor(centroid)
           k = min(int(centroid(1)),nxe-1) + 1
           j = min(int(centroid(2)),nye-1) + 1
           i = min(int(centroid(3)),nze-1) + 1
           WRITE(17,'(I8,2E12.4)') iel, efld(k,j,i), v
        END DO
        
        CLOSE(16)
        CLOSE(17)

     END IF
     
!------------------------------------------------------------------------------
! 10. Write new -rfem .dat file for RFEMSOLVE with updated np_types(m_nels)
!------------------------------------------------------------------------------

     IF( do_model ) THEN

        fname = model_job_name(1:INDEX(model_job_name," ")-1) // "-" &
             // instance_id_arg(1:LEN_TRIM(instance_id_arg)) // ".dat"
        OPEN (14, file=fname, status='REPLACE', action='WRITE')
        WRITE(14,*) m_element
        WRITE(14,*) m_mesh
        WRITE(14,*) m_partition
        WRITE(14,*) m_nels
        WRITE(14, '(6I8)') m_nels, m_nn, m_nr, m_nip, m_nod, m_loaded_nodes, m_fixed_freedoms
        WRITE(14, '(E14.6,I8,E14.6)') m_tol, m_limit, m_mises
        CLOSE(14)

     END IF
  
!------------------------------------------------------------------------------
! 10. Cleanup
!------------------------------------------------------------------------------

     DEALLOCATE(efld)
     IF( do_model ) THEN
        DEALLOCATE(m_coord,m_num)
     END IF

  CASE DEFAULT
     
     PRINT *
     PRINT *, "  *** Program aborted ***"
     PRINT *
     PRINT *, "  Wrong value given in variable: FIELD_TYPE"
     PRINT *, "  The only accepted value is:    SIM3DE"
     PRINT *
     
  END SELECT
  
  CLOSE(10)
  
END PROGRAM rfemfield
