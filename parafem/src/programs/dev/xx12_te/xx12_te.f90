PROGRAM xx12_te
!------------------------------------------------------------------------------ 
!      Program xx12_te three dimensional analysis of an elastic solid
!                        load control or displacement control; multiple
!                        material types; sequential version
!      
!      Author(s): Llion Evans, David Arregui and Lee Margetts
!
!      Cite DOI : 10.1016/j.nucmat.2015.05.058
!                 10.1007/s11831-014-9139-3
!                 
!------------------------------------------------------------------------------ 
                    ! need to change this to test on multiple cores  
!  USE mpi_wrapper   ! uncomment line to compile without MPI

  USE precision     ; USE global_variables ; USE mp_interface
  USE input         ; USE output           ; USE loading
  USE timing        ; USE maths            ; USE gather_scatter
  USE partition     ; USE elements         ; USE steering
  USE geometry      ; USE pcg              ; USE new_library

  USE, INTRINSIC :: ISO_C_BINDING ! to output C binary file
  
  IMPLICIT NONE

!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

 ! neq,ntot are now global variables - not declared

  INTEGER,PARAMETER     :: nodof=3,ndim=3,nst=6,nprops=8
  INTEGER               :: loaded_nodes,fixed_freedoms,iel,i,j,k,l,idx1,idx2
  INTEGER               :: iters,limit,nn,nr,nip,nod,nels,ndof,npes_pp
  INTEGER               :: node_end,node_start,nodes_pp,int_in
  INTEGER               :: argc,iargc,meshgen,partitioner,np_types
  INTEGER               :: fixed_freedoms_pp,fixed_freedoms_start
  INTEGER               :: nodecount_pp(1) ! count of local nodes
  INTEGER               :: nodecount    ! global count of nodes
  INTEGER               :: prog
  REAL(iwp)             :: e,v,det,tol,up,alpha,beta,tload
!  REAL(iwp)             :: test
  REAL(iwp)             :: mises        ! threshold for mises stress
  REAL(iwp),PARAMETER   :: zero    = 0.0_iwp
  REAL(iwp),PARAMETER   :: penalty = 1.0e20_iwp
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=15)     :: keyword
  CHARACTER(LEN=50)     :: program_name='xx12_te'
  CHARACTER(LEN=50)     :: fname,inst_in,job_name,label,instance_id,stepnum,val0_str
  CHARACTER(LEN=80)     :: cbuffer
  LOGICAL               :: converged = .false.

  !New variables
  !Temporal variable
  REAL(iwp)             :: constant,val0
       
  !Temperature value
  REAL(iwp)             :: gtemp

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:),val(:,:),disp_pp(:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:)
  REAL(iwp),ALLOCATABLE :: bee(:,:),storkm_pp(:,:,:),eld(:),eps(:),sigma(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:),xnew_pp(:)
  REAL(iwp),ALLOCATABLE :: u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),timest(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_tmp(:,:),eld_pp(:,:)!,tensor_pp(:,:,:)
  REAL(iwp),ALLOCATABLE :: valf(:),store_pp(:),prop(:,:)
  REAL(iwp),ALLOCATABLE :: fun(:),shape_integral_pp(:,:)
  REAL(iwp),ALLOCATABLE :: stress_integral_pp(:,:),stressnodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal_integral_pp(:,:),princinodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal(:),reacnodes_pp(:)  
  REAL(iwp),ALLOCATABLE :: ndscal_pp(:,:),tempres(:)
  INTEGER,  ALLOCATABLE :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)
  INTEGER,  ALLOCATABLE :: no(:),no_pp(:),no_pp_temp(:),sense(:),etype_pp(:)
    
  !Temperature variables
    
  REAL(iwp),ALLOCATABLE :: dtel(:),dtemp(:),etl(:),cte(:),teps(:)
  INTEGER(iwp),ALLOCATABLE :: num(:)
  REAL(iwp),ALLOCATABLE :: etl_pp(:,:)
  REAL(iwp),ALLOCATABLE :: tload_pp(:) 
 
!------------------------------------------------------------------------------
! 2a. Definition of variable names not listed in the 5th edition
!------------------------------------------------------------------------------

! Scalar reals

! cte       coefficient of thermal expansion in x, y, z direction


! Dynamic real arrays

! etl_pp    distributed element thermal forces array
! tload_pp  distributed global thermal forces vector
! dtemp     storage of all nodal temperature changes
! dtel      vectore storage for temperature nodal change for each element
! etl       element thermal force array for an individual element
! teps      

!------------------------------------------------------------------------------
! 3. Read job_name and instance_id from the command line. 
!    Read control data, mesh data, boundary and loading conditions. 
!------------------------------------------------------------------------------
  
  !-- To run in rfem_te mode set prog=11, for xx12_te prog=12   
  prog=12
  
  ALLOCATE(timest(20))
  timest    = zero 
  timest(1) = elap_time()
  
  CALL find_pe_procs(numpe,npes)
  argc = iargc()
  
  IF(prog==11)THEN
    IF( argc /= 2 ) THEN
      PRINT*
      PRINT*, "Usage:  xx12_te <model_name> <instance-id>"
      PRINT*
      PRINT*, "        program expects as input:"
      PRINT*, "          <model_name>-<instance-id>.d"
      PRINT*, "          <model_name>-<instance-id>.dat"
      PRINT*, "          <model_name>-<instance-id>.bnd"
      PRINT*, "          <model_name>-<instance-id>.fix"
      PRINT*, "          <model_name>-<instance-id>.mat"
      PRINT*
      PRINT*, "        and outputs:"
      PRINT*, "          <model_name>-<instance-id>.dis" 
      PRINT*, "          <model_name>-<instance-id>.pri" 
      PRINT*, "          <model_name>-<instance-id>.str" 
      PRINT*, "          <model_name>-<instance-id>.vms" 
      PRINT*, "          <model_name>-<instance-id>.rea" 
      PRINT*, "          <model_name>-<instance-id>.tnc" 
      PRINT*, "          <model_name>-<instance-id>.res" 
      PRINT*
      CALL job_name_error(numpe,program_name)
    END IF
    CALL GETARG(1, job_name) 
    CALL GETARG(2, instance_id)
    
    inst_in = job_name(1:LEN_TRIM(job_name)) // "-" // instance_id(1:LEN_TRIM(instance_id))
    CALL read_rfemsolve(inst_in,numpe,element,fixed_freedoms,limit,loaded_nodes, &
                meshgen,mises,nels,nip,nn,nod,np_types,nr,partitioner,tol)
    
    CALL calc_nels_pp(inst_in,nels,npes,numpe,partitioner,nels_pp)
  END IF
  IF(prog==12)THEN
    IF(argc /= 3) THEN
      IF(numpe==1)PRINT*
      IF(numpe==1)PRINT*, "Usage:  xx12_te <job_name> <stepnum> <val0>"
      IF(numpe==1)PRINT*
      CALL job_name_error(numpe,program_name)
    END IF
    CALL GETARG(1, job_name)
    CALL GETARG(2, stepnum)
    CALL GETARG(3, val0_str)
    READ(val0_str,*)val0
    IF(numpe==1)PRINT*,"Zero stress temperature = ",val0
    
    CALL read_rfemsolve(job_name,numpe,element,fixed_freedoms,limit,loaded_nodes, &
                meshgen,mises,nels,nip,nn,nod,np_types,nr,partitioner,tol)
    
!    inst_in = job_name(1:LEN_TRIM(job_name)) // "-" // instance_id(1:LEN_TRIM(instance_id))
!    fname = job_name(1:INDEX(job_name, " ")-1) // "-" // instance_id(1:LEN_TRIM(instance_id))
!    fname = job_name(1:INDEX(job_name, " ")-1)
    
    CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)
  END IF
  
  ndof = nod*nodof
  ntot = ndof
  
  ALLOCATE(g_num_pp(nod, nels_pp)) 
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp)) 
  ALLOCATE(rest(nr,nodof+1)) 
  ALLOCATE(etype_pp(nels_pp))
  ALLOCATE(prop(nprops,np_types))
  
  !Force Vector
  ALLOCATE(num(nod))
  ALLOCATE(etl(ndof))
  ALLOCATE(dtel(nod))
  ALLOCATE(dtemp(nn))
  ALLOCATE(teps(nst))
  
  ALLOCATE(etl_pp(ndof,nels_pp))
 
  g_num_pp   = 0
  rest       = 0
  etype_pp   = 0
  g_coord_pp = zero
  prop       = zero
  
  timest(2) = elap_time()
  
  IF(prog==11)THEN
    CALL read_elements(inst_in,iel_start,nn,npes,numpe,etype_pp,g_num_pp)
  END IF
  IF(prog==12)THEN
!    CALL read_elements(job_name,iel_start,nn,npes,numpe,etype_pp,g_num_pp)
    !---Open ascii ENSIGHT GOLD
    IF(numpe==1) PRINT *, "CALL read_g_num_pp_e"
    CALL read_g_num_pp_e(job_name,iel_start,nn,npes,numpe,g_num_pp)
    IF(numpe==1) PRINT *, "CALL read_etype_pp"
    CALL read_etype_pp(job_name,npes,numpe,etype_pp)
  END IF
  timest(3) = elap_time()
  IF(numpe==1) PRINT *, "READ_ELEMENTS COMPLETED"

! CALL read_g_num_pp(job_name,iel_start,nels,nn,numpe,g_num_pp)

  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  timest(4) = elap_time()
  
  IF(prog==11)THEN
    CALL read_g_coord_pp(inst_in,g_num_pp,nn,npes,numpe,g_coord_pp)
  END IF
  IF(prog==12)THEN
!    CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
    !---Open ascii ENSIGHT GOLD
    IF(numpe==1) PRINT *, "CALL read_g_coord_pp_e"
    CALL read_g_coord_pp_e(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  END IF
  timest(5) = elap_time()
  IF(numpe==1) PRINT *, "READ_G_COORD_PP COMPLETED"

  IF (nr>0) CALL read_rest(job_name,numpe,rest)
  timest(6) = elap_time()
  IF(numpe==1) PRINT *, "READ_REST COMPLETED"
  
  IF(prog==11) fname = inst_in(1:LEN_TRIM(inst_in)) // ".mat" ! Move to subroutine
  IF(prog==12) fname = job_name(1:INDEX(job_name, " ")-1) // ".mat"
  CALL read_materialValue(prop,fname,numpe,npes)
  
!------------------------------------------------------------------------------
! Nodal temperatures
!------------------------------------------------------------------------------

  !OPEN TEMPERATURE FILES
  !Read the temperature change at each node and assign to the dtemp vector
  
  ALLOCATE(cte(3))
  
  IF(prog==11)THEN
    OPEN(55,File='Temperature.txt',STATUS='OLD',ACTION='read')
    READ(55,*)keyword
    READ(55,*)
    
    DO i=1, nn
        READ(55,*)k,dtemp(i)
    END DO
    
    CLOSE(55) 
    
    bufsize = nn
    CALL MPI_BCAST(dtemp,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
  END IF
  
!!!!!!!!!!!!!!!!
!  fname = inst_in(1:INDEX(inst_in, " ")-1) // ".ensi.NDTTR-000056"
!  PRINT *, "fname = ",fname
!  IF(numpe==1)THEN
!    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
!!    OPEN(10,FILE='Tethra-elem-1.ensi.NDTTR-000056',STATUS='OLD',ACTION='READ')
!    READ(10,*)   cbuffer
!    READ(10,*)   cbuffer
!    READ(10,*)  int_in
!    READ(10,*)   cbuffer
!!    READ(10,'(A)')   cbuffer
!    DO i=1, nn
!      READ(10,*)  test
!      PRINT *, test
!    END DO
!    CLOSE(10)
!  END IF
!!!!!!!!!!!!!!!!
  
  IF(prog==12)THEN
    ALLOCATE(ndscal_pp(nod,nels_pp))
    ndscal_pp = zero
    CALL read_ensi_scalar_pn(job_name,g_num_pp,nn,npes,numpe,stepnum,ndscal_pp)
  
!    PRINT *, "ndscal_pp = "
!    PRINT *, ndscal_pp
  END IF
  
  PRINT *, "READ_MATERIALVALUE COMPLETED"
  
!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),principal(ndim),      &
           der(ndim,nod),deriv(ndim,nod),bee(nst,ntot),                       &
           storkm_pp(ntot,ntot,nels_pp),eld(ntot),eps(nst),sigma(nst),        &
           pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),                      &
           weights(nip),g_g_pp(ntot,nels_pp),fun(nod))

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------
  
  CALL rearrange(rest)
!  IF(prog==11) CALL rearrange(rest)
!  IF(prog==12 .AND. nr>0) CALL rearrange_2(rest)
  IF(numpe==1) PRINT *, "REARRANGE COMPLETED"
  
  g_g_pp = 0

!  IF(prog==11)THEN
    elements_1: DO iel = 1, nels_pp
!     CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
      CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
    END DO elements_1
!  END IF
!  IF(prog==12)THEN
!    ! When nr = 0, g_num_pp and g_g_pp are identical
!    IF(nr>0) THEN
!      elements_1b: DO iel = 1, nels_pp
!        CALL find_g4(g_num_pp(:,iel),g_g_pp(:,iel),rest)
!      END DO elements_1b
!      DEALLOCATE(rest)
!    ELSE
!      g_g_pp = g_num_pp
!    END IF
!  END IF
  
  !PRINT*,'g_num_pp'
  !PRINT*,g_num_pp

  IF(numpe==1) PRINT *, "FIND_G COMPLETED"

  neq = 0
  
  elements_2: DO iel = 1, nels_pp  
    i = MAXVAL(g_g_pp(:,iel))
    IF(i > neq) neq = i
  END DO elements_2  

  neq = MAX_INTEGER_P(neq)
  
  timest(7) = elap_time()

!------------------------------------------------------------------------------
! 6. Create interprocessor communication tables
!------------------------------------------------------------------------------

  CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)
 
  timest(8) = elap_time()

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  IF(numpe==1) PRINT *, "CREATED INTERPROCESSOR COMMUNICATION TABLES"

!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------
  
  ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),            &
           u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp))

  p_pp    = zero  ;  r_pp = zero  ;  x_pp = zero
  xnew_pp = zero  ;  u_pp = zero  ;  d_pp = zero  ; diag_precon_pp = zero

  !CHANGE array size
  
  ALLOCATE(tload_pp(neq_pp))
 
  
  timest(9) = elap_time()

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------

  points = zero
  
!  IF(numpe==1) PRINT *, "Reached 'CALL sample'"
  CALL sample(element,points,weights)

  teps      = zero
  storkm_pp = zero
  etl_pp    = zero
  cte       = zero

! The relationship between CTE and Young's modulus was removed for the Hartree
! summer school 2016. For more information of this code please search for the 
! journal paper: Spatial variability in the coefficient of thermal expansion 
! induces pre-service stresses in computer models of virgin Gilsocarbon bricks
! doi:10.1016/j.jnucmat.2015.05.058

!  constant  = 0.0435
!  IF(numpe==1) PRINT *, "Reached elements_3 loop"
  elements_3: DO iel=1,nels_pp
    
    cte (1)   = prop(8,etype_pp(iel))
    cte (2)   = prop(8,etype_pp(iel))
    cte (3)   = prop(8,etype_pp(iel))
    dee = zero
    
!Relationship between CTE and Young's modulus
    !e = constant/prop(1,etype_pp(iel))
    e = prop(6,etype_pp(iel))
    v = prop(7,etype_pp(iel))
    !e = 10000
    
    CALL deemat(dee,e,v)
    
    !Extraction of nodal temperature changes for each element
    IF(prog==11)THEN
      !Extraction of nodal numbers for each element
      num=g_num_pp(:,iel)
      dtel=dtemp(num)
    END IF
    IF(prog==12)THEN
      dtel=ndscal_pp(:,iel)-val0
    END IF
    
    etl=zero
    
    gauss_pts_1: DO i=1,nip
               
      CALL shape_fun(fun,points,i)
              
      !Temperature calculation at each integration point
      gtemp=dot_product(fun,dtel)
            
      !Calculation of the force vector
      ! ε = αΔT{1 1 1 0 0 0}
      teps(1:3)=gtemp*cte(1:3)
                        
      CALL shape_der(der,points,i)
      jac   = MATMUL(der,g_coord_pp(:,:,iel))
      det   = determinant(jac)
          
      CALL invert(jac)
      deriv = MATMUL(jac,der)
         
      CALL beemat(bee,deriv)
      
      storkm_pp(:,:,iel)   = storkm_pp(:,:,iel) +                             &
                             MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) *         &
                             det*weights(i)
    
      !Calculation of thermal gradient force vector - [B]T[D]{ε} 
      etl=etl+MATMUL(MATMUL(TRANSPOSE(bee),dee),teps)*det*weights(i)
          
    END DO gauss_pts_1
     
    !Storage and distribution of thermal gradient force vector
    
    etl_pp(:,iel)=etl_pp(:,iel)+etl

    
  END DO elements_3
  
  timest(10) = elap_time()
  
  IF(numpe==1) PRINT*, "COMPLETED ELEMENT STIFFNESS INTEGRATION AND STORAGE"
    
  CALL scatter(tload_pp,etl_pp)
  

!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------
  
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  diag_precon_tmp = zero
 
  elements_4: DO iel = 1,nels_pp 
    DO i = 1,ndof
      diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + storkm_pp(i,i,iel)
    END DO
  END DO elements_4
  
  CALL scatter(diag_precon_pp,diag_precon_tmp)

  DEALLOCATE(diag_precon_tmp)

  timest(11) = elap_time()

  IF(numpe==1) PRINT *, "BUILT THE DIAGONAL PRECONDITIONER"

!------------------------------------------------------------------------------
! 10. Read in fixed nodal displacements and assign to equations
!------------------------------------------------------------------------------

  IF(fixed_freedoms > 0) THEN

    ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),valf(fixed_freedoms),    &
             no_pp_temp(fixed_freedoms),sense(fixed_freedoms))
    
    node = 0 ; no = 0 ; no_pp_temp = 0 ; sense = 0 ; valf = zero

    CALL read_fixed(job_name,numpe,node,sense,valf)
    CALL find_no(node,rest,sense,no)
    CALL reindex(ieq_start,no,no_pp_temp,fixed_freedoms_pp,                   &
                             fixed_freedoms_start,neq_pp)

    ALLOCATE(no_pp(fixed_freedoms_pp),store_pp(fixed_freedoms_pp))

    no_pp    = 0
    store_pp = zero
    no_pp    = no_pp_temp(1:fixed_freedoms_pp)

    DEALLOCATE(node,no,sense,no_pp_temp,rest)

  END IF

  IF(fixed_freedoms == 0) fixed_freedoms_pp = 0

!  DEALLOCATE(rest)

  IF(numpe==1) PRINT *, "READ FIXED NODAL DISPLACEMENTS"

!------------------------------------------------------------------------------
! 11. Read in loaded nodes and get starting r_pp
!------------------------------------------------------------------------------
 
  IF(loaded_nodes > 0) THEN

    ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes))
    
    val  = zero ; node = 0

    CALL read_loads(job_name,numpe,node,val)
    CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))

    !CHECK
    tload = SUM_P(r_pp(1:))

    DEALLOCATE(node,val)

  ELSE

    tload = zero

  END IF
  
  DEALLOCATE(g_g_pp)
  
  timest(12) = elap_time()

!------------------------------------------------------------------------------
! 12. Invert the preconditioner.
!     If there are fixed freedoms, first apply a penalty
!------------------------------------------------------------------------------

  IF(fixed_freedoms_pp > 0) THEN
    DO i = 1, fixed_freedoms_pp
       j                 = no_pp(i) - ieq_start + 1
       diag_precon_pp(j) = diag_precon_pp(j) + penalty
       store_pp(i)       = diag_precon_pp(j)
    END DO
  END IF

  diag_precon_pp = 1._iwp/diag_precon_pp
 
!------------------------------------------------------------------------------
! 13. Initiallize preconditioned conjugate gradient
!------------------------------------------------------------------------------

  IF(fixed_freedoms_pp > 0) THEN
    DO i = 1,fixed_freedoms_pp
       j       = no_pp(i) - ieq_start + 1
       k       = fixed_freedoms_start + i - 1
       r_pp(j) = store_pp(i) * valf(k)
    END DO
  END IF
  
  !CHANGE 
  
  !AQUI ES
  r_pp = tload_pp - r_pp
  
  d_pp  = diag_precon_pp*r_pp
  p_pp  = d_pp
  x_pp  = zero

!------------------------------------------------------------------------------
! 14. Preconditioned conjugate gradient iterations
!------------------------------------------------------------------------------

  iters = 0

  iterations: DO 
    iters    = iters + 1
    u_pp     = zero
    pmul_pp  = zero
    utemp_pp = zero
    
    CALL gather(p_pp,pmul_pp)
    
    elements_5: DO iel=1,nels_pp
      utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
    END DO elements_5
    
   CALL scatter(u_pp,utemp_pp)

    IF(fixed_freedoms_pp > 0) THEN
      DO i = 1, fixed_freedoms_pp
        j       = no_pp(i) - ieq_start + 1
        u_pp(j) = p_pp(j) * store_pp(i)
      END DO
    END IF

    up      = DOT_PRODUCT_P(r_pp,d_pp)
    alpha   = up/DOT_PRODUCT_P(p_pp,u_pp)
    xnew_pp = x_pp + p_pp*alpha
    r_pp    = r_pp - u_pp*alpha
    d_pp    = diag_precon_pp*r_pp
    beta    = DOT_PRODUCT_P(r_pp,d_pp)/up
    p_pp    = d_pp + p_pp*beta  

    !IF(numpe==1) PRINT *, "ITERATION", iters

    CALL checon_par(xnew_pp,tol,converged,x_pp)    
    IF(converged.OR.iters==limit)EXIT

  END DO iterations

  DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,storkm_pp,pmul_pp) 
  
  timest(13) = elap_time()

!------------------------------------------------------------------------------
! 15. Print out displacements, stress, principal stress and reactions
!------------------------------------------------------------------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  
  IF(prog==11)THEN
    IF(numpe==1) THEN
      fname = inst_in(1:LEN_TRIM(inst_in)) // ".dis"
      OPEN(24, file=fname, status='replace', action='write')
      fname = inst_in(1:LEN_TRIM(inst_in)) // ".str"
      OPEN(25, file=fname, status='replace', action='write')
      fname = inst_in(1:LEN_TRIM(inst_in)) // ".pri"
      OPEN(26, file=fname, status='replace', action='write')
      fname = inst_in(1:LEN_TRIM(inst_in)) // ".vms"
      OPEN(27, file=fname, status='replace', action='write')
      fname = inst_in(1:LEN_TRIM(inst_in)) // ".rea"
      OPEN(28, file=fname, status='replace', action='write')
    END IF
  END IF
  IF(prog==12)THEN
    IF(numpe==1) THEN
      fname=job_name(1:INDEX(job_name," ")-1)//".ensi.DIS-"//stepnum
      OPEN(24, file=fname, status='replace', action='write')
      fname=job_name(1:INDEX(job_name," ")-1)//".ensi.STR-"//stepnum
!      fname = job_name(1:LEN_TRIM(job_name)) // ".str"
      OPEN(25, file=fname, status='replace', action='write')
      fname=job_name(1:INDEX(job_name," ")-1)//".ensi.PRI-"//stepnum
!      fname = job_name(1:LEN_TRIM(job_name)) // ".pri"
      OPEN(26, file=fname, status='replace', action='write')
      fname=job_name(1:INDEX(job_name," ")-1)//".ensi.VMS-"//stepnum
!      fname = job_name(1:LEN_TRIM(job_name)) // ".vms"
      OPEN(27, file=fname, status='replace', action='write')
      fname=job_name(1:INDEX(job_name," ")-1)//".ensi.REA-"//stepnum
!      fname = job_name(1:LEN_TRIM(job_name)) // ".rea"
      OPEN(28, file=fname, status='replace', action='write')
    END IF
  END IF
  
!------------------------------------------------------------------------------
! 16a. Displacements
!------------------------------------------------------------------------------

  ALLOCATE(eld_pp(ntot,nels_pp)) 
  eld_pp = zero
  CALL gather(xnew_pp(1:),eld_pp)
  DEALLOCATE(xnew_pp)

  ALLOCATE(disp_pp(nodes_pp*ndim))
  disp_pp = zero

  label   = "*DISPLACEMENT"

  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,              &
                     node_start,node_end,eld_pp,disp_pp,1)
  IF(prog==11)THEN
    CALL write_nodal_variable(label,24,1,nodes_pp,npes,numpe,ndim,disp_pp)
    IF(numpe==1) CLOSE(24)
  END IF
  IF(prog==12)THEN
    !---Write file: outputs-ENSIGHT GOLD
    IF(numpe==1)THEN
      WRITE(24,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
      WRITE(24,'(A)') "part"
      WRITE(24,'(I1)') 1
      WRITE(24,'(A)') "coordinates"
    END IF
    
    ALLOCATE(tempres(nodes_pp))
!    tempres = zero
    DO i=1,ndim
      tempres = zero
      DO j=1,nodes_pp
        !tempres(nodes_pp*(i-1)+j)=disp_pp(ndim*(j-1)+i)
        tempres(j)=disp_pp(ndim*(j-1)+i)
      END DO
      CALL dismsh_ensi_p(24,1,nodes_pp,npes,numpe,nodof,tempres)
    END DO
    DEALLOCATE(tempres)
    IF(numpe==1) CLOSE (24)
  END IF
  
  DEALLOCATE(disp_pp)
  
  IF(numpe==1) PRINT *, "OUTPUT DISPLACEMENTS"

!------------------------------------------------------------------------------
! 16b. Stresses
!------------------------------------------------------------------------------

  ALLOCATE(shape_integral_pp(nod,nels_pp))
  ALLOCATE(stress_integral_pp(nod*nst,nels_pp))
  ALLOCATE(stressnodes_pp(nodes_pp*nst))
  ALLOCATE(principal_integral_pp(nod*nodof,nels_pp))
  ALLOCATE(princinodes_pp(nodes_pp*nodof))
  ALLOCATE(reacnodes_pp(nodes_pp*nodof))
  
  shape_integral_pp     = zero
  stress_integral_pp    = zero
  stressnodes_pp        = zero
  principal_integral_pp = zero  
  princinodes_pp        = zero
  reacnodes_pp          = zero
  utemp_pp              = zero

  v = 0
  
! The relationship between CTE and Young's modulus was removed for the Hartree
! summer school 2016. For more information of this code please search for the 
! journal paper: Spatial variability in the coefficient of thermal expansion 
! induces pre-service stresses in computer models of virgin Gilsocarbon bricks
! doi:10.1016/j.jnucmat.2015.05.058


  DO iel = 1,nels_pp
              
!    cte (1)   = prop(1,etype_pp(iel))
!    cte (2)   = prop(1,etype_pp(iel))
!    cte (3)   = prop(1,etype_pp(iel))
!    v   = prop(2,etype_pp(iel))
!    !e = constant/prop(1,etype_pp(iel))
!    e = 10000
    
    cte (1) = prop(8,etype_pp(iel))
    cte (2) = prop(8,etype_pp(iel))
    cte (3) = prop(8,etype_pp(iel))
    e       = prop(6,etype_pp(iel))
    v       = prop(7,etype_pp(iel))
    
    dee = zero
    CALL deemat(dee,e,v)

    DO i = 1,nip
      CALL shape_der(der,points,i)
      jac   = MATMUL(der,g_coord_pp(:,:,iel))
      det   = DETERMINANT(jac) 
      CALL invert(jac)
      deriv = MATMUL(jac,der)
      CALL beemat(bee,deriv)
      eps   = MATMUL(bee,eld_pp(:,iel))
      sigma = MATMUL(dee,eps)
      CALL PRINCIPALSTRESS3D(sigma,principal)
      utemp_pp(:,iel) = utemp_pp(:,iel) +                                    &
                        MATMUL(TRANSPOSE(bee),sigma)*det*weights(i)

      CALL shape_fun(fun,points,i)

      DO j = 1,nod
        idx1 = (j-1)*nst
        idx2 = (j-1)*nodof
        shape_integral_pp(j,iel) = shape_integral_pp(j,iel) +                 &
                                   fun(j)*det*weights(i)
        DO k = 1,nst
          stress_integral_pp(idx1+k,iel) = stress_integral_pp(idx1+k,iel) +   &
                                           fun(j)*sigma(k)*det*weights(i)
        END DO
        DO k = 1,nodof
          principal_integral_pp(idx2+k,iel) = principal_integral_pp(idx2+k,iel) + &
                                              fun(j)*principal(k)*det*weights(i)
        END DO
      END DO
    END DO !gauss
  END DO !elements
    
!------------------------------------------------------------------------------
! 16c. Stress
!------------------------------------------------------------------------------

  label = "*STRESS"

  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,            &
                        node_start,node_end,shape_integral_pp,                &
                        stress_integral_pp,stressnodes_pp)
!  CALL write_nodal_variable(label,25,1,nodes_pp,npes,numpe,nst,               &
!                            stressnodes_pp)
                            
  IF(prog==11)THEN
    CALL write_nodal_variable(label,25,1,nodes_pp,npes,numpe,nst,               &
                              stressnodes_pp)
    IF(numpe==1) CLOSE(25)
  END IF
  IF(prog==12)THEN
    !---Write file: outputs-ENSIGHT GOLD
    IF(numpe==1)THEN
      WRITE(25,'(A)') "Alya Ensight Gold --- Tensor per-node variable file"
      WRITE(25,'(A)') "part"
      WRITE(25,'(I1)') 1
      WRITE(25,'(A)') "coordinates"
    END IF
    
    ALLOCATE(tempres(nodes_pp))
    DO i=1,nst
      tempres = zero
      DO j=1,nodes_pp
        tempres(j)=stressnodes_pp(nst*(j-1)+i)
      END DO
      CALL dismsh_ensi_p(25,1,nodes_pp,npes,numpe,nodof,tempres)
    END DO
    DEALLOCATE(tempres)
    IF(numpe==1) CLOSE (25)
  END IF
  
  DEALLOCATE(stress_integral_pp,stressnodes_pp)
  
  IF(numpe==1) PRINT *, "OUTPUT STRESS TENSOR"
  
!------------------------------------------------------------------------------
! 16d. Principal stress
!------------------------------------------------------------------------------
  
  label = "*PRINCIPAL STRESS"
  
  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,          &
                        node_start,node_end,shape_integral_pp,                &
                        principal_integral_pp,princinodes_pp)
!  CALL write_nodal_variable(label,26,1,nodes_pp,npes,numpe,nodof,             &
!                            princinodes_pp)
  
  IF(prog==11)THEN
    CALL write_nodal_variable(label,26,1,nodes_pp,npes,numpe,nodof,             &
                              princinodes_pp)
    IF(numpe==1) CLOSE(26)
  END IF
  IF(prog==12)THEN
    !---Write file: outputs-ENSIGHT GOLD
    IF(numpe==1)THEN
      WRITE(26,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
      WRITE(26,'(A)') "part"
      WRITE(26,'(I1)') 1
      WRITE(26,'(A)') "coordinates"
    END IF
    
    ALLOCATE(tempres(nodes_pp))
    DO i=1,ndim
      tempres = zero
      DO j=1,nodes_pp
        tempres(j)=princinodes_pp(ndim*(j-1)+i)
      END DO
      CALL dismsh_ensi_p(26,1,nodes_pp,npes,numpe,nodof,tempres)
    END DO
    DEALLOCATE(tempres)
    IF(numpe==1) CLOSE (26)
  END IF
  
  DEALLOCATE(principal_integral_pp)
  
  IF(numpe==1) PRINT *, "OUTPUT PRINCIPAL STRESS VECTOR"
  
!------------------------------------------------------------------------------
! 16e. Von Mises stress (rho_v)
!      rho_v = sqrt( ( (rho1-rho2)^2 + (rho2-rho3)^2 + (rho1-rho3)^2 ) / 2 )
!------------------------------------------------------------------------------
  
  label = "*MISES STRESS"
  
  DO i = 1,nodes_pp
    j = ((i-1)*nodof)+1
    k = j + 1
    l = j + 2
    princinodes_pp(j) = SQRT(((princinodes_pp(j)-princinodes_pp(k)) **2 +     &
                              (princinodes_pp(k)-princinodes_pp(l)) **2 +     &
                              (princinodes_pp(j)-princinodes_pp(l)) **2)      &
                              * 0.5_iwp)
    princinodes_pp(k:l) = zero
  END DO

  IF(prog==11)THEN
    CALL write_nodal_variable(label,27,1,nodes_pp,npes,numpe,nodof,             &
                              princinodes_pp)
    IF(numpe==1) CLOSE(27)
  END IF
  IF(prog==12)THEN
    !---Write file: outputs-ENSIGHT GOLD
    IF(numpe==1)THEN
      WRITE(27,'(A)') "Alya Ensight Gold --- Scalar per-node variable file"
      WRITE(27,'(A)') "part"
      WRITE(27,'(I1)') 1
      WRITE(27,'(A)') "coordinates"
    END IF
    
    ALLOCATE(tempres(nodes_pp))
!    DO i=1,ndim
      tempres = zero
      DO j=1,nodes_pp
        tempres(j)=princinodes_pp(ndim*(j-1)+1)
!        tempres(j)=princinodes_pp(j)
      END DO
      CALL dismsh_ensi_p(27,1,nodes_pp,npes,numpe,nodof,tempres)
!    END DO
    DEALLOCATE(tempres)
    IF(numpe==1) CLOSE (27)
  END IF
  
  ! count number of nodes with value above threshold

  nodecount_pp(1) = 0 ; nodecount = 0
 
  DO i = 1,nodes_pp
    IF(princinodes_pp(((i-1)*3)+1) > mises) nodecount_pp(1) = nodecount_pp(1) + 1
  END DO

  nodecount = ISUM_P(nodecount_pp) 
                              
  DEALLOCATE(princinodes_pp)

  IF(numpe==1) PRINT *, "OUTPUT VON MISES SCALAR"

!------------------------------------------------------------------------------
! 16f. Reactions
!------------------------------------------------------------------------------

  label = "*NODAL REACTIONS"
  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,             &
                     node_start,node_end,utemp_pp,reacnodes_pp,0)
  
  IF(prog==11)THEN
    CALL write_nodal_variable(label,28,1,nodes_pp,npes,numpe,nodof,           &
                              reacnodes_pp)
    IF(numpe==1) CLOSE(28)
  END IF
  IF(prog==12)THEN
    !---Write file: outputs-ENSIGHT GOLD
    IF(numpe==1)THEN
      WRITE(28,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
      WRITE(28,'(A)') "part"
      WRITE(28,'(I1)') 1
      WRITE(28,'(A)') "coordinates"
    END IF
    
    ALLOCATE(tempres(nodes_pp))
    DO i=1,ndim
      tempres = zero
      DO j=1,nodes_pp
        tempres(j)=reacnodes_pp(ndim*(j-1)+i)
      END DO
      CALL dismsh_ensi_p(28,1,nodes_pp,npes,numpe,nodof,tempres)
    END DO
    DEALLOCATE(tempres)
    IF(numpe==1) CLOSE (28)
  END IF
  
  DEALLOCATE(reacnodes_pp,shape_integral_pp)
  
  timest(14) = elap_time()
  
  IF(numpe==1) PRINT *, "OUTPUT NODAL REACTIONS VECTOR"
  
!------------------------------------------------------------------------------
! 16g. ENSIGHT GOLD Case file
!------------------------------------------------------------------------------
  
  IF(prog==12)THEN
    IF(numpe==1) THEN
      fname=job_name(1:INDEX(job_name," ")-1)// "-te-" //stepnum(1:INDEX(stepnum, " ")-1)//".ensi.case"
      OPEN(29, file=fname, status='replace', action='write')
      
      WRITE(29,'(A)') "#"
      WRITE(29,'(A)') "# ParaFEM generated post-process file"
      WRITE(29,'(A)') "# Ensight Gold Format"
      WRITE(29,'(A)') "#"
      WRITE(29,'(2A/A)') "# Problem name: ",job_name(1:INDEX(job_name, " ")-1),"#"
      WRITE(29,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
      WRITE(29,'(2A/A)')   "model:  ",job_name(1:INDEX(job_name, " ")-1)//'.ensi.geo',"VARIABLE"
      WRITE(29,'(2A)')     "scalar per element:   MaterialID       ",            &
                          job_name(1:INDEX(job_name, " ")-1)//'.ensi.MATID'
      WRITE(29,'(3A)')     "scalar per node:      Temperature      ",            &
                          job_name(1:INDEX(job_name, " ")-1)//'.ensi.NDTTR-',stepnum(1:INDEX(stepnum, " ")-1)
      WRITE(29,'(3A)')     "vector per node:      Displacements    ",            &
                          job_name(1:INDEX(job_name, " ")-1)//'.ensi.DIS-',stepnum(1:INDEX(stepnum, " ")-1)
      WRITE(29,'(3A)')     "tensor symm per node: Stress           ",            &
                          job_name(1:INDEX(job_name, " ")-1)//'.ensi.STR-',stepnum(1:INDEX(stepnum, " ")-1)
      WRITE(29,'(3A)')     "vector per node:      PrincipalStress  ",            &
                          job_name(1:INDEX(job_name, " ")-1)//'.ensi.PRI-',stepnum(1:INDEX(stepnum, " ")-1)
      WRITE(29,'(3A)')     "scalar per node:      VonMises         ",            &
                          job_name(1:INDEX(job_name, " ")-1)//'.ensi.VMS-',stepnum(1:INDEX(stepnum, " ")-1)
      WRITE(29,'(3A)')     "vector per node:      NodalReactions   ",            &
                          job_name(1:INDEX(job_name, " ")-1)//'.ensi.REA-',stepnum(1:INDEX(stepnum, " ")-1)
      
      CLOSE (29)
      
      PRINT *, "OUTPUT ENSIGHT GOLD CASE FILE"
    END IF
  END IF
  
!------------------------------------------------------------------------------
! 17. Output performance data
!------------------------------------------------------------------------------
  
  IF(prog==11)THEN
    CALL WRITE_RFEMSOLVE(fixed_freedoms,iters,inst_in,loaded_nodes,mises,neq,   &
                         nn,nodecount,npes,nr,numpe,timest,tload)
  END IF
  IF(prog==12)THEN
    CALL WRITE_RFEMSOLVE(fixed_freedoms,iters,job_name,loaded_nodes,mises,neq,   &
                         nn,nodecount,npes,nr,numpe,timest,tload)
  END IF
 
  CALL shutdown() 
 
END PROGRAM xx12_te
