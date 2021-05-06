PROGRAM XX9       
!------------------------------------------------------------------------------ 
!      Program XX.9 Three dimensional analysis of an elastic solid using load
!                   control or displacement control. GPU support using CUDA.
!                   Special case with single element stiffness matrix and 
!                   matrix-matrix kernal in PCG.
!                   See program P121.
!------------------------------------------------------------------------------ 
                                 
  USE precision     ; USE global_variables ; USE mp_interface
  USE input         ; USE output           ; USE loading
  USE timing        ; USE maths            ; USE gather_scatter
  USE partition     ; USE elements         ; USE steering        ; USE pcg
  USE new_library   ; USE iso_c_binding    ;

  IMPLICIT NONE

!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

 ! neq,ntot are now global variables - not declared

  INTEGER,PARAMETER     :: nodof=3,ndim=3,nst=6
  INTEGER               :: loaded_nodes,fixed_freedoms,iel,i,j,k,l,idx1,idx2
  INTEGER               :: iters,limit,nn,nr,nip,nod,nels,ndof,npes_pp
  INTEGER               :: node_end,node_start,nodes_pp
  INTEGER               :: argc,iargc,meshgen,partitioner
  INTEGER               :: fixed_freedoms_pp,fixed_freedoms_start
  REAL(iwp)             :: e,v,det,tol,up,up0,up1,alpha,beta,tload
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp
  REAL(iwp),PARAMETER   :: penalty = 1.0e20_iwp
  CHARACTER(LEN=3)      :: xpu
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: program_name='xx9'
  CHARACTER(LEN=50)     :: fname,job_name,label
  CHARACTER(LEN=80)     :: cbuffer
  LOGICAL               :: converged = .false.
  LOGICAL               :: output_stress = .false.

  ! GPU related variables etc
  ! -------------------------
  integer :: cublas_alloc, cublas_free
  integer :: cublas_set_matrix, cublas_get_matrix
  integer :: cublas_init, cublas_shutdown

  integer :: status
  integer :: ndof_per_element
  integer(kind=c_size_t) :: device_matrix
  integer(kind=c_size_t) :: device_lhs_vectors
  integer(kind=c_size_t) :: device_rhs_vectors

  !RZ Pointers for  p_pp and u_pp to move them onto the GPU
  integer(kind=c_size_t) :: device_p_pp
  integer(kind=c_size_t) :: device_u_pp
 
  logical :: use_gpu = .false.
  character(len=20, kind=c_char) :: op

  real(iwp) :: t_rawcomp
  real(iwp) :: t_comp
  real(iwp) :: t_start1
  real(iwp) :: t_start2

  !RZ Getting cublas dot product and set vector functions
  real(iwp) :: cublas_Ddot
  integer :: cublas_set_vector
  
  !Timer function
  integer :: m, n, ans
  real :: startT, endT, execTime
  

  ! Interface to function to set device
  interface
     integer(c_int) function set_gpu(device_id) bind(C)
       
       use iso_c_binding
       
       integer(c_int) :: device_id
     end function set_gpu
     
     integer(c_int) function sync_gpu() bind(C)
       use iso_c_binding
     end function sync_gpu

  end interface
 
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:),val(:,:),disp_pp(:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:)
  REAL(iwp),ALLOCATABLE :: bee(:,:),km(:,:),eld(:),eps(:),sigma(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:),xnew_pp(:)
  REAL(iwp),ALLOCATABLE :: u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),timest(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_tmp(:,:),eld_pp(:,:),tensor_pp(:,:,:)
  REAL(iwp),ALLOCATABLE :: valf(:),store_pp(:)
  REAL(iwp),ALLOCATABLE :: fun(:),shape_integral_pp(:,:)
  REAL(iwp),ALLOCATABLE :: stress_integral_pp(:,:),stressnodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal_integral_pp(:,:),princinodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal(:),reacnodes_pp(:),tempres(:)  
  INTEGER,  ALLOCATABLE :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)
  INTEGER,  ALLOCATABLE :: no(:),no_pp(:),no_pp_temp(:),sense(:)

!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions. 
!------------------------------------------------------------------------------
 
  ALLOCATE(timest(20))
  timest    = zero 
  timest(1) = elap_time()
  
  CALL find_pe_procs(numpe,npes)
  argc = iargc()
  IF (argc /= 1) CALL job_name_error(numpe,program_name)
  CALL GETARG(1, job_name) 

  CALL read_xx3(job_name,numpe,e,element,limit,loaded_nodes,                  &
                 meshgen,nels,nip,nn,nod,nr,partitioner,tol,v,xpu)

  fixed_freedoms = 0

  IF(xpu=='gpu') use_gpu = .true.

  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)

  ndof = nod*nodof
  ntot = ndof
 
  ALLOCATE(g_num_pp(nod, nels_pp)) 
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp)) 
  ALLOCATE(rest(nr,nodof+1)) 
 
  g_num_pp  = 0
  g_coord_pp= zero
  rest      = 0

  timest(2) = elap_time()
  
  CALL read_g_num_pp(job_name,iel_start,nn,npes,numpe,g_num_pp)
  timest(3) = elap_time()

  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  timest(4) = elap_time()

  CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  timest(5) = elap_time()

  CALL read_rest(job_name,numpe,rest)
  timest(6) = elap_time()
    
!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),principal(ndim),      &
           der(ndim,nod),deriv(ndim,nod),bee(nst,ntot),                       &
           km(ntot,ntot),eld(ntot),eps(nst),sigma(nst),                       &
           pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),                      &
           weights(nip),g_g_pp(ntot,nels_pp),fun(nod))

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------

  CALL rearrange(rest)
  
  g_g_pp = 0

  elements_1: DO iel = 1, nels_pp
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
!   CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
  END DO elements_1

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
!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------
  
  ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),            &
           u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp))

  p_pp    = zero  ;  r_pp = zero  ;  x_pp = zero
  xnew_pp = zero  ;  u_pp = zero  ;  d_pp = zero  ; diag_precon_pp = zero

  timest(9) = elap_time()

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!    Special case: Assume all elements have the same size and material 
!    properties and therefore the same stiffness matrix
!    Assume first element on each processor has same stiffness matrix
!    This is not strictly true if we take into account numerical precision.
!------------------------------------------------------------------------------

  dee = zero
  km  = zero
  iel = 1

  CALL deemat(dee,e,v)
  CALL sample(element,points,weights)

  gauss_pts_1: DO i=1,nip
    CALL shape_der(der,points,i)
    jac   = MATMUL(der,g_coord_pp(:,:,iel))
    det   = determinant(jac)
    CALL invert(jac)
    deriv = MATMUL(jac,der)
    CALL beemat(bee,deriv)
    km    = km + MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) * det * weights(i)
  END DO gauss_pts_1
  
  timest(10) = elap_time()

!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------
  
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  diag_precon_tmp = zero
 
  elements_4: DO iel = 1,nels_pp 
    DO i = 1,ndof
      diag_precon_tmp (i,iel) = diag_precon_tmp(i,iel) + km(i,i)
    END DO
  END DO elements_4

  CALL scatter(diag_precon_pp,diag_precon_tmp)

  DEALLOCATE(diag_precon_tmp)

  timest(11) = elap_time()

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

    DEALLOCATE(node,no,sense,no_pp_temp)

  END IF

  IF(fixed_freedoms == 0) fixed_freedoms_pp = 0

  DEALLOCATE(rest)

!------------------------------------------------------------------------------
! 11. Read in loaded nodes and get starting r_pp
!------------------------------------------------------------------------------
 
  IF(loaded_nodes > 0) THEN

    ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes))
    
    val  = zero ; node = 0

    CALL read_loads(job_name,numpe,node,val)
    CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))

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

  d_pp  = diag_precon_pp*r_pp
  p_pp  = d_pp
  x_pp  = zero
  up0   = DOT_PRODUCT_P(r_pp,d_pp)

  ! Code to set up the gpu 
  if (use_gpu) then
     
     ! Execute this loop iff gpu not in exclusive mode
     if (.false.) then
        status = set_gpu(numpe-1)
        if (status > 0) then
           print *, "gpu memory failed to allocate!"
           stop
        end if
     end if

     ! Initialize CUBLAS
     status = cublas_init
     if (status .ne. 0) then
        print *, "GPU failed to initialise!"
        status = cublas_shutdown
        stop
     end if

     ndof_per_element = 3 * nod

     ! Allocate memory on the gpu for coefficient matrix
     status = cublas_alloc( &
          ndof_per_element**2, &
          sizeof(0.0d0), &
          device_matrix)
     if (status .ne. 0) then
        print *, "GPU memory failed to allocate!"
        status = cublas_shutdown
        stop
     end if

     ! Allocate memory for matrix of lhs vectors
     status = cublas_alloc( &
          nels_pp*ndof_per_element, &
          sizeof(0.0d0), &
          device_lhs_vectors)
     if (status .ne. 0) then
        print *, "GPU memory failed to allocate!"
        status = cublas_shutdown
        stop
     end if

     ! Allocate memory for matrix of rhs vectors
     status =  cublas_alloc( &
          nels_pp*ndof_per_element, &
          sizeof(0.0d0), &
          device_rhs_vectors)
     if (status .ne. 0) then
        print *, "GPU memory failed to allocate!"
        status = cublas_shutdown
        stop
     end if
     
     !RZ Allocate memory for p_pp and u_pp on GPU
     
     status = cublas_alloc( &
          nels_pp*ndof_per_element, &
          sizeof(0.0d0), &
          device_p_pp)
     if (status .ne. 0) then 
        print *,"GPU memory failed to allocate p_pp!"
        status = cublas_shutdown
        stop
     end if   


     status = cublas_alloc( &
          nels_pp*ndof_per_element, &
          sizeof(0.0d0), &
          device_u_pp)
     if (status .ne. 0) then
        print *, "GPU memory failed to allocate u_pp!"
        status = cublas_shutdown
        stop
     end if

     ! Copy coefficient matrix to the gpu
     status = cublas_set_matrix( &
          ndof_per_element, &
          ndof_per_element, &
          sizeof(0.0d0), &
          km, &
          ndof_per_element, &
          device_matrix, &
          ndof_per_element)
     if (status .ne. 0) then
        print *, "Failed to copy data to gpu!"
        status = cublas_shutdown
        stop
     end if

  end if

!------------------------------------------------------------------------------
! 14. Preconditioned conjugate gradient iterations
!------------------------------------------------------------------------------

  iters = 0
  t_rawcomp = 0
  t_comp = 0

  iterations: DO 
    iters    = iters + 1
    u_pp     = zero
    pmul_pp  = zero
    utemp_pp = zero
    u_pp     = zero

   CALL gather(p_pp,pmul_pp)

    t_start1 = elap_time()

    ! matmul version
    if (.not. use_gpu) then
   
       utemp_pp = MATMUL(km,pmul_pp)

       ! gpu version
    else
       ! Copy lhs vectors to gpu
       status = cublas_set_matrix( &
            ndof_per_element, &
            nels_pp, &
            sizeof(0.0d0), &
            pmul_pp, &
            ndof_per_element, &
            device_lhs_vectors, &
            ndof_per_element)
       if (status .ne. 0) then
          print *, "Failed to copy data to gpu!"
          status = cublas_shutdown
          stop
       end if

    !RZ Copying p_pp and u_pp data over to GPU initially 
    !   
    !   status = cublas_set_vector( &
    !            neq_pp, &
    !            sizeof(0.d0), &
    !            p_pp, &
    !            1, &
    !            device_p_pp, &
    !            1)
    !
    !   status = cublas_set_vector( &
    !            neq_pp, &
    !            sizeof(0.d0), &
    !            u_pp, &
    !            1, &
    !            device_u_pp, &
    !            1)

       t_start2 = elap_time() 

       ! Call CUBLAS 
       alpha = 1.d0
       beta = 0.d0
       op = "N"
       call cublas_dgemm( &
            trim(op)//c_null_char, &
            trim(op)//c_null_char, &
            ndof_per_element, &
            nels_pp, &
            ndof_per_element, &
            alpha, &
            device_matrix, &
            ndof_per_element, &
            device_lhs_vectors, &
            ndof_per_element, &
            beta, &
            device_rhs_vectors, &
            ndof_per_element)

       ! Set to true to measure raw kernel execution time
       if (.true.) then

          status = sync_gpu()
          if (status .ne. 0) then
             print *, "Failed to sync gpu!"
             status = cublas_shutdown
             stop
          end if
       end if

       ! Accumulate time for computation
       t_rawcomp = t_rawcomp + (elap_time() - t_start2)

       ! Copy result vectors back from gpu
       status = cublas_get_matrix( &
            ndof_per_element, & 
            nels_pp, &
            sizeof(0.0d0), &
            device_rhs_vectors, &
            ndof_per_element, & 
            utemp_pp, &
            ndof_per_element)
       if (status .ne. 0) then
          print *, "Failed to copy data from gpu!"
          status = cublas_shutdown
          stop
       end if

    end if

    ! Accumulate total time for matrix multiply
    t_comp = t_comp + (elap_time() - t_start1)

    CALL scatter(u_pp,utemp_pp)

    IF(fixed_freedoms_pp > 0) THEN
      DO i = 1, fixed_freedoms_pp
        j       = no_pp(i) - ieq_start + 1
        u_pp(j) = p_pp(j) * store_pp(i)
      END DO
    END IF
                                       
    !RZ Copying p_pp and u_pp over for use in iteration                                                                                                 

    !   status = cublas_set_vector( &
    !            neq_pp, &
    !            sizeof(0.d0), &
    !            p_pp, &
    !            1, &
    !            device_p_pp, &
    !            1)
    !
    !   status = cublas_set_vector( &
    !            neq_pp, &
    !            sizeof(0.d0), &
    !            u_pp, &
    !            1, &
    !            device_u_pp, &
    !            1)                                  

!   up      = DOT_PRODUCT_P(r_pp,d_pp)
          
    !RZ Alpha and dot product on CPU and output
    if (.not. use_gpu) then

      timest(19) = elap_time()

      alpha   = up0/DOT_PRODUCT_P(p_pp,u_pp)
    
      timest(20) = timest(20) + (elap_time()-timest(19))

    else
    !RZ Alpha and dot product on GPU and output
    !RZ Copying p_pp and u_pp over for use in iteration                                                                                                 

       timest(17) = elap_time()

       status = cublas_set_vector( &
                neq_pp, &
                sizeof(0.d0), &
                p_pp, &
                1, &
                device_p_pp, &
                1)

       status = cublas_set_vector( &
                neq_pp, &
                sizeof(0.d0), &
                u_pp, &
                1, &
                device_u_pp, &
                1)                                 
 
      timest(18) = timest(18) + (elap_time()-timest(17))

      timest(19) = elap_time()

      alpha    = up0/cublas_Ddot(neq_pp,device_p_pp,1,device_u_pp,1)

      timest(20) = timest(20) + (elap_time()-timest(19))
    
   ! WRITE (*,*) numpe,alpha    
    endif

    xnew_pp = x_pp + p_pp*alpha
    r_pp    = r_pp - u_pp*alpha
    d_pp    = diag_precon_pp*r_pp
    up1     = DOT_PRODUCT_P(r_pp,d_pp)
!   beta    = DOT_PRODUCT_P(r_pp,d_pp)/up0
    beta    = up1/up0
    up0     = up1
    p_pp    = d_pp + p_pp*beta  

    CALL checon_par(xnew_pp,tol,converged,x_pp)    
    IF(converged.OR.iters==limit)EXIT

  END DO iterations

  
  write (*,*) "Total time in matrix-vector multiply:" , t_comp
  if (use_gpu) then
     write(*,*) &
     "Time in matrix-vector multiply excluding vector transfer", &
     t_rawcomp
  end if

  DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,km,pmul_pp) 

  ! Tidy up GPU
  status = cublas_free(device_matrix)
  status = cublas_free(device_lhs_vectors)
  status = cublas_free(device_rhs_vectors)
  status = cublas_free(device_matrix)
  status = cublas_shutdown

  timest(13) = elap_time()

!------------------------------------------------------------------------------
! 15. Print out displacements
!------------------------------------------------------------------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  
  IF(numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ")-1)//".bin.ensi.DISPL-000001"
    OPEN(24, file=fname, status='replace', action='write',                    & 
         form='unformatted',access='stream')
    cbuffer = "Alya Ensight Gold --- Vector per-node variable file"
    WRITE(24) cbuffer
    cbuffer = "part"        ;  WRITE(24) cbuffer
    WRITE(24) int(1,kind=c_int)
    cbuffer = "coordinates" ; WRITE(24) cbuffer
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

  ALLOCATE(tempres(nodes_pp))
  tempres = zero

  DO i=1,ndim
    tempres = zero
    DO j=1,nodes_pp
      k=i+(ndim*(j-1))
      tempres(j)=disp_pp(k)
    END DO
    CALL dismsh_ensi_pb2(24,1,nodes_pp,npes,numpe,1,tempres)
  END DO

! CALL write_nodal_variable(label,24,1,nodes_pp,npes,numpe,ndim,disp_pp)

  DEALLOCATE(disp_pp)
  DEALLOCATE(tempres)

  IF(numpe==1) CLOSE(24)

!------------------------------------------------------------------------------
! 16b. Print out stress, principal stress and reactions
!------------------------------------------------------------------------------

  IF(output_stress) THEN

  IF(numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ")-1) // ".str"
    OPEN(25, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".pri"
    OPEN(26, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".vms"
    OPEN(27, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".rea"
    OPEN(28, file=fname, status='replace', action='write')
  END IF

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

  DO iel = 1,nels_pp
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
  CALL write_nodal_variable(label,25,1,nodes_pp,npes,numpe,nst,               &
                            stressnodes_pp)
                            
  DEALLOCATE(stress_integral_pp,stressnodes_pp)

  IF(numpe==1) CLOSE(25)

!------------------------------------------------------------------------------
! 16d. Principal stress
!------------------------------------------------------------------------------
  
  label = "*PRINCIPAL STRESS"
  
  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,          &
                        node_start,node_end,shape_integral_pp,                &
                        principal_integral_pp,princinodes_pp)
  CALL write_nodal_variable(label,26,1,nodes_pp,npes,numpe,nodof,             &
                            princinodes_pp)
                            
  DEALLOCATE(principal_integral_pp)

  IF(numpe==1) CLOSE(26)

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

  CALL write_nodal_variable(label,27,1,nodes_pp,npes,numpe,nodof,             &
                            princinodes_pp)
                            
  DEALLOCATE(princinodes_pp)

  IF(numpe==1) CLOSE(27)

!------------------------------------------------------------------------------
! 16f. Reactions
!------------------------------------------------------------------------------

  label = "*NODAL REACTIONS"
  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,             &
                     node_start,node_end,utemp_pp,reacnodes_pp,0)
  CALL write_nodal_variable(label,28,1,nodes_pp,npes,numpe,nodof,             &
                            reacnodes_pp)
  DEALLOCATE(reacnodes_pp,shape_integral_pp)

  IF(numpe==1) CLOSE(28)

!------------------------------------------------------------------------------

  END IF

!------------------------------------------------------------------------------
! 16g. End timer
!------------------------------------------------------------------------------

  timest(14) = elap_time()
  timest(15) = t_rawcomp
  timest(16) = t_comp
   
!------------------------------------------------------------------------------
! 17. Output performance data
!------------------------------------------------------------------------------

  CALL WRITE_XX9(fixed_freedoms,iters,job_name,loaded_nodes,neq,nn,npes,nr,  &
                  numpe,timest,tload,use_gpu)
 
  CALL shutdown() 
 
END PROGRAM XX9
