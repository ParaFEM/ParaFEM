PROGRAM xx3        
!------------------------------------------------------------------------------ 
!      Program XX.3 Three dimensional analysis of an elastic solid using load
!                   control or displacement control. GPU support using CUDA.
!                   See program P121.
!------------------------------------------------------------------------------ 
                                 
  USE precision     ; USE global_variables ; USE mp_interface
  USE input         ; USE output           ; USE loading
  USE timing        ; USE maths            ; USE gather_scatter
  USE partition     ; USE elements         ; USE steering        ; USE pcg
  USE new_library

  USE iso_c_binding
  
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
  REAL(iwp)             :: e,v,det,tol,up,alpha,beta,tload
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp
  REAL(iwp),PARAMETER   :: penalty = 1.0e20_iwp
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: program_name='xx3'
  CHARACTER(LEN=50)     :: fname,job_name,label
  LOGICAL               :: converged = .false.
  INTEGER               :: dsize = sizeof(0.0d0)   ! Size of C double

  ! GPU related variables etc
  ! -------------------------
  ! Pointers to device memory
  type (c_ptr) :: device_matrix
  type (c_ptr) :: device_lhs_vector
  type (c_ptr) :: device_rhs_vector

  logical :: use_gpu = .true.
  integer :: status
  integer :: ndof_per_element
  
  real(iwp) :: t_rawcomp
  real(iwp) :: t_comp
  real(iwp) :: t_start1
  real(iwp) :: t_start2

  interface
     integer(c_int) function set_gpu(device_id) bind(C)
       
       use iso_c_binding
       
       integer(c_int) :: device_id
     end function set_gpu

     integer(c_int) function allocate_memory_on_gpu( &
          n_elements, &
          element_size, &
          device_pointer) bind(C)
       
       use iso_c_binding
       
       integer(c_int) :: n_elements
       integer(c_int) :: element_size
       type (c_ptr) :: device_pointer
     end function allocate_memory_on_gpu

     
     integer(c_int) function free_memory_on_gpu( &
          device_pointer) bind(C)

       use iso_c_binding

       type (c_ptr) :: device_pointer
     end function free_memory_on_gpu

     
     integer(c_int) function copy_3d_data_to_gpu( & 
          n_elements, &
          element_size, &
          host_data, &
          device_pointer) bind(C, name='copy_data_to_gpu')
       
       use iso_c_binding
       
       integer(c_int) :: n_elements
       integer(c_int) :: element_size
       real(c_double) :: host_data(*)
       type (c_ptr) :: device_pointer
     end function copy_3d_data_to_gpu
     
     
     integer(c_int) function copy_2d_data_to_gpu( & 
          n_elements, &
          element_size, &
          host_data, &
          device_pointer) bind(C, name='copy_data_to_gpu')
       
       use iso_c_binding
       
       integer(c_int) :: n_elements
       integer(c_int) :: element_size
       real(c_double) :: host_data(*)
       type (c_ptr) :: device_pointer
     end function copy_2d_data_to_gpu
     

     integer(c_int) function copy_data_from_gpu( &
          n_elements, &
          element_size, &
          host_data, &
          device_pointer) bind(C)

       use iso_c_binding

       integer(c_int) :: n_elements
       integer(c_int) :: element_size
       real(c_double) :: host_data(*)
       type (c_ptr) :: device_pointer
     end function copy_data_from_gpu

     
     integer(c_int) function matrix_vector_multiplies( &
          n_mat, &
          n_row, &
          n_col, &
          d_lhs_vector, &
          d_matrix, &
          d_rhs_vector) bind(C)
       
       use iso_c_binding

       integer(c_int) :: n_mat
       integer(c_int) :: n_row
       integer(c_int) :: n_col
       type (c_ptr) :: d_lhs_vector
       type (c_ptr) :: d_matrix
       type (c_ptr) :: d_rhs_vector
     end function matrix_vector_multiplies
  end interface

  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:),val(:,:),disp_pp(:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:)
  REAL(iwp),ALLOCATABLE :: bee(:,:),storkm_pp(:,:,:),eld(:),eps(:),sigma(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:),xnew_pp(:)
  REAL(iwp),ALLOCATABLE :: u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),timest(:)
  REAL(iwp),ALLOCATABLE :: diag_precon_tmp(:,:),eld_pp(:,:),tensor_pp(:,:,:)
  REAL(iwp),ALLOCATABLE :: valf(:),store_pp(:)
  REAL(iwp),ALLOCATABLE :: fun(:),shape_integral_pp(:,:)
  REAL(iwp),ALLOCATABLE :: stress_integral_pp(:,:),stressnodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal_integral_pp(:,:),princinodes_pp(:)
  REAL(iwp),ALLOCATABLE :: principal(:),reacnodes_pp(:)  
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

! CALL read_p121(job_name,numpe,e,element,fixed_freedoms,limit,loaded_nodes, &
!                meshgen,nels,nip,nn,nod,nr,partitioner,tol,v)
 
  fixed_freedoms = 0 
  CALL read_p121(job_name,numpe,e,element,limit,loaded_nodes,                &
                 meshgen,nels,nip,nn,nod,nr,partitioner,tol,v)

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
  
  CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)
  timest(3) = elap_time()

! CALL read_g_num_pp(job_name,iel_start,nels,nn,numpe,g_num_pp)

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
           storkm_pp(ntot,ntot,nels_pp),eld(ntot),eps(nst),sigma(nst),        &
           pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),                      &
           weights(nip),g_g_pp(ntot,nels_pp),fun(nod))

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------

  CALL rearrange(rest)
  
  g_g_pp = 0

  elements_1: DO iel = 1, nels_pp
!   CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
    CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
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
!------------------------------------------------------------------------------

  dee = zero
  CALL deemat(dee,e,v)
  CALL sample(element,points,weights)
 
  storkm_pp       = zero
 
  elements_3: DO iel=1,nels_pp
    gauss_pts_1: DO i=1,nip
      CALL shape_der(der,points,i)
      jac   = MATMUL(der,g_coord_pp(:,:,iel))
      det   = determinant(jac)
      CALL invert(jac)
      deriv = MATMUL(jac,der)
      CALL beemat(deriv,bee)
      storkm_pp(:,:,iel)   = storkm_pp(:,:,iel) +                             &
                             MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) *         &
                             det*weights(i)   
    END DO gauss_pts_1
  END DO elements_3
  
  timest(10) = elap_time()

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

!------------------------------------------------------------------------------
! 10. Read in fixed nodal displacements and assign to equations
!------------------------------------------------------------------------------

  IF(fixed_freedoms > 0) THEN

    ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),valf(fixed_freedoms),    &
             no_pp_temp(fixed_freedoms),sense(fixed_freedoms))

    node = 0 ; no = 0 ; no_pp_temp = 0 ; sense = 0 ; valf = zero

    CALL read_fixed(job_name,numpe,node,sense,valf)
    CALL find_no(node,rest,sense,no)
!   CALL reindex_fixed_nodes(ieq_start,no,no_pp_temp,fixed_freedoms_pp,       &
!                            fixed_freedoms_start,neq_pp)
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
     
     ndof_per_element = 3 * nod

     ! Allocate memory on the gpu for matrices
     status = allocate_memory_on_gpu( &
          nels_pp*ndof_per_element**2, &
          dsize, &
          device_matrix)
     if (status > 0) then
        print *, "gpu memory failed to allocate!"
        stop
     end if

     ! Allocate memory for lhs vectors
     status = allocate_memory_on_gpu( &
          nels_pp*ndof_per_element, &
          dsize, &
          device_lhs_vector)
     if (status > 0) then
        print *, "gpu memory failed to allocate!"
        stop
     end if

     ! Allocate memory for rhs vectors
     status =  allocate_memory_on_gpu( &
          nels_pp*ndof_per_element, &
          dsize, &
          device_rhs_vector)
     if (status > 0) then
        print *, "gpu memory failed to allocate!"
        stop
     end if
     
     ! Copy matrix data to the gpu
     status = copy_3d_data_to_gpu( &
          nels_pp*ndof_per_element**2, &
          dsize, &
          storkm_pp(:,:,:), &
          device_matrix)
     if (status > 0) then
        print *, "Failed to copy data to gpu!"
        stop
     end if

   end if
  
!------------------------------------------------------------------------------
! 14. Preconditioned conjugate gradient iterations
!------------------------------------------------------------------------------

  iters = 0

  t_comp = 0
  t_rawcomp = 0

  iterations: DO
   
    iters    = iters + 1
    u_pp     = zero
    pmul_pp  = zero
    utemp_pp = zero
 
    CALL gather(p_pp,pmul_pp)

    t_start1 = elap_time()
    
    ! matmul version
    if (.not. use_gpu) then

       elements_5: DO iel=1,nels_pp
          utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
       END DO elements_5

       ! gpu version    
    else

       ! Copy lhs vectors to gpu
       status = copy_2d_data_to_gpu( &
            nels_pp*ndof_per_element, &
            dsize, &
            pmul_pp(:,1:nels_pp), &
            device_lhs_vector)
       if (status > 0) then
          print *, "Failed to copy data to gpu!"
          stop
       end if

       t_start2 = elap_time()
             
       ! Call matrix-vector multiply kernel
       status = matrix_vector_multiplies( &
            nels_pp, &
            ndof_per_element, &
            ndof_per_element, &
            device_lhs_vector, &
            device_matrix, &
            device_rhs_vector)
       if (status > 0) then
          print *, "gpu failed to multiply!"
          stop
       end if

       ! Accumulate time for computation
       t_rawcomp = t_rawcomp + (elap_time() - t_start2)

       ! Copy result vector back from gpu
       status = copy_data_from_gpu( &
            nels_pp*ndof_per_element, & 
            dsize, &
            utemp_pp(:,1:nels_pp), &
            device_rhs_vector)
       if (status > 0) then
          print *, "Failed to copy data from gpu!"
          stop
       end if

    end if

    ! Accumulate time for matrix multiply
    t_comp = t_comp + (elap_time() - t_start1)

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

    CALL checon_par(xnew_pp,tol,converged,x_pp)    
    IF(converged.OR.iters==limit)EXIT

  END DO iterations

  write(*,*) "Total time in matrix-vector multiply:", t_comp
  if (use_gpu) then
     write(*,*) &
     "Time in matrix-vector multiply excluding vector transfer", &
     t_rawcomp
  end if

  DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,storkm_pp,pmul_pp) 
 
  ! Deallocate memory on GPU
  if (use_gpu) then
     status = free_memory_on_gpu(device_matrix)
     if (status > 0) then
        print *, "gpu memory failed to deallocate device_matrix!"
        stop
     end if
     
     status = free_memory_on_gpu(device_lhs_vector)
     if (status > 0) then
        print *, "gpu memory failed to deallocate device_lhs_vector!"
        stop
     end if

     status =  free_memory_on_gpu(device_rhs_vector)
     if (status > 0) then
        print *, "gpu memory failed to deallocate device_rhs_vector!"
        stop
     end if
  end if

  timest(13) = elap_time()

!------------------------------------------------------------------------------
! 15. Print out displacements, stress, principal stress and reactions
!------------------------------------------------------------------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  
  IF(numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ")-1)//".dis"
    OPEN(24, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".str"
    OPEN(25, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".pri"
    OPEN(26, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".vms"
    OPEN(27, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".rea"
    OPEN(28, file=fname, status='replace', action='write')
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
  CALL write_nodal_variable(label,24,1,nodes_pp,npes,numpe,ndim,disp_pp)

  DEALLOCATE(disp_pp)

  IF(numpe==1) CLOSE(24)

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

  DO iel = 1,nels_pp
    DO i = 1,nip
      CALL shape_der(der,points,i)
      jac   = MATMUL(der,g_coord_pp(:,:,iel))
      det   = DETERMINANT(jac) 
      CALL invert(jac)
      deriv = MATMUL(jac,der)
      CALL beemat(deriv,bee)
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

  timest(14) = elap_time()
   
!------------------------------------------------------------------------------
! 17. Output performance data
!------------------------------------------------------------------------------

  CALL WRITE_P121(fixed_freedoms,iters,job_name,loaded_nodes,neq,nn,npes,nr,  &
                  numpe,timest,tload)

  if (use_gpu) then
     print*, "gpu acceleration used in solver"
  else
     print*, "gpu acceleration not used in solver"
  end if
 
  CALL shutdown() 
 
END PROGRAM XX3
