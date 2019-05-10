PROGRAM XX10       
!------------------------------------------------------------------------------ 
!      Program XX.10 Three dimensional analysis of an elastic solid using load
!                    control or displacement control. GPU support using OpenCL.
!                    Special case with single element stiffness matrix and 
!                    matrix-matrix kernal in PCG.
!                    See program P121.
!------------------------------------------------------------------------------ 
                                 
  USE precision     ; USE global_variables ; USE mp_interface
  USE input         ; USE output           ; USE loading
  USE timing        ; USE maths            ; USE gather_scatter
  USE steering      ; USE new_library
  
  ! C types for GPU code
  ! --------------------
  USE, INTRINSIC :: ISO_C_BINDING
  
  IMPLICIT NONE

  ! C functions for GPU code
  ! ------------------------
  INTEGER, EXTERNAL :: init_opencl
  INTEGER, EXTERNAL :: stop_opencl
  INTEGER, EXTERNAL :: allocate_memory_on_gpu
  INTEGER, EXTERNAL :: free_memory_on_gpu
  INTEGER, EXTERNAL :: copy_data_to_gpu
  INTEGER, EXTERNAL :: copy_data_from_gpu
  INTEGER, EXTERNAL :: compile_kernel_from_file
  INTEGER, EXTERNAL :: matrix_matrix_multiplies_1d
  INTEGER, EXTERNAL :: matrix_matrix_multiplies_2d
  INTEGER, EXTERNAL :: blas_dgemm
  INTEGER, EXTERNAL :: c_compare_vecs

  ! GPU Code config
  ! ---------------
  INTEGER       :: mult_method           ! Read in at the command line
                                         ! Method to use to do the matrix-matrix mult
                                         ! 0 - CPU: original code
                                         ! 1 - GPU: our own matmul 1D kernel (naive)
                                         ! 2 - GPU: our own matmul 2D kernel (naive)
                                         ! 3 - GPU: AMD clBlas

  LOGICAL       :: check_gpu = .false.   ! Compare GPU results to intrinsic matmul results (slow)
  INTEGER       :: gpu_device = 1        ! Set to 0 to use OpenCL on CPU (AMD only)

  ! Our own kernel source file (must be in current directory when executing xx10)
  CHARACTER(LEN=32) :: srcfilename = "xx10.cl"//CHAR(0)
  CHARACTER(LEN=32), dimension(2) :: kernelnames = (/ &
       "MatMatMultiply_col_1d"//CHAR(0),      &
       "MatMatMultiply_col_2d"//CHAR(0)       &
  /)

  ! End of GPU Code config
  ! ----------------------

  ! GPU related variables
  ! ---------------------
  ! Matrices are square so could use one var for x and y vectors but code here is general.
  ! X,Y refers to the matrix, vector multiplication of the form Y = M.X
  ! All sizes except dsize are in terms of number of elements, not bytes.
  LOGICAL       :: use_amdblas             ! Will be set later. Don't set here.
  INTEGER       :: dsize = sizeof(0.0d0)   ! Size of matrix, vector elements (i.e, C double)
  INTEGER       :: matsize                 ! Size of the single matrix
  INTEGER       :: vecsize_lhs             ! Size of the entire Y vector array
  INTEGER       :: vecsize_rhs             ! Size of the entire X vector array
  INTEGER       :: memflag                 ! GPU mem type: 0 Read only, 1 write only, 2 read and write
  INTEGER       :: status                  ! OpenCL return flag
  REAL(iwp)     :: misc_timers(4)          ! Time code around GPU work and matmuls

  ! Pointers to device memory (store a cl_mem as a void*)
  type (c_ptr)  :: device_matrix      = C_NULL_PTR
  type (c_ptr)  :: device_lhs_vectors = C_NULL_PTR
  type (c_ptr)  :: device_rhs_vectors = C_NULL_PTR

!------------------------------------------------------------------------------ 
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

 ! neq,ntot are now global variables - not declared

  INTEGER,PARAMETER     :: nodof=3,ndim=3,nst=6
  INTEGER               :: loaded_nodes,fixed_freedoms,iel,i,j,k,l,idx1,idx2
  INTEGER               :: iters,limit,nn,nr,nip,nod,nels,ndof,npes_pp
  INTEGER               :: node_end,node_start,nodes_pp
  INTEGER               :: argc,iargc,nlen,meshgen,partitioner
  INTEGER               :: fixed_freedoms_pp,fixed_freedoms_start
  REAL(iwp)             :: e,v,det,tol,up,alpha,beta,tload
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp
  REAL(iwp),PARAMETER   :: penalty = 1.0e20_iwp
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: program_name='xx10'
  CHARACTER(LEN=50)     :: fname,job_name,label
  LOGICAL               :: converged = .false.

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
  REAL(iwp),ALLOCATABLE :: principal(:),reacnodes_pp(:)  
  INTEGER,  ALLOCATABLE :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)
  INTEGER,  ALLOCATABLE :: no(:),no_pp(:),no_pp_temp(:),sense(:)

  ! For OpenCL-testing
  REAL(iwp),ALLOCATABLE :: check_utemp_pp(:,:)

!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions. 
!------------------------------------------------------------------------------
 
  ALLOCATE(timest(20))
  timest    = zero 
  timest(1) = elap_time()
  
  CALL find_pe_procs(numpe,npes)

  CALL getname2(job_name,nlen,mult_method)

  CALL read_p121(job_name,numpe,e,element,limit,loaded_nodes, &
                 meshgen,nels,nip,nn,nod,nr,partitioner,tol,v)

  fixed_freedoms = 0 ! fixed_freedoms not read in via READ_P121

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

! For C-testing
  IF ( mult_method > 0 .AND. check_gpu ) THEN
     ALLOCATE( check_utemp_pp(ntot,nels_pp) )
     print *, "Will check GPU results against CPU results. This is SLOW!"
  END IF

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
      diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + km(i,i)
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
  misc_timers = zero

  ! Set up GPU
  IF ( mult_method > 0 ) THEN

     ! Flags for type of matmul based on method
     use_amdblas = (mult_method==3)

     print *, "Initializing GPU in process ", numpe, " Using Blas: ", use_amdblas
     misc_timers(1) = elap_time();
     
     ! Matrix and arrays-of-vectors sizes
     matsize     = ntot**2
     vecsize_lhs = ntot*nels_pp
     vecsize_rhs = ntot*nels_pp

     ! Create an OpenCL context on the device
     status = init_opencl( gpu_device, numpe-1, use_amdblas )
     !status = init_opencl( gpu_device, -1, use_amdblas )
     if ( status /= 1 ) then
        print *, "Failed to init OpenCL"
        stop
     end if

     ! Allocate read-only device memory (kernel won't write) for matrix data
     memflag = 0
     status = allocate_memory_on_gpu( matsize, dsize, memflag, device_matrix );
     if ( status /= 1 ) then
        print *, "Failed to allocate device memory (matrix)"
        stop
     end if

     ! Allocate read-only device memory (kernel won't write) for rhs vector 
     memflag = 0
     status = allocate_memory_on_gpu( vecsize_rhs, dsize, memflag, device_rhs_vectors )
     if ( status /= 1 ) then
        print *, "Failed to allocate device memory (rhs vector)"
        stop
     end if

     ! Allocate at least write (and possibly read) memory for lhs vector (the result vector)
     IF ( use_amdblas ) THEN
        memflag = 2     ! AMD clBlas can read and write the the result vector
     ELSE
        memflag = 1     ! Our own kernel only writes to the result vector
     END IF
     status = allocate_memory_on_gpu( vecsize_lhs, dsize, memflag, device_lhs_vectors )
     IF ( status /= 1 ) THEN
        print *, "Failed to allocate device memory (lhs vector)"
        stop
     END IF

     ! Copy coefficient matrix km to the GPU
     status = copy_data_to_gpu( matsize, dsize, km, device_matrix )
     IF ( status /= 1 ) THEN
        print *, "Error copying entire storkm_pp array of matrices to the GPU"
        stop
     END IF

     ! Compile our own kernels
     IF ( use_amdblas == 0 ) THEN
        status = compile_kernel_from_file( srcfilename, kernelnames(mult_method), C_NULL_PTR );
        IF ( status /= 1 ) THEN
           print *, "Error compiling kernel ", kernelnames(mult_method)
           stop
        END IF
     END IF

     misc_timers(2) = elap_time() - misc_timers(1);
     write(*,*) 'Time for GPU init + host-to-device copy: ', misc_timers(2)

  END IF
!------------------------------------------------------------------------------
! 14. Preconditioned conjugate gradient iterations
!------------------------------------------------------------------------------

  iters = 0
  
  iterations: DO 
    iters    = iters + 1
    u_pp     = zero
    pmul_pp  = zero
    utemp_pp = zero
    u_pp     = zero

    CALL gather(p_pp,pmul_pp)

    IF ( mult_method == 0 ) THEN
       misc_timers(3) = elap_time();

       ! Original cpu version
       utemp_pp = MATMUL(km,pmul_pp)

       ! Accumulate time for only the matmul code
       misc_timers(4) = misc_timers(4) + (elap_time() - misc_timers(3));

    ELSE IF ( mult_method >= 1 .OR. mult_method <= 3 ) THEN

       ! gpu versions
       ! ------------

       ! Time GPU kernel and data upload/download
       misc_timers(1) = elap_time();

       ! Transfer the current RHS vectors to device, do matmul, transfer result vectors back
       status = copy_data_to_gpu(vecsize_rhs, dsize, pmul_pp, device_rhs_vectors )

       ! Time only the GPU kernel
       misc_timers(3) = elap_time()
       IF ( mult_method == 1 ) THEN
          ! Multiply using our own 1D kernel
          status = status + matrix_matrix_multiplies_1d( ntot, nels_pp, ntot, &
                                                         device_matrix, device_rhs_vectors, device_lhs_vectors )
       ELSE IF ( mult_method == 2 ) THEN
          ! Multiply using our own 2D kernel
          status = status + matrix_matrix_multiplies_2d( ntot, nels_pp, ntot, &
                                                         device_matrix, device_rhs_vectors, device_lhs_vectors )
       ELSE
          ! Multiply using AMD Blas
          status = status + blas_dgemm( ntot, nels_pp, ntot, device_matrix, device_rhs_vectors, device_lhs_vectors )
       END IF
       misc_timers(4) = misc_timers(4) + (elap_time()-misc_timers(3));

       status = status + copy_data_from_gpu(vecsize_lhs, dsize, device_lhs_vectors, utemp_pp )
       
       misc_timers(2) = misc_timers(2) + (elap_time()-misc_timers(1));

       IF ( status /= 3 ) THEN
          write(*,*) 'Error in OpenCL calls: ', status, ' of 3 succeeded'
          stop
       END IF
                 
       ! Compare GPU matmul results to CPU results
       IF( check_gpu ) THEN             
          check_utemp_pp = MATMUL(km,pmul_pp)
          IF( c_compare_vecs(vecsize_lhs, check_utemp_pp, utemp_pp) == 0 ) THEN
             write(*,*) 'Iteration ', iters, ' GPU/CPU result vectors differ'
          END IF
       END IF

    ELSE
       write(*,*) 'Must set mult_method to 0,1,2,3'
       stop                
    END IF

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

  write(*,*) 'Total iters time CPU matmul() or GPU kernel: ', (misc_timers(4))
  write(*,*) 'Total iters time GPU kernel+transfer:        ', (misc_timers(2))

  DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,km,pmul_pp) 

  ! Cleanup GPU
  IF ( mult_method > 0 ) THEN
     IF( check_gpu ) THEN
        DEALLOCATE( check_utemp_pp )
     END IF
     status = free_memory_on_gpu( device_matrix )
     status = status + free_memory_on_gpu( device_lhs_vectors )
     status = status + free_memory_on_gpu( device_rhs_vectors )
     status = status + stop_opencl()
     IF( status /= 4 ) THEN
        write(*,*) 'Error in OpenCL cleanup calls: ', status, ' of 4 succeeded'
        stop
     END IF
  END IF
  
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

  CALL WRITE_XX10(fixed_freedoms,iters,job_name,loaded_nodes,mult_method,     &
                  neq,nn,npes,nr,numpe,timest,tload)
 
  CALL shutdown() 
 
END PROGRAM XX10
