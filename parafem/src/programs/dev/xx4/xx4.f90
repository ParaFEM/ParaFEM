PROGRAM xx4        
!------------------------------------------------------------------------------ 
!      Program XX.4 Three dimensional analysis of an elastic solid using load
!                   control or displacement control. GPU support using OpenCL.
!                   See program P121.
!------------------------------------------------------------------------------ 
                                 
  USE precision     ; USE global_variables ; USE mp_interface
  USE input         ; USE output           ; USE loading
  USE timing        ; USE maths            ; USE gather_scatter
  USE steering      ; USE new_library


  ! From xx16 on 27march
  USE omp_lib
  
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
  INTEGER, EXTERNAL :: matrix_vector_multiplies
  INTEGER, EXTERNAL :: blas_dgemv
  INTEGER, EXTERNAL :: mem_copy_blas_dgemv_async
  INTEGER, EXTERNAL :: c_compare_vecs

  ! GPU Code config
  ! ---------------
  INTEGER       :: mult_method           ! Read in at command line
                                         ! Method to use to do the matrix-vector mult
                                         ! 0 - CPU: original code 
                                         ! 1 - CPU: voxels as 8 node bricks
                                         ! 2 - GPU: our own matmul kernel 1 (naive)
                                         ! 3 - GPU: our own matmul kernel 2 (local mem)
                                         ! 4 - GPU: AMD clBlas (synchronously)
                                         ! 5 - GPU: AMD clBlas (asynchronously, experimental)
                                         

  LOGICAL       :: check_gpu = .false.   ! Compare GPU results to intrinsic matmul results (slow)
  INTEGER       :: gpu_device = 1        ! Set to 0 to use OpenCL on CPU (AMD only)

  ! Our own kernel source file (must be in current directory when executing xx4)
  CHARACTER(LEN=32) :: srcfilename = "xx4.cl"//CHAR(0)

  ! The mult_method indexes in to this array of kernel names (currently for methods 1 and 2)
  CHARACTER(LEN=32), dimension(2) :: kernelnames = (/ &
       "MultiMatVecMultiply_col_naive"//CHAR(0),      &
       "MultiMatVecMultiply_col_local"//CHAR(0)       &
  /)

  ! End of GPU Code config
  ! ----------------------

  ! GPU related variables
  ! ---------------------
  ! Matrices are square so could use one var for x and y vectors but code here is general.
  ! X,Y refers to the matrix, vector multiplication of the form Y = M.X
  ! All sizes except dsize are in terms of number of elements, not bytes.
  LOGICAL       :: use_kernels             ! Types of mat mul according to method above
  LOGICAL       :: use_amdblas             ! Will be set later. Don't set here.
  INTEGER       :: dsize = sizeof(0.0d0)   ! Size of matrix, vector elements (i.e, C double)
  INTEGER       :: matsize                 ! Size of one matrix
  INTEGER       :: vecsize_lhs             ! Size of one Y vector
  INTEGER       :: vecsize_rhs             ! Size of one X vector
  INTEGER       :: msize                   ! Temp working size of matrix
  INTEGER       :: vsize_lhs               ! Temp working size of Y vector
  INTEGER       :: vsize_rhs               ! Temp working size of X vector
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
  INTEGER               :: argc,iargc,meshgen,partitioner
  INTEGER               :: fixed_freedoms_pp,fixed_freedoms_start
  REAL(iwp)             :: e,v,det,tol,up,alpha,beta,tload
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp
  REAL(iwp),PARAMETER   :: penalty = 1.0e20_iwp
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: program_name='xx4'
  CHARACTER(LEN=50)     :: fname,job_name,label,argv
  LOGICAL               :: converged = .false.

 ! VARIABLES FROM XX16 on 27march not sure about argv
  INTEGER               :: testcase
  INTEGER               :: nlen,section,fnum
  REAL(iwp),PARAMETER   :: one = 1.0_iwp
  REAL(iwp)             :: gigaflops,gigaflops_1,gigaflops_2
  REAL(iwp)             :: flop,flop_1,flop_2
  REAL(iwp)             :: DDOT ! BLAS MKL
  CHARACTER(LEN=6)      :: ch
  CHARACTER(LEN=80)     :: cbuffer
  CHARACTER(LEN=80)     :: message

 ! Switches for various PCG options XX16 27march
  LOGICAL               :: io_binary    = .false.
  LOGICAL               :: sys_storkm   = .false.
  LOGICAL               :: vox_storkm   = .false.
  LOGICAL               :: write_stress = .false.
 
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
  REAL(iwp),ALLOCATABLE :: timest_min(:),timest_max(:),temp(:),pad1(:),pad2(:),km(:,:)
  REAL(iwp),ALLOCATABLE :: vstorkm_pp(:,:)    ! store km as vector
  REAL(iwp),ALLOCATABLE :: vox_storkm_pp(:,:) ! minimal storage of km
  REAL(iwp),ALLOCATABLE :: vtemp(:)           ! temporary array
  REAL(iwp),ALLOCATABLE :: ptemp(:)           ! temporary array
  INTEGER,ALLOCATABLE   :: iv(:,:)            ! index for vox value

  ! For OpenCL-testing
  REAL(iwp),ALLOCATABLE :: check_utemp_pp(:,:)
  
!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions. 
!------------------------------------------------------------------------------
  
  ALLOCATE(timest(24))
  timest    = zero 
  timest(1) = elap_time()
  
  CALL find_pe_procs(numpe,npes)

! argc = iargc()
! IF (argc /= 1) CALL job_name_error(numpe,program_name)
! CALL GETARG(1, job_name) 

  CALL getname2(argv,nlen,mult_method)

! CALL read_p121(job_name,numpe,e,element,fixed_freedoms,limit,loaded_nodes, &
!                meshgen,nels,nip,nn,nod,nr,partitioner,tol,v)

  CALL read_p121(argv,numpe,e,element,limit,loaded_nodes, &
                 meshgen,nels,nip,nn,nod,nr,partitioner,tol,v)

  fixed_freedoms = 0 ! fixed_freedoms not read in via READ_P121

  CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)

  ndof = nod*nodof
  ntot = ndof

  ALLOCATE(g_num_pp(nod, nels_pp)) 
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp)) 
  ALLOCATE(rest(nr,nodof+1)) 
 
  g_num_pp  = 0
  g_coord_pp= zero
  rest      = 0

  timest(2) = elap_time()
  
  CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
  timest(3) = elap_time()

  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  timest(4) = elap_time()

  CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
  timest(5) = elap_time()

  CALL read_rest(argv,numpe,rest)
  timest(6) = elap_time()

  IF(numpe==1) PRINT *, "Read input data"
    
!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),principal(ndim),      &
           der(ndim,nod),deriv(ndim,nod),bee(nst,ntot),                       &
           eld(ntot),eps(nst),sigma(nst),                                     &
           weights(nip),g_g_pp(ntot,nels_pp),fun(nod))

! For C-testing
  IF ( check_gpu ) THEN
     ALLOCATE( check_utemp_pp(ntot,nels_pp) )
     print *, "Will check GPU results against CPU results. This is SLOW!"
  END IF

  IF(numpe==1) PRINT *, "Allocated dynamic arrays used in main program"

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

  IF(numpe==1) PRINT *, "Computed steering array"

!------------------------------------------------------------------------------
! 6. Create interprocessor communication tables
!------------------------------------------------------------------------------

  CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)
 
  timest(8) = elap_time()

  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

  IF(numpe==1) PRINT *, "Created interprocessor communication tables"

!------------------------------------------------------------------------------
! 7. Read in loaded nodes and get starting r_pp
!------------------------------------------------------------------------------
 
  IF(loaded_nodes > 0) THEN

    ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes))
    ALLOCATE(r_pp(neq_pp))
  
    val  = zero
    r_pp = zero
    node = 0

    CALL read_loads(argv,numpe,node,val)
    CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))

    tload = SUM_P(r_pp(1:))

    DEALLOCATE(node,val)

  ELSE

    tload = zero

  END IF
  
  DEALLOCATE(g_g_pp)
  
  IF(numpe==1) PRINT *, "Read in loads and/or fixed displacements"

  timest(9) = elap_time()

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------

  dee = zero
  CALL deemat(dee,e,v)
  CALL sample(element,points,weights)
 
  IF(mult_method==1) THEN
    IF (nod==8) THEN
      ALLOCATE(km(ntot,ntot))
      ALLOCATE(iv(ntot,ntot))
      ALLOCATE(vox_storkm_pp(10,nels_pp)) 
      ALLOCATE(vtemp(ntot)) 
      ALLOCATE(ptemp(ntot)) 
      km            = zero
      vox_storkm_pp = zero
      vtemp         = zero
      iv            = 0
      iv(1,:)  = (/1,2,2,3,2,4,8,9,7,10,9,9,3,4,2,5,4,4,6,7,7,8,7,9/)
      iv(2,:)  = (/2,1,2,2,3,4,4,5,4,4,3,2,9,10,9,9,8,7,7,6,7,7,8,9/)
      iv(3,:)  = (/2,2,1,9,9,10,7,9,8,4,2,3,2,4,3,9,7,8,7,7,6,4,4,5/)
      iv(4,:)  = (/3,2,9,1,2,7,10,9,4,8,9,2,5,4,9,3,4,7,8,7,4,6,7,2/)
      iv(5,:)  = (/2,3,9,2,1,7,4,3,7,4,5,9,9,8,2,9,10,4,7,8,4,7,6,2/)
      iv(6,:)  = (/4,4,10,7,7,1,9,7,3,2,4,8,4,2,8,7,9,3,9,9,5,2,2,6/)
      iv(7,:)  = (/8,4,7,10,4,9,1,7,2,3,7,4,6,2,7,8,2,9,3,9,2,5,9,4/)
      iv(8,:)  = (/9,5,9,9,3,7,7,1,7,7,3,9,2,6,2,2,8,4,4,10,4,4,8,2/)
      iv(9,:)  = (/7,4,8,4,7,3,2,7,1,9,4,10,7,2,6,4,9,5,2,9,3,9,2,8/)
      iv(10,:) = (/10,4,4,8,4,2,3,7,9,1,7,7,8,2,4,6,2,2,5,9,9,3,9,7/)
      iv(11,:) = (/9,3,2,9,5,4,7,3,4,7,1,2,2,8,9,2,6,7,4,8,7,4,10,9/)
      iv(12,:) = (/9,2,3,2,9,8,4,9,10,7,2,1,9,4,5,2,7,6,4,7,8,7,4,3/)
      iv(13,:) = (/3,9,2,5,9,4,6,2,7,8,2,9,1,7,2,3,7,4,8,4,7,10,4,9/)
      iv(14,:) = (/4,10,4,4,8,2,2,6,2,2,8,4,7,1,7,7,3,9,9,5,9,9,3,7/)
      iv(15,:) = (/2,9,3,9,2,8,7,2,6,4,9,5,2,7,1,9,4,10,7,4,8,4,7,3/)
      iv(16,:) = (/5,9,9,3,9,7,8,2,4,6,2,2,3,7,9,1,7,7,10,4,4,8,4,2/)
      iv(17,:) = (/4,8,7,4,10,9,2,8,9,2,6,7,7,3,4,7,1,2,9,3,2,9,5,4/)
      iv(18,:) = (/4,7,8,7,4,3,9,4,5,2,7,6,4,9,10,7,2,1,9,2,3,2,9,8/)
      iv(19,:) = (/6,7,7,8,7,9,3,4,2,5,4,4,8,9,7,10,9,9,1,2,2,3,2,4/)
      iv(20,:) = (/7,6,7,7,8,9,9,10,9,9,8,7,4,5,4,4,3,2,2,1,2,2,3,4/)
      iv(21,:) = (/7,7,6,4,4,5,2,4,3,9,7,8,7,9,8,4,2,3,2,2,1,9,9,10/)
      iv(22,:) = (/8,7,4,6,7,2,5,4,9,3,4,7,10,9,4,8,9,2,3,2,9,1,2,7/)
      iv(23,:) = (/7,8,4,7,6,2,9,8,2,9,10,4,4,3,7,4,5,9,2,3,9,2,1,7/)
      iv(24,:) = (/9,9,5,2,2,6,4,2,8,7,9,3,9,7,3,2,4,8,4,4,10,7,7,1/)
    ELSE
      PRINT *, ""
      PRINT *, "The elements in the model have ",nod," nodes"
      PRINT *, "The voxel special case only works with 8 node bricks"
      PRINT *, "Analysis aborted"
      PRINT *, ""
      CALL shutdown()
    END IF
    vox_storkm_pp       = zero 
    DO iel=1,nels_pp
      km = zero
      DO i=1,nip
        CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
        det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
        CALL beemat(bee,deriv)
        km(1:ntot,1:ntot)=km(1:ntot,1:ntot) +                             &
                    MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)   
      END DO

      vox_storkm_pp(1,iel)  = km(1,1)  
      vox_storkm_pp(2,iel)  = km(1,2)  
      vox_storkm_pp(3,iel)  = km(1,4)  
      vox_storkm_pp(4,iel)  = km(1,6)  
      vox_storkm_pp(5,iel)  = km(1,16)  
      vox_storkm_pp(6,iel)  = km(1,19)  
      vox_storkm_pp(7,iel)  = km(1,9)  
      vox_storkm_pp(8,iel)  = km(1,7)  
      vox_storkm_pp(9,iel)  = km(1,8)  
      vox_storkm_pp(10,iel) = km(1,10)  

    END DO
  ELSE
    ALLOCATE(storkm_pp(ntot,ntot,nels_pp))
    storkm_pp=zero
    elements_3: DO iel=1,nels_pp
      gauss_pts_1: DO i=1,nip
        CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
        det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
        CALL beemat(bee,deriv)
        storkm_pp(1:ntot,1:ntot,iel)=storkm_pp(1:ntot,1:ntot,iel) +       &
        MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)   
      END DO gauss_pts_1
    END DO elements_3
  END IF

  timest(10) = elap_time()
  
  IF(numpe==1) PRINT *, "Calculated element stiffness arrays"

!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------
  
 ALLOCATE(diag_precon_pp(neq_pp))
 ALLOCATE(diag_precon_tmp(ntot,nels_pp))

 diag_precon_pp  = zero
 diag_precon_tmp = zero
 
 IF(mult_method==1) THEN
   DO iel = 1,nels_pp
     diag_precon_tmp(:,iel) = diag_precon_tmp(:,iel) + vox_storkm_pp(1,iel)
   END DO
 ELSE   
    elements_4: DO iel = 1,nels_pp 
     DO i = 1,ndof
      diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + storkm_pp(i,i,iel)
     END DO
    END DO elements_4   
 END IF
 
 CALL scatter(diag_precon_pp,diag_precon_tmp)

 DEALLOCATE(diag_precon_tmp)
 
 timest(11) = elap_time()

 IF(numpe==1) PRINT *, "Built the preconditioner"
 
!------------------------------------------------------------------------------
! 10. Read in fixed nodal displacements and assign to equations
!------------------------------------------------------------------------------

  IF(fixed_freedoms > 0) THEN

    ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),valf(fixed_freedoms),    &
             no_pp_temp(fixed_freedoms),sense(fixed_freedoms))

    node = 0 ; no = 0 ; no_pp_temp = 0 ; sense = 0 ; valf = zero

    CALL read_fixed(argv,numpe,node,sense,valf)
    CALL find_no(node,rest,sense,no)
!   CALL reindex_fixed_nodes(ieq_start,no,no_pp_temp,fixed_freedoms_pp,       &
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

  timest(12) = elap_time()

  IF(numpe==1) PRINT *, "Read in fixed displacements"

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
! 13. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------
  
  ALLOCATE(p_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),                         &
           u_pp(neq_pp),d_pp(neq_pp))

  p_pp    = zero  ;  x_pp = zero
  xnew_pp = zero  ;  u_pp = zero  ;  d_pp = zero

!------------------------------------------------------------------------------
! 14. Initiallize preconditioned conjugate gradient
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
  IF ( mult_method > 1 ) THEN

     timest(13) = elap_time()

     ! Flags for type of matmul based on method
     use_kernels = (mult_method==2 .OR. mult_method==3)
     use_amdblas = (mult_method==4 .OR. mult_method==5)

     print *, "Initializing GPU in process ", numpe
     misc_timers(1) = elap_time()
     
     ! Notation: vec_lhs = M . vec_rhs
     ! Single matrix and vector size
     matsize     = ntot**2
     vecsize_lhs = ntot
     vecsize_rhs = ntot

     ! Create an OpenCL context on the device
     status = init_opencl( gpu_device, numpe-1, use_amdblas )
     !status = init_opencl( gpu_device, -1, use_amdblas )
     if ( status /= 1 ) then
        print *, "Failed to init OpenCL"
        stop
     end if

     ! Set working sizes of device mem alloc according to matmul method
     IF ( use_kernels ) THEN
        ! Use own own kernel: try to transfer all arrays and vectors to device.
        msize     = nels_pp * matsize
        vsize_lhs = nels_pp * vecsize_lhs
        vsize_rhs = nels_pp * vecsize_rhs
     ELSE
        ! Use AMD clBlas methods. Can only transfer single matrices to device.
        msize     = matsize
        vsize_lhs = vecsize_lhs
        vsize_rhs = vecsize_rhs
     END IF

     ! Allocate read-only device memory (kernel won't write) for matrix data
     memflag = 0
     status = allocate_memory_on_gpu( msize, dsize, memflag, device_matrix )
     if ( status /= 1 ) then
        print *, "Failed to allocate device memory (matrix)"
        stop
     end if

     ! Allocate read-only device memory (kernel won't write) for rhs vector 
     memflag = 0
     status = allocate_memory_on_gpu( vsize_rhs, dsize, memflag, device_rhs_vectors)
     if ( status /= 1 ) then
        print *, "Failed to allocate device memory (rhs vector)"
        stop
     end if

     ! Allocate at least write (and possibly read) memory for lhs vector (the result vector)
     IF ( use_kernels ) THEN
        memflag = 1     ! Our own kernel only writes to the result vector
     ELSE
        memflag = 2     ! AMD clBlas can read and write the the result vector
     END IF
     status = allocate_memory_on_gpu( vsize_lhs, dsize, memflag, device_lhs_vectors)
     IF ( status /= 1 ) THEN
        print *, "Failed to allocate device memory (lhs vector)"
        stop
     END IF

     ! Compile our own kernels, upload all matrix data
     IF ( use_kernels ) THEN

        ! Use our own kernel (no build options used)
!        status = compile_kernel_from_file( srcfilename, kernelnames(mult_method), C_NULL_PTR )
!        IF ( status /= 1 ) THEN
!          print *, "Error compiling kernel ", kernelnames(mult_method)
!          stop
!        END IF

        ! Use our own kernel (no build options used)
        status = compile_kernel_from_file( srcfilename, kernelnames(mult_method-1), C_NULL_PTR )
        IF ( status /= 1 ) THEN
           print *, "Error compiling kernel ", kernelnames(mult_method-1)
           stop
        END IF

        timest(14) = elap_time()

        ! Our kernels want all data on the device so transfer it now
        status = copy_data_to_gpu(msize, dsize, storkm_pp, device_matrix )
        IF ( status /= 1 ) THEN
           print *, "Error copying entire storkm_pp array of matrices to the GPU."
           stop
        END IF

        timest(15) = elap_time()

     END IF

     misc_timers(2) = elap_time() - misc_timers(1)
     write(*,*) 'Time for GPU init + host-to-device copy: ', misc_timers(2)

  END IF
  
!------------------------------------------------------------------------------
! 14. Preconditioned conjugate gradient iterations
!------------------------------------------------------------------------------

  ALLOCATE(pmul_pp(ntot,nels_pp))
  ALLOCATE(utemp_pp(ntot,nels_pp))                      
  
  iters = 0

  iterations: DO 
    iters    = iters + 1
    u_pp     = zero
    pmul_pp  = zero
    utemp_pp = zero
    CALL gather(p_pp,pmul_pp)
    timest(16) = elap_time()
    misc_timers(3) = elap_time()
    
    IF(mult_method == 0) THEN
      ! Original cpu version    
      timest(16) = elap_time()
      elements_5: DO iel=1,nels_pp
        utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
      END DO elements_5
      timest(17) = timest(17) + (elap_time() - timest(16))
    ELSE IF(mult_method == 1) THEN
      timest(16) = elap_time()
!     DO iel=1,nels_pp
!       DO j=1,ntot
!         vtemp(i) = zero
!         DO i=1,ntot
!           k = iv(j,i)
!           vtemp(i) = vox_storkm_pp(k,iel)
!         END DO
!         utemp_pp(j,iel) = DOT_PRODUCT(vtemp,pmul_pp(:,iel))
!       END DO
!     END DO

      ! no memory copy
      ! --------------
      DO iel=1,nels_pp
        DO j=1,ntot
          utemp_pp(j,iel)=DOT_PRODUCT(vox_storkm_pp(iv(:,j),iel),pmul_pp(:,iel))
        END DO
      END DO

      timest(17) = timest(17) + (elap_time() - timest(16))
       ! Accumulate time for only the matmul code
       misc_timers(4) = misc_timers(4) + (elap_time() - misc_timers(3))
    ELSE IF(use_kernels) THEN
       ! gpu versions
       ! ------------

       ! Time GPU kernel and data upload/download
       misc_timers(1) = elap_time()
       timest(18) = elap_time()

          ! Transfer the current RHS vectors to device, do matmul, transfer result vectors back
          timest(18) = elap_time()
          status =          copy_data_to_gpu(vsize_rhs, dsize, pmul_pp, device_rhs_vectors )
          misc_timers(3) = elap_time()
          timest(19) = timest(19) + (elap_time() - timest(18))
          timest(20) = elap_time()

          status = status + matrix_vector_multiplies( nels_pp, vecsize_lhs, vecsize_rhs, &
                                                      device_matrix, device_rhs_vectors, device_lhs_vectors )
          misc_timers(4) = misc_timers(4) + (elap_time() - misc_timers(3))
          timest(21) = timest(21) + (elap_time() - timest(20))

          timest(22) = elap_time()
          status = status + copy_data_from_gpu(vsize_lhs, dsize, device_lhs_vectors, utemp_pp )
          timest(19) = timest(19) + (elap_time() - timest(22))

          IF ( status /= 3 ) THEN
             write(*,*) 'Error in OpenCL calls: ', status, ' of 3 succeeded'
             stop
          END IF
          

       ELSE IF ( mult_method == 4 ) THEN

          ! Copy a matrix/vector to GPU, use CL Blas, copy result vector back

          elements_5_clblas: DO iel=1,nels_pp
                
             ! Transfer current matrix and vector to GPU (using a blocking write)
             status =          copy_data_to_gpu(matsize,     dsize, storkm_pp(:,:,iel), device_matrix )
             status = status + copy_data_to_gpu(vecsize_rhs, dsize, pmul_pp(:,iel),     device_rhs_vectors )
                
             ! Multiply using AMD Blas
             misc_timers(3) = elap_time()
             status = status + blas_dgemv( vecsize_lhs, vecsize_rhs, device_matrix, device_rhs_vectors, device_lhs_vectors )
             misc_timers(4) = misc_timers(4) + (elap_time() - misc_timers(3))
             
             ! Transfer vector back from GPU (using a blocking read)
             status = status + copy_data_from_gpu(vecsize_lhs, dsize, device_lhs_vectors, utemp_pp(:,iel) )
                
             IF( status /= 4 ) THEN
                write(*,*) 'Error in OpenCL (clblas_sync) calls: ', status, ' of 4 succeeded'
                stop
             END IF
          END DO elements_5_clblas

       ELSE IF ( mult_method == 5 ) THEN
          
          ! The transfer matrices and vector asynchronously and use CL BLAS to multiply.
          status = mem_copy_blas_dgemv_async( nels_pp, vecsize_lhs, vecsize_rhs, storkm_pp, pmul_pp, utemp_pp, &
                                              device_matrix, device_rhs_vectors, device_lhs_vectors )
          IF ( status /= 1 ) THEN
             print *, 'Error in mem_copy_blas_dgemv_async'
             stop
          END IF
      
       ELSE
          write(*,*) 'Must set mult_method to 0,1,2,3,4,5'
          stop                
       END IF
       
       ! Only time the GPU work (roughly)
       misc_timers(2) = misc_timers(2) + (elap_time()-misc_timers(1))

     IF(mult_method > 1) THEN
       ! Compare GPU matmul results to CPU results
       IF( check_gpu ) THEN             
          elements_5_check_gpu: DO iel=1,nels_pp
             check_utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel), pmul_pp(:,iel))
             IF( c_compare_vecs(vecsize_lhs, check_utemp_pp(:,iel), utemp_pp(:,iel)) == 0 ) THEN
                write(*,*) 'Iteration ', iters, ' Matrix/vector ', iel, ' Fortran and C vectors NOT equal'
             END IF
          END DO elements_5_check_gpu
       END IF
       
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

  IF(numpe==1) PRINT *, "Solved equations"

  write(*,*) 'Total iters time CPU matmul() or GPU kernel: ', (misc_timers(4))
  write(*,*) 'Total iters time GPU kernel+transfer:        ', (misc_timers(2))
  
  IF(mult_method==1) THEN
   DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,vox_storkm_pp,pmul_pp)
  ELSE
   DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,storkm_pp,pmul_pp) 
  END IF

  ! Cleanup GPU
  IF ( mult_method > 1 ) THEN
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
  
  timest(23) = elap_time()
  
!------------------------------------------------------------------------------
! 15. Print out displacements, stress, principal stress and reactions
!------------------------------------------------------------------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  
  IF(numpe==1) THEN
    fname = argv(1:INDEX(argv, " ")-1)//".dis"
    OPEN(24, file=fname, status='replace', action='write')
    fname = argv(1:INDEX(argv, " ")-1) // ".str"
    OPEN(25, file=fname, status='replace', action='write')
    fname = argv(1:INDEX(argv, " ")-1) // ".pri"
    OPEN(26, file=fname, status='replace', action='write')
    fname = argv(1:INDEX(argv, " ")-1) // ".vms"
    OPEN(27, file=fname, status='replace', action='write')
    fname = argv(1:INDEX(argv, " ")-1) // ".rea"
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

 IF(write_stress == .false.) THEN

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

END IF

  timest(24) = elap_time()

!------------------------------------------------------------------------------
! 17. Output performance data
!------------------------------------------------------------------------------

  CALL WRITE_XX4(fixed_freedoms,iters,argv,loaded_nodes,mult_method,      &
                 neq,nn,npes,nr,numpe,timest,tload,vox_storkm) 

! CALL WRITE_XX4(fixed_freedoms,iters,argv,loaded_nodes,neq,nn,npes,nr,   &
!                numpe,timest,tload)

  CALL shutdown() 
 
END PROGRAM XX4
