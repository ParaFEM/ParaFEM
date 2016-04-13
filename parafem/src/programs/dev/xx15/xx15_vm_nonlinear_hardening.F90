PROGRAM xx15_vm_nonlinear_hardening
!------------------------------------------------------------------------------
!   program xx15:  finite strain elasto-plastic analysis with Newton-Raphson
!------------------------------------------------------------------------------

  ! Choice of solvers
  USE choose_solvers
  ! PETSc interface and modules.  PETSc will always use MPI (unless PETSc itself
  ! has been compiled without MPI)
  USE parafem_petsc
  USE PRECISION
  USE GLOBAL_VARIABLES
  USE MP_INTERFACE
  USE INPUT
  USE OUTPUT
  USE GATHER_SCATTER
  USE PARTITION
  USE MATHS
  USE TIMING
  USE LARGE_STRAIN
  USE plasticity_xx15
  
  IMPLICIT NONE

  ! PETSc types
#include <petsc/finclude/petscdef.h>
 
  INTEGER :: nels, nn, nr, nip, nodof=3, nod, nst=6, loaded_nodes, nn_pp,     &
   nf_start, fmt=1, i, j, k, l, ndim=3, iters, limit, iel, nn_start,          &
   num_load_steps, iload, igauss, dimH, inewton, jump, npes_pp, partitioner=1,&
   argc, iargc, limit_2, fixed_nodes, numfix_pp, fixdim, writetimes=0,        &
   nodes_pp, node_start, node_end, idx1, idx2, yield_ip_pp, unload_ip_pp,     &
   yield_tot, unload_tot

  REAL(iwp) :: det, tol, maxdiff, tol2, detF, detFinc, energy, energy1, rn0,  &
   timest(40), detF_mean, initial_guess, lambda, lambda_total, lambda_prev,   &
   energy_prev, energy_prev_prev, min_inc, next_output, temp_inc, dw, e, v,   &
   max_disp, max_disp_inc, det_max, det_min, sum_force

  REAL(iwp), PARAMETER :: tol_increment=0.000001_iwp, tol_val=0.00000001_iwp
  REAL(iwp), PARAMETER :: zero=0.0_iwp, one=1._iwp, half=0.5_iwp
  INTEGER, PARAMETER :: statevar_num=10
  
  REAL(iwp) :: E_mat(3,3), S_mat(3,3), C_mat(3,3), jacFtransinv(3,3),         &
   jacFtrans(3,3), sum_strain(3,3), sum_stress(3,3), sum_strain_pp(3,3),      &
   sum_stress_pp(3,3)
  
  INTEGER :: iter, prev_iter, max_inc, number_output, counter_output

  CHARACTER(len=15) :: element
  CHARACTER(len=50) :: text, fname_base, fname
 
  LOGICAL :: converged, timewrite=.TRUE., flag=.FALSE., print_output=.FALSE., &
   tol_inc=.FALSE., lambda_inc=.TRUE., flag_load=.TRUE., noncon_flag=.FALSE.

  REAL(iwp) :: sum_force_2_pp, sum_force_2

! Default solvers are ParaFEM.  parafem_solvers is defined in choose_solvers
  CHARACTER(len=50) :: solvers = parafem_solvers 

  ! PETSc variables
  CHARACTER(len=1024) :: p_fname
  LOGICAL             :: p_exist
  PetscErrorCode      :: p_ierr
  ! The PETSc objects cannot be initialised here because PETSC_NULL_OBJECT is a
  ! common-block-object and not a constant.
  Vec                 :: p_x,p_b,p_r
  Mat                 :: p_A
  KSP                 :: p_ksp
  PetscScalar         :: p_pr_n2,p_r_n2,p_b_n2
  PetscInt,ALLOCATABLE    :: p_rows(:),p_cols(:)
  PetscScalar,ALLOCATABLE :: p_values(:)
  PetscScalar,POINTER :: p_varray(:)
  PetscInt            :: p_its,p_nnz
  PetscReal           :: p_rtol
  DOUBLE PRECISION    :: p_info(MAT_INFO_SIZE)
  KSPConvergedReason  :: p_reason
  CHARACTER(LEN=p_max_string_length) :: p_description

  !-------------------------- dynamic arrays-----------------------------------
  REAL(iwp), ALLOCATABLE:: points(:,:), coord(:,:), weights(:), xnew_pp(:),   &
   diag_precon_pp(:), r_pp(:), bee(:,:), load_value(:,:), g_coord_pp(:,:),    &
   diag_precon_tmp(:,:), storekm_pp(:,:,:), disp(:,:), g_coord(:,:),          &
   val_pp(:), disp_pp(:,:), res_pp(:), upd_coord(:,:), fext_pp(:), auxm(:,:), &
   fextpiece_pp(:), deltax_pp(:), fint_pp(:), kmat_elem(:,:), jacF(:,:),      &
   xnewel_pp(:,:),   derivFtran(:,:), derivF(:,:), storefint_pp(:,:),         &
   deeF(:,:), jacFinv(:,:), beeF(:,:), sigma1C(:), fixed_value(:),            &
   fixval_pp(:), fixvalpiece_pp(:), jacFinc(:,:), statev(:), sigma(:,:),      &
   lnstrainelas_mat(:,:,:), lnstrainelas(:), statevar(:,:,:), xnewnodes_pp(:),&
   sigma1C_mat(:,:,:), xnewelinc_pp(:,:), auxm_inc(:,:), sigma1C_temp(:,:,:), &
   deltax_pp_temp(:), statevar_con(:,:,:), sigma1C_mat_con(:,:,:),            &
   lnstrainelas_mat_con(:,:,:), detFincmod_vec(:), geeF(:,:), geeFT(:,:),     &
   xnew_pp_previous(:), xnewel_pp_previous(:,:), auxm_previous(:,:),          &
   fextprev_pp(:), fixvalprev_pp(:), fixvaltot_pp(:), xnewprev_pp(:),         &
   value_shape(:), shape_integral_pp(:,:), stress_integral_pp(:,:),           &
   stressnodes_pp(:), strain_integral_pp(:,:), strainnodes_pp(:),             &
   principal_integral_pp(:,:), princinodes_pp(:), principal(:),               &
   reacnodes_pp(:), stiffness_mat_con(:,:,:,:)

  INTEGER, ALLOCATABLE  :: num(:), g_num(:,:), g_num_pp(:,:), g_g_pp(:,:),    &
   load_node(:), rest(:,:), nf_pp(:,:), no_pp(:), comp(:,:), fixed_node(:),   &
   fixed_dof(:), fixelem_pp(:), fixdof_pp(:), unload_pp(:,:)
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! 1. Input and initialisation
  !----------------------------------------------------------------------------

  timest(1) = ELAP_TIME()

  CALL FIND_PE_PROCS(numpe,npes)

  argc = IARGC()
!     Output:  argc i: Number of arguments in the command line,
!                      ./par131 arg1 arg2 arg3 ... argn

!  If I have forgotten to write the name of the file in the command
!  line, the program is stopped
  IF (argc < 1) THEN
    IF (numpe == 1) THEN
      WRITE(*,*) "Need name of filename_base!!"
    END IF
    CALL SHUTDOWN
    STOP
  END IF

! Input:  1: The first argument in the command line (arg1) 
  IF (argc >= 1) THEN
    CALL GETARG(1,fname_base)
  END IF

  fname = fname_base(1:INDEX(fname_base," ")-1) // ".dat"
  CALL READ_DATA_XX7(fname,numpe,nels,nn,nr,loaded_nodes,fixed_nodes,          &
                     nip,limit,tol,e,v,nod,num_load_steps,jump,tol2)

! Input:  2: The second argument in the command line (arg2)
  IF (argc >= 2) THEN
    CALL GETARG(2,solvers)
  END IF

  IF (.NOT. solvers_valid(solvers)) THEN
    IF (numpe == 1) THEN
      WRITE(*,*) "Solvers can be " // solvers_list()
    END IF
    CALL SHUTDOWN
    STOP
  END IF

  IF(solvers == petsc_solvers) THEN
    IF(numpe == 1) THEN
      p_fname = TRIM(fname_base) // ".petsc"
      INQUIRE(file=TRIM(p_fname),exist=p_exist)
      IF (.NOT. p_exist) THEN
        p_fname = ""
      END IF
    END IF
    CALL MPI_BCAST(p_fname,LEN(p_fname),MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  END IF
    
  IF (nels < npes) THEN
    IF (numpe==1) THEN
      WRITE(*,*)"Error: fewer elements than processors"
    END IF
    CALL SHUTDOWN
    STOP
  END IF

  IF(numpe==1) THEN
    fname = fname_base(1:INDEX(fname_base," ")-1) // ".res"
    OPEN (11, file=fname, status='replace', action='write')
  END IF

  IF (nod==8) THEN
    element='hexahedron'
    dimH=8
  ELSE IF (nod==4) THEN
    element='tetrahedron'
    dimH=4
  END IF
  
  ntot = nod * nodof

  CALL CALC_NODES_PP(nn,npes,numpe,node_end,node_start,nodes_pp)

!--------------------------------------------------------------------------
! 1a. Start up PETSc after MPI has been started
!--------------------------------------------------------------------------
  IF (solvers == petsc_solvers) THEN
    CALL PetscInitialize(p_fname,p_ierr)
    p_x   = PETSC_NULL_OBJECT
    p_b   = PETSC_NULL_OBJECT
    p_r   = PETSC_NULL_OBJECT
    p_A   = PETSC_NULL_OBJECT
    p_ksp = PETSC_NULL_OBJECT
  END IF

!------------------------------------------------------------------------------
! 2. Get integration Gauss points and weights in the element
!------------------------------------------------------------------------------
  
  ALLOCATE(points(ndim,nip), weights(nip))
  CALL GET_GAUSS_POINTS(element,points,weights)

!------------------------------------------------------------------------------
! 3. Import and distribute mesh
!------------------------------------------------------------------------------

  CALL CALC_NELS_PP(fname_base,nels,npes,numpe,partitioner,nels_pp)
  
  timest(2) = ELAP_TIME()

  ALLOCATE(g_num_pp(nod, nels_pp)) 
  
  fname = fname_base(1:INDEX(fname_base," ")-1) // ".d"
  CALL READ_ELEMENTS_2(fname,npes,nn,numpe,g_num_pp)
  
  timest(3) = ELAP_TIME()

  CALL CALC_NN_PP(g_num_pp,nn_pp,nn_start)
  
  timest(4) = ELAP_TIME()

  ALLOCATE(g_coord_pp(ndim, nn_pp))
  CALL READ_NODES(fname,nn,nn_start,numpe,g_coord_pp)
  
  timest(5) = ELAP_TIME()
  
  ALLOCATE(rest(nr,nodof+1))
  fname = fname_base(1:INDEX(fname_base," ")-1) // ".bnd"
  CALL READ_RESTRAINTS(fname,numpe,rest)
  
  timest(6) = ELAP_TIME()

!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in the main program  
!------------------------------------------------------------------------------

  ALLOCATE(coord(nod,ndim))
  ALLOCATE(upd_coord(nod,ndim))
  ALLOCATE(bee(nst,ntot))
  ALLOCATE(num(nod))
  ALLOCATE(g_g_pp(ntot,nels_pp))
  ALLOCATE(load_value(ndim,loaded_nodes))
  ALLOCATE(load_node(loaded_nodes))
  ALLOCATE(nf_pp(nodof,nn_pp))
  ALLOCATE(storekm_pp(ntot,ntot,nels_pp))
  ALLOCATE(kmat_elem(ntot,ntot))
  ALLOCATE(xnewel_pp(ntot,nels_pp))
  ALLOCATE(xnewel_pp_previous(ntot,nels_pp))
  ALLOCATE(comp(nod,ndim))
  ALLOCATE(jacF(ndim,ndim))
  ALLOCATE(jacFinc(ndim,ndim))
  ALLOCATE(lnstrainelas_mat_con(nels_pp,nip,nst))
  ALLOCATE(lnstrainelas(nst))
  ALLOCATE(statev(statevar_num))
  ALLOCATE(statevar_con(nels_pp,nip,statevar_num))
  ALLOCATE(sigma(ndim,ndim))
  ALLOCATE(sigma1C_mat_con(nels_pp,nip,nst))
  ALLOCATE(auxm(nod,ndim))
  ALLOCATE(jacFinv(ndim,ndim))
  ALLOCATE(derivFtran(nod,ndim))
  ALLOCATE(derivF(ndim,nod))
  ALLOCATE(beeF(nst,ntot))
  ALLOCATE(geeF((ndim*ndim),ntot))
  ALLOCATE(geeFT(ntot,(ndim*ndim)))
  ALLOCATE(sigma1C(nst))
  ALLOCATE(storefint_pp(ntot,nels_pp))
  ALLOCATE(deeF((ndim*ndim),(ndim*ndim)))
  ALLOCATE(principal(ndim))
  ALLOCATE(value_shape(nod))
  ALLOCATE(xnewelinc_pp(ntot,nels_pp))
  ALLOCATE(auxm_inc(nod,ndim))
  ALLOCATE(auxm_previous(nod,ndim))
  ALLOCATE(stiffness_mat_con(nels_pp,nip,(ndim*ndim),(ndim*ndim)))
  ALLOCATE(unload_pp(nels_pp,nip))

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of 
!    equations to solve. 
!------------------------------------------------------------------------------

  CALL REST_NF(nn_start,rest,nf_pp)

  DO iel = 1,nels_pp
    CALL NUM_TO_G2(g_num_pp(:,iel), nf_pp, g_g_pp(:,iel), nn_start)
  END DO

  CALL CALC_NEQ(nn,rest,neq)
  
  timest(7) = ELAP_TIME()

!------------------------------------------------------------------------------
! 6. Create interprocessor communications tables
!------------------------------------------------------------------------------  
  
  CALL CALC_NPES_PP(npes,npes_pp)

  CALL CALC_NEQ_PP

  CALL MAKE_GGL(npes_pp,npes,g_g_pp)

  timest(8) = ELAP_TIME()
  
!------------------------------------------------------------------------------
! 7. Read and distribute essential boundary conditions
!------------------------------------------------------------------------------

  numfix_pp = 0

  IF (fixed_nodes>0) THEN

    ALLOCATE(fixed_node(fixed_nodes))
    ALLOCATE(fixed_dof(fixed_nodes))
    ALLOCATE(fixed_value(fixed_nodes))

    CALL READ_FIXED(fname_base,numpe,fixed_node,fixed_dof,fixed_value)

    IF (element=='hexahedron') THEN
      fixdim=4
    ELSE IF (element=='tetrahedron') THEN
      fixdim=20
    END IF

    ALLOCATE(fixelem_pp(fixdim*fixed_nodes))
    ALLOCATE(fixdof_pp(fixdim*fixed_nodes))
    ALLOCATE(fixval_pp(fixdim*fixed_nodes))
    ALLOCATE(fixvalpiece_pp(fixdim*fixed_nodes))
    ALLOCATE(fixvalprev_pp(fixdim*fixed_nodes))
    ALLOCATE(fixvaltot_pp(fixdim*fixed_nodes))
    
    fixelem_pp      = 0
    fixdof_pp       = 0
    fixval_pp       = zero
    fixvalpiece_pp  = zero
    fixvalprev_pp   = zero
    fixvaltot_pp    = zero
    
    DO i = 1,fixed_nodes
      DO iel = 1,nels_pp
        DO j = 1,nod
          IF (fixed_node(i)==g_num_pp(j,iel)) THEN
            numfix_pp = numfix_pp + 1
            fixelem_pp(numfix_pp) = iel
            fixdof_pp(numfix_pp) = (j-1)*ndim + fixed_dof(i)
            fixval_pp(numfix_pp) = fixed_value(i)
          END IF
        END DO
      END DO
    END DO
    
    max_disp = MAXVAL(ABS(fixed_value))

    DEALLOCATE(fixed_dof, fixed_value)

  END IF
  
  timest(9) = ELAP_TIME()

  !----------------------------------------------------------------------------
  ! 8. Read and distribute natural boundary conditions
  !----------------------------------------------------------------------------

  ALLOCATE(fextpiece_pp(0:neq_pp))
  ALLOCATE(fext_pp(0:neq_pp))
  ALLOCATE(fextprev_pp(0:neq_pp))
  
  fextpiece_pp = zero
  fext_pp      = zero
  fextprev_pp  = zero

  IF (loaded_nodes>0) THEN

    CALL READ_LOADS(fname_base,numpe,load_node,load_value)

    CALL LOAD_2(nn_start,g_num_pp,load_node,load_value,nf_pp,fext_pp(1:))
   
    DEALLOCATE(load_node,load_value)

  END IF
  
  timest(10) = ELAP_TIME()

  !----------------------------------------------------------------------------
  ! 9. Allocate arrays dimensioned by neq_pp
  !----------------------------------------------------------------------------

  ALLOCATE(r_pp(0:neq_pp), xnew_pp(0:neq_pp), diag_precon_pp(0:neq_pp))
  ALLOCATE(res_pp(0:neq_pp), deltax_pp(0:neq_pp), fint_pp(0:neq_pp)) 
  ALLOCATE(deltax_pp_temp(0:neq_pp), xnew_pp_previous(0:neq_pp))
  ALLOCATE(xnewprev_pp(0:neq_pp))

  r_pp             = zero
  xnew_pp          = zero
  diag_precon_pp   = zero
  res_pp           = zero
  deltax_pp        = zero
  fint_pp          = zero
  xnew_pp_previous = zero
  xnewprev_pp      = zero

  ALLOCATE(diag_precon_tmp(ntot,nels_pp))

!---------------------------------------------------------------------------
! 9a. Set up PETSc
!---------------------------------------------------------------------------
  IF (solvers == petsc_solvers) THEN
    ! PETSc 64-bit indices and 64-bit reals.  In most (all?) places, passing a
    ! 32-bit integer where an intent(in) 64-integer is required is safe, because
    ! the PETSc Fortran-C interface de-references the pointer that is actually
    ! passed, even though it does not specify any intent in the Fortran
    ! interface.  This is not safe when passing arrays, they need to be copied
    ! to PetscInt arrays.  And for safety the same should be done for the
    ! PetscScalar arrays.
    CALL MatCreate(PETSC_COMM_WORLD,p_A,p_ierr)
    CALL MatSetSizes(p_A,neq_pp,neq_pp,PETSC_DETERMINE,PETSC_DETERMINE,p_ierr)
    CALL MatSetType(p_A,MATAIJ,p_ierr)
    !- Block size fixed to 1 for just now - this cannot be set to nodof until
    !- the restraints are handled block-wise in ParaFEM.  CALL
    !- MatSetBlockSize(p_A,nodof,p_ierr)
    
    ! Find an approximate number of zeroes per row for the matrix size
    ! pre-allocation.
    CALL p_row_nnz(nodof,ndim,nod,p_over_allocation,p_nnz)
    CALL MatSeqAIJSetPreallocation(p_A,p_nnz,PETSC_NULL_INTEGER,p_ierr)
    CALL MatMPIAIJSetPreallocation(p_A,p_nnz,PETSC_NULL_INTEGER,              &
      p_nnz,PETSC_NULL_INTEGER,p_ierr)
    ! If the allocation is too small, PETSc will produce reams of information
    ! and not construct the matrix properly.  We output some information at the
    ! end if p_over_allocation should be increased.
    CALL MatSetOption(p_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,p_ierr)
    
    ! RHS vector.  For this particular order (create, set size, set type) the
    ! allocation is done by set type.
    CALL VecCreate(PETSC_COMM_WORLD,p_b,p_ierr)
    CALL VecSetSizes(p_b,neq_pp,PETSC_DECIDE,p_ierr)
    CALL VecSetType(p_b,VECSTANDARD,p_ierr)
    !- Block size fixed to 1 for just now - this cannot be set to nodof until the
    !- restraints are handled block-wise in ParaFEM.
    !- CALL VecSetBlockSize(p_b,nodof,p_ierr)

    ! Solution vector
    CALL VecDuplicate(p_b,p_x,p_ierr) 
    
    ! Krylov solver data structures
    CALL KSPCreate(PETSC_COMM_WORLD,p_ksp,p_ierr)
    CALL KSPSetOperators(p_ksp,p_A,p_A,p_ierr)
  
    ! Arrays of PETSc types to be proof against changes in index and scalar
    ! sizes.
    ALLOCATE(p_rows(ntot),p_cols(ntot),p_values(ntot*ntot))
  END IF

  !----------------------------------------------------------------------------
  ! 10. Initialise the solution vector to 0.0
  !----------------------------------------------------------------------------

  ! Vector comp to compute F (gradient of deformation)
  DO i = 1,nod
    DO j = 1,ndim
      comp(i,j) = (i-1)*ndim + j
    END DO
  END DO

  ! Initialise the converged stresses, elastic strains and the state variables
  sigma1C_mat_con      = zero
  lnstrainelas_mat_con = zero
  statevar_con         = zero
  stiffness_mat_con    = zero

  ! Limit for Newton-Raphson iterations
  limit_2=10
  
  ! Limit for the PCG iteration
  limit=12500

  ! Establish a maximum number of increments
  max_inc=100

  ! Establish the minimum time increment
  min_inc=0.008_iwp

  ! Establish the number of times the output is printed
  number_output=1
  counter_output=1
  
  ! Set the initial guess (normalized)
  initial_guess  = 0.01_iwp
  lambda         = initial_guess
  lambda_total   = zero
  lambda_prev    = zero
  next_output=one/number_output
  
  ! Initiate counters
  unload_ip_pp=0
  yield_ip_pp=0
  
  timest(11) = ELAP_TIME()

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!----------------------- Start Load LOOP --------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  sum_force=zero

  DO iload = 1,max_inc
    converged = .FALSE.

    ! Increase the stage control by 50% if the number of iterations of the two 
    ! previously converged full increments (not affected by the output  
    ! requests) are below 6
    IF ((iter<6).AND.(prev_iter<6).AND.(iload>2)) THEN
      IF(numpe==1) THEN
        WRITE(*,*) 'The load increment is increased by 50%'
      END IF
      lambda=1.5*lambda
    END IF

    ! Make sure it is not larger than unity
    !IF (counter_output>number_output) THEN
    !  GOTO 200
    !END IF

    100 CONTINUE

    IF (lambda_total>=(1-tol_increment)) THEN
    !IF (lambda_total>=(0.1_iwp-tol_increment)) THEN
      EXIT
    END IF

    
    
    
    
    
    
    !IF ((iload==9).AND.(flag_load)) THEN
    !  lambda=0.01_iwp
    !  flag_load=.FALSE.
    !END IF
    
    
    
    
    
    
    
    
    
    ! Update the load increments
    lambda_prev=lambda_total
    lambda_total=lambda_total+lambda

    !IF (tol_inc) THEN
    !  lambda=0.5*min_inc
    !  tol_inc=.FALSE.
    !END IF

    ! Make sure that we print the corresponding output
    !IF (lambda_total>(next_output-tol_increment)) THEN
    !  WRITE(*,*) 'The load increment exceeds the next output request time ',  &
    !   next_output
    !  WRITE(*,*) 'Adjusting the load increment by ',(lambda-(lambda-          &
    !   (lambda_total-next_output)))
    !  WRITE(*,*) 'This output request is the number ',counter_output
    !  print_output=.TRUE.
    !  lambda_inc=.FALSE.
    !  prev_inc=lambda
    !  lambda_total=lambda_total-lambda
    !  lambda=next_output-lambda_prev
    !  lambda_total=lambda_total+lambda
    !  next_output=FLOAT(counter_output+1)/number_output
    !  counter_output=counter_output+1
    !  IF (lambda<min_inc) THEN
    !    tol_inc=.TRUE.
    !  END IF
    !END IF

    ! Exit if the time increment is less that the specified minimum
    IF (lambda<min_inc) THEN
      IF(numpe==1) THEN
        WRITE(*,*) 'The load increment is too small'
      END IF
      EXIT
    END IF

    ! Display the incremental, total and previous load increments
    IF(numpe==1) THEN
      WRITE(*,*) 'For the increment ',iload,' the data is as following:'
      WRITE(*,*) 'The current load increment is ',lambda
      WRITE(*,*) 'The previous total load multiplier is ',lambda_prev
      WRITE(*,*) 'The total load multiplier is ',lambda_total
    END IF

    ! Calculate the displacement increment, total displacement, and previously 
    ! converged displacement
    fixvalpiece_pp(1:)=fixval_pp(1:)*lambda
    fixvaltot_pp(1:)=fixval_pp(1:)*lambda_total
    fixvalprev_pp(1:)=fixval_pp(1:)*lambda_prev

    !IF (print_output) THEN
    !  !lambda=(1.0/number_output)
    !  temp_inc=lambda
    !  lambda=prev_inc
    !END IF

!------------------------------------------------------------------------------
!----------------------- Start Newton-Raphson iterations ----------------------
!------------------------------------------------------------------------------

    deltax_pp_temp = zero
    inewton = 0
    
    iterations: DO
      inewton = inewton + 1

      storefint_pp = zero
      
      CALL GATHER(xnew_pp(1:),xnewel_pp)
      CALL GATHER(deltax_pp_temp(1:),xnewelinc_pp)
      CALL GATHER(xnew_pp_previous(1:),xnewel_pp_previous)

      DO i = 1,numfix_pp
        xnewel_pp_previous(fixdof_pp(i),fixelem_pp(i))=fixvalprev_pp(i)
      END DO
        
      IF (inewton>1 .AND. numfix_pp>0) THEN
        DO i = 1,numfix_pp
          xnewel_pp(fixdof_pp(i),fixelem_pp(i)) = fixvaltot_pp(i)
        END DO
        DO i = 1,numfix_pp
          xnewelinc_pp(fixdof_pp(i),fixelem_pp(i)) = fixvalpiece_pp(i)
        END DO
      END IF

      IF (iload>1 .AND. inewton==1 .AND. numfix_pp>0) THEN
        DO i = 1,numfix_pp
          xnewel_pp(fixdof_pp(i),fixelem_pp(i))=fixvalprev_pp(i)
        END DO
      END IF      

!------------------------------------------------------------------------------
! 11. Element stiffness integration and storage
!------------------------------------------------------------------------------
      
      timest(12) = ELAP_TIME()
      
      storekm_pp = zero
      IF (solvers == petsc_solvers) THEN
        CALL MatZeroEntries(p_A,p_ierr)
      END IF
      
      DO iel = 1,nels_pp

        DO i = 1,nod
          num(i) = g_num_pp(i,iel) - nn_start + 1
        END DO

        DO i = 1,ndim
          auxm(:,i) = xnewel_pp(comp(:,i),iel)
        END DO

        DO i = 1,ndim
          auxm_inc(:,i) = xnewelinc_pp(comp(:,i),iel)
        END DO

        DO i = 1,ndim
          auxm_previous(:,i) = xnewel_pp_previous(comp(:,i),iel)
        END DO 
        
        coord = TRANSPOSE(g_coord_pp(:,num))
        upd_coord=coord+auxm_previous
        
        timest(13) = ELAP_TIME()

        DO igauss = 1,nip

          ! Initialise the state variables to the same converged value of     
          ! the last time increment
          statev(:)=statevar_con(iel,igauss,:)
          lnstrainelas(:)=lnstrainelas_mat_con(iel,igauss,:)

          ! Calculates the total and incremental deformation gradients
          CALL DEFGRA(igauss,auxm,coord,points,det,detF,beeF,geeF,jacF,ndim,  &
           nod)
          CALL DEFGRAINC(igauss,auxm_inc,upd_coord,points,jacFinc,ndim,nod)
          
          timest(14) = ELAP_TIME()

          
          
          
          IF ((iel==1).AND.(igauss==1)) THEN
            det_max=det
            det_min=det
          END IF
          
          IF (det>=det_max) THEN
            det_max=det
          END IF
          
          IF (det<=det_min) THEN
            det_min=det
          END IF
          

          
          
          CALL PLASTICITY(deeF,jacF,jacFinc,sigma1C,statev,lnstrainelas,sigma, &
                          detF,statevar_num,iel,igauss,noncon_flag,            &
                          umat_vm_nonlinear_hardening)
           
          timest(15) = ELAP_TIME()

          ! During the first Newton-Raphson iteration, retrieve the previous  
          ! converged stress and stiffness tensor
          IF ((iload>1).AND.(inewton==1)) THEN
            deeF(:,:)=stiffness_mat_con(iel,igauss,:,:)
            sigma1C(:)=sigma1C_mat_con(iel,igauss,:)
          END IF
          
          dw = det*weights(igauss)
          
          ! Calculate the internal force vector of the element
          storefint_pp(:,iel)=storefint_pp(:,iel) + MATMUL(TRANSPOSE(beeF),   &
           sigma1C)*dw
           
          ! Calculate the stiffness tensor of the element
          geeFT = TRANSPOSE(geeF)
          storekm_pp(:,:,iel)=storekm_pp(:,:,iel) + (MATMUL(MATMUL(geeFT,     &
           deeF),geeF)*dw)

        END DO

        IF (solvers == petsc_solvers) THEN
          ! 1/ Use ubound(g_g_pp,1) instead of ntot? 2/ The following depends on
          ! g_g_pp holding 0 for erased rows/columns, and MatSetValues ignoring
          ! negative indices (PETSc always uses zero-based indexing). 3/ PETSc
          ! uses C array order, so a transpose is needed.
          p_rows   = g_g_pp(:,iel) - 1
          p_cols   = p_rows
          p_values = RESHAPE(TRANSPOSE(storekm_pp(:,:,iel)),SHAPE(p_values))
          CALL MatSetValues(p_A,ntot,p_rows,ntot,p_cols,p_values,ADD_VALUES,   &
                            p_ierr)
        END IF
      END DO

      IF (solvers == petsc_solvers) THEN
        CALL MatAssemblyBegin(p_A,MAT_FINAL_ASSEMBLY,p_ierr)
        CALL MatAssemblyEnd(p_A,MAT_FINAL_ASSEMBLY,p_ierr)
        
        CALL MatGetInfo(p_A,MAT_GLOBAL_SUM,p_info,p_ierr)
        IF(numpe==1)THEN
          IF(p_info(MAT_INFO_MALLOCS)/=0.0)THEN
            WRITE(*,'(A,I0,A)') "The matrix assembly required ",               &
              NINT(p_info(MAT_INFO_MALLOCS)),                                  &
              " mallocs.  Increase p_over_allocation to speed up "             &
              //"the assembly."
          END IF
        END IF
      END IF
      
      ! Time is accumulated over all load increments
      timest(16) = timest(16) + (ELAP_TIME() - timest(12))
      
!------------------------------------------------------------------------------
! 12. Build and invert the preconditioner (ParaFEM only)
!------------------------------------------------------------------------------

      IF (solvers == parafem_solvers) THEN
        timest(17) = ELAP_TIME()
        
        diag_precon_tmp = zero
        
        DO iel = 1,nels_pp
          DO k = 1,ntot 
            diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel) + storekm_pp(k,k,iel)
          END DO
        END DO
        
        diag_precon_pp(:) = zero
        CALL SCATTER(diag_precon_pp(1:),diag_precon_tmp)
        
        diag_precon_pp(1:) = one/diag_precon_pp(1:)
        diag_precon_pp(0)  = zero
      END IF

!------------------------------------------------------------------------------
! 13. Initialize PCG
!------------------------------------------------------------------------------
      
      timest(18) = ELAP_TIME()
      
      ! During the first Newton-Raphson iteration, the incremental 
      ! displacements are applied through a linear mapping of these 
      ! displacements to the internal force vector
      IF (inewton==1 .AND. numfix_pp>0) THEN
        DO i = 1,numfix_pp
          DO j = 1,ntot
            IF (g_g_pp(j,fixelem_pp(i))>0) THEN
              storefint_pp(j,fixelem_pp(i)) = storefint_pp(j,fixelem_pp(i)) + &
                fixvalpiece_pp(i)*storekm_pp(j,fixdof_pp(i),fixelem_pp(i))
            END IF
          END DO
        END DO
      END IF
      
      fint_pp = zero
      CALL SCATTER(fint_pp(1:),storefint_pp)

      fextpiece_pp(1:) = lambda_total*fext_pp(1:)
      r_pp(1:) = fextpiece_pp(1:) - fint_pp(1:)       !Residual
      r_pp(0) = zero

      ! Compute maxdiff of residual 
      maxdiff =  MAXABSVAL_P(r_pp(1:))

      ! Normalise residual vector and stiffness matrix for pcg
      IF (maxdiff == zero) THEN
        IF(numpe==1) THEN
          WRITE(*,*) "maxdiff = zero and now exiting loop"
        END IF
        EXIT
      END IF

      IF (solvers == petsc_solvers) THEN
        ! The general tolerance and maximum number of linear solver iterations
        ! come from the value in the ParaFEM control file.  Note that relative
        ! tolerance for PETSc is for the preconditioned residual.
        CALL KSPSetTolerances(p_ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, &
                              limit,p_ierr)
        ! KSP type, per-KSP tolerances (rtol, abstol, dtol, maxits), KSP
        ! options, PC type, PC options are set in the xx*.petsc file.  Those
        ! options are used to set up the preconditioned Krylov solver.  If there
        ! are several KSP types to be chosen from, then each one will be
        ! bracketed by -prefix_push and -prefix_pop.  For example
        !
        ! -prefix_push abc1_
        !   -ksp_type minres
        ! -prefix_pop
        !
        ! in the xx*.petsc file and 
        !
        ! CALL KSPSetOptionsPrefix(p_ksp,"abc1_",p_ierr)
        !
        ! before KSPSetFromOptions before using the 'abc1_' solver.  Thus you
        ! can switch for CG to GMRES during a simulation.
        CALL KSPSetFromOptions(p_ksp,p_ierr)

        ! load vector
        CALL VecGetArrayF90(p_b,p_varray,p_ierr)
        ! This is OK as long as PetscScalars not smaller than ParaFEM reals.  There
        ! should be a test for sizes of PetscScalars (and PetscInts) and ParaFEM reals
        ! (and indices).
        p_varray = r_pp(1:)
        CALL VecRestoreArrayF90(p_b,p_varray,p_ierr)
      END IF

!------------------------------------------------------------------------------
!----------------- Solve using preconditioned Krylov solver -------------------
!------------------------------------------------------------------------------
      
      deltax_pp = zero
      res_pp    = r_pp

      IF (solvers == parafem_solvers) THEN
        CALL PCG_VER1(inewton,limit,tol,storekm_pp,r_pp(1:),                   &
                      diag_precon_pp(1:),rn0,deltax_pp(1:),iters)
      ELSE IF (solvers == petsc_solvers) THEN
        ! Solution vector
        ! For non-linear solves, the previous solution will be used as an initial
        ! guess: copy the ParaFEM solution vector to PETSc.
        CALL VecGetArrayF90(p_x,p_varray,p_ierr)
        ! This is OK as long as PetscScalars are not smaller than ParaFEM reals.
        ! There should be a test for sizes of PetscScalars (and PetscInts) and ParaFEM
        ! reals (and indices).
        p_varray = deltax_pp(1:)
        CALL VecRestoreArrayF90(p_x,p_varray,p_ierr)
        CALL KSPSetInitialGuessNonzero(p_ksp,PETSC_TRUE,p_ierr)
        CALL KSPSolve(p_ksp,p_b,p_x,p_ierr)

        ! Preconditioned norm tolerance
        CALL KSPGetTolerances(p_ksp,p_rtol,PETSC_NULL_REAL,PETSC_NULL_REAL,    &
                              PETSC_NULL_INTEGER,p_ierr)
        ! True residual L2 norm
        CALL VecDuplicate(p_b,p_r,p_ierr)
        CALL KSPBuildResidual(p_ksp,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,p_r,p_ierr)
        CALL VecNorm(p_r,NORM_2,p_r_n2,p_ierr)
        CALL VecDestroy(p_r,p_ierr)
        ! L2 norm of load
        CALL VecNorm(p_b,NORM_2,p_b_n2,p_ierr)
        
        CALL KSPGetIterationNumber(p_ksp,p_its,p_ierr)
        CALL KSPGetConvergedReason(p_ksp,p_reason,p_ierr)
        CALL p_describe_reason(p_reason,p_description)
        
        ! Copy PETSc solution vector to ParaFEM
        CALL VecGetArrayF90(p_x,p_varray,p_ierr)
        ! This is OK as long as ParaFEM reals not smaller than PetscScalars.  There
        ! should be a test for sizes of PetscScalars (and PetscInts) and ParaFEM reals
        ! (and indices).
        deltax_pp(1:) = p_varray
        CALL VecRestoreArrayF90(p_x,p_varray,p_ierr)

        IF(numpe == 1)THEN
          WRITE(11,'(A,I0,A)')                                                 &
            "The reason for convergence was ",p_reason," "//TRIM(p_description)
          WRITE(11,'(A,I0)')                                                   &
            "The number of iterations to convergence was ",p_its
          WRITE(11,'(A,E17.7)')                                                &
            "The preconditioned relative error tolerance was ",p_rtol
          WRITE(11,'(A,E17.7)')                                                &
            "The relative error ||b-Ax||/||b|| was           ",p_r_n2/p_b_n2
        END IF
      END IF

      IF (numpe==1) THEN
        WRITE(91,*)iload,inewton,iters
        CALL FLUSH(91)
      END IF

      timest(19) = timest(19) + (ELAP_TIME() - timest(18))

      ! Total displacements
      xnew_pp(1:) = xnew_pp(1:) + deltax_pp(1:)
      xnew_pp(0) = zero
      
      ! Incremental displacements corresponding to the time increment
      deltax_pp_temp(1:) = deltax_pp_temp(1:) + deltax_pp(1:)
      deltax_pp_temp(0) = zero

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

! Check convergence for Newton-Raphson iterations 
      energy_prev_prev=energy_prev
      energy_prev=energy
      energy = ABS(DOT_PRODUCT_P(res_pp(1:),deltax_pp(1:)))
      IF (inewton==1) THEN
        energy1 = energy
      END IF
      IF (numpe==1) THEN
        WRITE(*,2000)iload,inewton,energy,energy/energy1
      END IF
      IF (inewton>1) THEN
        IF ((energy/energy1)<=tol2) THEN
          converged = .TRUE.
        END IF 
      END IF

      ! The time increment is cut in half if one of the following situations
      ! happen:
      ! - The maximum number of iterations has been met
      ! - The current iteration diverged (the current energy is larger than
      ! the two previous energies)
      ! - There is some numerical instability which causes some NaN values
      IF ((inewton==limit_2).OR.(((energy_prev_prev/energy1)<(energy_prev/    &
       energy1)).AND.((energy_prev/energy1)<(energy/energy1)).AND.            &
        (inewton>2)).OR.(ISNAN(energy))) THEN
        
        IF (numpe==1) THEN
          WRITE(*,*) 'The load increment is cut in half'
        END IF
        !IF (print_output) THEN
        !  lambda_total=lambda_total-temp_inc
        !ELSE
        lambda_total=lambda_total-lambda
        !END IF
        
      
        lambda=0.25_iwp*lambda
        

        
        xnew_pp=xnewprev_pp
        flag=.FALSE.
        IF (print_output) THEN
          counter_output=counter_output-1
          next_output=FLOAT(counter_output)/number_output
          print_output=.FALSE.
        END IF
        
        ! Set the number of iterations to a number bigger than six to reset the  
        ! iteration counter
        prev_iter=7
        iter=7
        GOTO 100
      END IF
 
      ! After convergence, a last "iteration" is needed to calculate the fully
      ! converged values (some FE algorithms omit this step, but we perform
      ! it because we are interested in a fully precise solution)
      IF (converged) THEN 

        timest(20) = ELAP_TIME()
      
        ! If the current increment is not an output request increment, save the  
        ! current and previous increment number of iterations to allow for   
        ! for lambda increment to be increased
        IF (.NOT.(lambda_inc)) THEN
          IF (iload>1) THEN
            prev_iter=iter
          END IF
          iter=inewton
          lambda_inc=.TRUE.
        END IF

        xnewprev_pp=xnew_pp

        CALL GATHER(xnew_pp(1:),xnewel_pp)
        CALL GATHER(deltax_pp_temp(1:),xnewelinc_pp)
        CALL GATHER(xnew_pp_previous(1:),xnewel_pp_previous)

        DO i = 1,numfix_pp
          xnewel_pp_previous(fixdof_pp(i),fixelem_pp(i))=fixvalprev_pp(i)
        END DO
        
        DO i = 1,numfix_pp
          xnewel_pp(fixdof_pp(i),fixelem_pp(i)) = fixvaltot_pp(i)
        END DO

        DO i = 1,numfix_pp
          xnewelinc_pp(fixdof_pp(i),fixelem_pp(i)) = fixvalpiece_pp(i)
        END DO

        ! Homogenized strain and stress start
        sum_stress_pp = zero
        sum_strain_pp = zero
        
        DO iel = 1,nels_pp

          DO i = 1,nod
            num(i) = g_num_pp(i,iel) - nn_start + 1
          END DO

          DO i = 1,ndim
            auxm(:,i) = xnewel_pp(comp(:,i),iel)
          END DO

          DO i = 1,ndim
            auxm_inc(:,i) = xnewelinc_pp(comp(:,i),iel)
          END DO

          DO i = 1,ndim
            auxm_previous(:,i) = xnewel_pp_previous(comp(:,i),iel)
          END DO 
          
          coord = TRANSPOSE(g_coord_pp(:,num))
          upd_coord=coord+auxm_previous

          timest(21) = ELAP_TIME()

          DO igauss = 1,nip

            statev(:)=statevar_con(iel,igauss,:)
            lnstrainelas(:)=lnstrainelas_mat_con(iel,igauss,:)

            CALL DEFGRA(igauss,auxm,coord,points,det,detF,beeF,geeF,jacF,ndim,&
             nod)

            CALL DEFGRAINC(igauss,auxm_inc,upd_coord,points,jacFinc,ndim,nod)

            CALL PLASTICITY(deeF,jacF,jacFinc,sigma1C,statev,lnstrainelas,     &
                            sigma,detF,statevar_num,iel,igauss,noncon_flag,    &
                            umat_vm_nonlinear_hardening)
             
            ! Homogenized stress and strain                
            ! Calculation of the Green-Lagrange strain
            ! ------------------------------------------------------
            C_mat=MATMUL(TRANSPOSE(jacF),jacF)
            E_mat=C_mat
            DO i=1,3
              E_mat(i,i)=E_mat(i,i)-one
            END DO
            E_mat=half*E_mat
             
            ! Calculation of the second Piola-Kirchhoff stress
            ! ------------------------------------------------------
            CALL invert_matrix_3x3(jacF,jacFinv)
            jacftrans=TRANSPOSE(jacF)
            CALL invert_matrix_3x3 (jacftrans,jacFtransinv)
            S_mat=detF*MATMUL(jacFinv,MATMUL(sigma,jacFtransinv))
             
            ! ------------------------------------------------------
            sum_stress_pp=sum_stress_pp+S_mat*det*weights(igauss)
            sum_strain_pp=sum_strain_pp+E_mat*det*weights(igauss)                
            ! ------------------------------------------------------       
             
            ! Check for unloading
            IF ((statev(1)>tol).AND.(statev(1)<=(statevar_con(iel,igauss,1)+  &
             tol_val))) THEN
              IF (unload_pp(iel,igauss)==0) THEN
                unload_pp(iel,igauss)=1
                unload_ip_pp=unload_ip_pp+1
              END IF
            END IF
             
            IF ((statev(1)>tol).AND.(iload==max_inc)) THEN
              yield_ip_pp=yield_ip_pp+1
            END IF

            ! Save the variables
            statevar_con(iel,igauss,:)=statev(:)
            lnstrainelas_mat_con(iel,igauss,:)=lnstrainelas(:)
            sigma1C_mat_con(iel,igauss,:)=sigma1C(:)
            stiffness_mat_con(iel,igauss,:,:)=deeF(:,:)
          END DO
        END DO
        
        ! Save the converged displacement
        xnew_pp_previous=xnew_pp

        EXIT
      END IF

    END DO iterations

!------------------------------------------------------------------------------
!------------------------- End Newton-Raphson iterations ----------------------
!------------------------------------------------------------------------------

    IF (numpe==1) THEN
      WRITE(11,'(a,i3,a,f12.4,a,i4,a)') "Time after load step ",iload,": ", &
      ELAP_TIME() - timest(1),"     (",inewton," iterations )"
    END IF

    !IF (iload==1) THEN
    !  ALLOCATE(disp_pp(ndim,nn_pp))
    !END IF

    !IF (iload==num_load_steps) THEN
    !  DEALLOCATE(diag_precon_tmp)
    !END IF

!------------------------------------------------------------------------------
!-----------------------------print out results -------------------------------
!------------------------------------------------------------------------------
    IF (numpe==1) THEN
      IF (iload==1) THEN
        ! Displacement and total strain
        fname = fname_base(1:INDEX(fname_base, " ")-1) // "_dis.res"
        OPEN(24, file=fname, status='replace', action='write')
        fname = fname_base(1:INDEX(fname_base, " ")-1) // "_est.res"
        OPEN(27, file=fname, status='replace', action='write')
        
        !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_str.res"
        !OPEN(25, file=fname, status='replace', action='write')
        fname = fname_base(1:INDEX(fname_base, " ")-1) // "_rea.res"
        OPEN(26, file=fname, status='replace', action='write')
        !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_fint.res"
        !OPEN(28, file=fname, status='replace', action='write')
                
        ! Homogenized stress and strain
        fname = fname_base(1:INDEX(fname_base, " ")-1) // "_hom_stress.res"
        OPEN(29, file=fname, status='replace', action='write')
        fname = fname_base(1:INDEX(fname_base, " ")-1) // "_hom_strain.res"
        OPEN(30, file=fname, status='replace', action='write')
        
        ! Load and displacement
        fname = fname_base(1:INDEX(fname_base, " ")-1) // "_disp_load.res"
        OPEN(31, file=fname, status='replace', action='write')
      END IF
    END IF

    ! Homogenized stress and strain
    sum_stress = zero
    sum_strain = zero

    !CALL MPI_ALLREDUCE(sum_stress_pp,sum_stress,9,MPI_REAL8,MPI_SUM,          &
    !                   MPI_COMM_WORLD,ier) 

    !CALL MPI_ALLREDUCE(sum_strain_pp,sum_strain,9,MPI_REAL8,MPI_SUM,          &
    !                   MPI_COMM_WORLD,ier) 

    IF(numpe==1) WRITE(29,*) sum_stress
    IF(numpe==1) WRITE(30,*) sum_strain
    
    ! Calculate the load
    IF (iload==1) THEN
      ALLOCATE(reacnodes_pp(nodes_pp*nodof))
    END IF
    
    CALL SCATTER_NODES(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,node_start,&
     node_end,storefint_pp,reacnodes_pp,0)
     
    max_disp_inc = max_disp*lambda_total

    sum_force_2_pp=zero
    DO i=1,nodes_pp
      DO j=1,fixed_nodes
        IF (node_start+i-1 == fixed_node(j)) THEN
          sum_force_2_pp=sum_force_2_pp+reacnodes_pp((i-1)*3+2)
        END IF
      END DO
    END DO
    CALL MPI_REDUCE(sum_force_2_pp,sum_force_2,1,MPI_REAL8,MPI_SUM,0,          &
                    MPI_COMM_WORLD,ier)

    IF (numpe==1) THEN
      WRITE(31,*) max_disp_inc,sum_force_2
    END IF
     
    !CALL WRITE_LOAD_DISP(31,nodes_pp,npes,numpe,reacnodes_pp,max_disp_inc,        &
    ! load_node,iload)
     
    IF (numpe==1) THEN
      CALL FLUSH(29)
      CALL FLUSH(30)
      CALL FLUSH(31)
      CALL FLUSH(26)
    END IF
    
!-----print out displacements, stress, principal stress and reactions -------
    !IF (print_output) THEN
    !IF ((iload==1).OR.(lambda_total>=(1-tol_increment)))  THEN
    !IF ((iload==1).OR.(iload==8)) THEN
    !IF (iload==8) THEN
      

      
      
      writetimes = writetimes + 1
      IF(timewrite) THEN
        timest(4) = ELAP_TIME( )
      END IF
      
      ALLOCATE(xnewnodes_pp(nodes_pp*nodof))
      ALLOCATE(shape_integral_pp(nod,nels_pp))
      !ALLOCATE(stress_integral_pp(nod*nst,nels_pp))
      ALLOCATE(strain_integral_pp(nod*nst,nels_pp))
      !ALLOCATE(stressnodes_pp(nodes_pp*nst))
      ALLOCATE(strainnodes_pp(nodes_pp*nst))
      !ALLOCATE(reacnodes_pp(nodes_pp*nodof))
      
      CALL GATHER(xnew_pp(1:),xnewel_pp)
      IF (numfix_pp > 0) THEN
        DO i = 1,numfix_pp
          xnewel_pp(fixdof_pp(i),fixelem_pp(i)) = fixvaltot_pp(i)
        END DO
      END IF

      shape_integral_pp     = zero
      !stress_integral_pp    = zero
      strain_integral_pp    = zero

      DO iel = 1, nels_pp
        DO i = 1, nod
          num(i) = g_num_pp(i,iel) - nn_start + 1
        END DO

        coord = TRANSPOSE(g_coord_pp(:,num))
        DO i = 1,ndim
          auxm(:,i) = xnewel_pp(comp(:,i),iel)
        END DO

        DO igauss = 1,nip
          CALL SHAPE_FUNCTIONS(igauss,points,value_shape)
          dw = det * weights(igauss)
!         DO i = 1,nod
            shape_integral_pp(:,iel) = shape_integral_pp(:,iel) + &
                                       value_shape*dw
!         END DO
        END DO

        !DO igauss = 1,nip
        !  CALL SHAPE_FUNCTIONS(igauss,points,value_shape)
        !  dw = det * weights(igauss)
        !  DO i = 1,nod
        !    idx1 = (i-1)*nst 
        !    DO j = 1,nst
        !      stress_integral_pp(idx1+j,iel) = stress_integral_pp(idx1+j,iel) +&
        !                value_shape(i)*sigma1C_mat_con(j,igauss,iel)*dw
        !    END DO
        !  END DO
        !END DO

        DO igauss = 1,nip
          CALL SHAPE_FUNCTIONS(igauss,points,value_shape)
          dw = det * weights(igauss)
          DO i = 1,nod
            idx1 = (i-1)*nst 
            DO j=1,nst
              strain_integral_pp(idx1+j,iel) = strain_integral_pp(idx1+j,iel) +&
                        value_shape(i)*lnstrainelas_mat_con(iel,igauss,j)*dw
            END DO
          END DO
        END DO
      END DO

      text = "*DISPLACEMENT"
      CALL SCATTER_NODES(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp, &
              node_start,node_end,xnewel_pp,xnewnodes_pp,1)
      CALL WRITE_NODAL_VARIABLE(text,24,iload,nodes_pp,npes,numpe,nodof, &
                                xnewnodes_pp)
	  DEALLOCATE(xnewnodes_pp)

!      text = "*STRESS"
      !CALL NODAL_PROJECTION(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,  &
      ! node_start,node_end,shape_integral_pp,stress_integral_pp,stressnodes_pp)
      !CALL WRITE_NODAL_VARIABLE(text,25,iload,nodes_pp,npes,numpe,nst,   &
      !                          stressnodes_pp)
      !DEALLOCATE(stress_integral_pp,stressnodes_pp)

      text = "*NODAL REACTIONS"
      CALL SCATTER_NODES(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp, &
              node_start,node_end,storefint_pp,reacnodes_pp,0)
      CALL WRITE_NODAL_VARIABLE(text,26,iload,nodes_pp,npes,numpe,nodof, &
                                reacnodes_pp)
      !DEALLOCATE(reacnodes_pp)

      text = "*ELASTIC STRAIN"
      CALL NODAL_PROJECTION(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,  &
       node_start,node_end,shape_integral_pp,strain_integral_pp,strainnodes_pp)
      CALL WRITE_NODAL_VARIABLE(text,27,iload,nodes_pp,npes,numpe,nst,   &
                                strainnodes_pp)
                                
      DEALLOCATE(strain_integral_pp,strainnodes_pp)
      DEALLOCATE(shape_integral_pp)

      IF(timewrite) THEN
        timest(5) = elap_time( )
        timewrite = .FALSE.
      END IF

      ! Outputting unloading
      CALL MPI_ALLREDUCE(yield_ip_pp,yield_tot,1,MPI_INT,MPI_SUM,             &
       MPI_COMM_WORLD,ier) 

      CALL MPI_ALLREDUCE(unload_ip_pp,unload_tot,1,MPI_INT,MPI_SUM,           &
       MPI_COMM_WORLD,ier) 
      
      IF (numpe==1) THEN
        WRITE(*,*) 'The total number of yielded integration points is ',      &
         yield_tot
        WRITE(*,*) 'The total number of unloaded integration points is ',     &
         unload_tot
        WRITE(*,*) 'The percentage of yielded over total is ',                &
         FLOAT(yield_tot)/FLOAT(nels*8)
        WRITE(*,*) 'The percentage of unloaded over total is ',               &
         FLOAT(unload_tot)/FLOAT(nels*8)
        WRITE(*,*) 'The percentage of unloaded over yielded is ',             &
         FLOAT(unload_tot)/FLOAT(yield_tot*8)
      END IF
      
      print_output=.false.

    !END IF  !printing
  END DO !iload
  
  
  
  
  

  
  
  
  
  IF(numpe==1) THEN
    !CLOSE(24)
    !CLOSE(25)
    !CLOSE(26)
    !CLOSE(27)
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------- End Load LOOP --------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  IF (solvers == petsc_solvers) THEN
    DEALLOCATE(p_rows,p_cols,p_values)
    CALL KSPDestroy(p_ksp,p_ierr)
    CALL VecDestroy(p_x,p_ierr)
    CALL VecDestroy(p_b,p_ierr)
    CALL MatDestroy(p_A,p_ierr)
  END IF

  IF (numpe==1) THEN
    WRITE(11,'(a,i5,a)') "This job ran on ",npes," processors"
    WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nels," elements and ",&
                           neq," equations"
    WRITE(11,*) "Time after the mesh   :", timest(2) - timest(1)
    WRITE(11,*) "Time after boundary conditions  :", timest(3) - timest(1)
    WRITE(11,*) "This analysis took  :", elap_time( ) - timest(1)
    WRITE(11,*) "Time to write results (each time) :", timest(5) - timest(4)
    WRITE(11,*) "Time inside the load loop  :", elap_time( ) - timest(3) - &
                                        writetimes*(timest(5)-timest(4))
!   CALL FLUSH(11)
    CLOSE(11)
  END IF

!   Formats
  2000 FORMAT(' Energy  ',i3,1p,i3,1p,e25.15,1p,e25.15) 

  IF (numpe==1) THEN
    WRITE(*,*) 'The simulation is finished'
  END IF

!---------------------------------- shutdown ----------------------------------

  IF (fixed_nodes > 0) THEN
    DEALLOCATE(fixed_node)
  END IF

  IF (solvers == petsc_solvers) THEN
    CALL PetscFinalize(p_ierr)
  END IF
  CALL SHUTDOWN()

CONTAINS
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE umat_vm_nonlinear_hardening(stress,statev,ddsdde,stran,dstran,    &
                                         ntens,statevar_num,iel,igauss,        &
                                         noncon_flag)
    ! This subroutine returns the updated stress, strain and the tangent 
    ! operator for the integrated constitutive law
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ntens
    INTEGER, INTENT(IN) :: statevar_num, iel, igauss
    REAL(iwp), INTENT(OUT) :: stress(:), ddsdde(:,:)
    REAL(iwp), INTENT(INOUT) :: stran(:), statev(:), dstran(:)
    LOGICAL, INTENT(INOUT)   :: noncon_flag

    integer :: i, j, iter

    real(iwp) :: bulk_mod, shear_mod, scalar_term, smises, syield, eqplas, &
     tr, plastic_mul, e, nu, h, syield0, tol, sdev_norm, smises_trial, delta, tangent, &
     res_der, yield_fun, stress_inf
    real(iwp) :: stress_dev(6), proj_tensor(6,6), flow_flow(6,6), unit_tensor(6,6), &
     ixi_tensor(6,6), flow(6), term(6,6), unit_sym(6,6)
     
    real(iwp), PARAMETER :: zero=0._iwp, half=0.5_iwp, one=1._iwp, two=2._iwp, three=3._iwp, six=6._iwp
    integer, PARAMETER :: max_iter=25
   
    ! Assign the tolerance to check converge of the newton-raphson
    tol=0.000001_iwp

    ! Assign material properties
    e=206900._iwp
    nu=0.29_iwp
    h=129.24_iwp
    syield0=450._iwp
    stress_inf=715._iwp
    delta=16.93_iwp

    ! Recover the previous equivalent plastic strain
    eqplas=statev(1)

    ! Setting the fourth and second-order tensors to zero
    proj_tensor=zero
    flow_flow=zero
    unit_tensor=zero
    ixi_tensor=zero
    ddsdde=zero
    stress=zero

    ! Assigning values to the shear and bulk modulus
    bulk_mod=e/(three*(one-two*nu))
    shear_mod=e/(two*(one+nu))

    ! Compute the linear elastic isotropic stiffness matrix
    scalar_term=e/((one+nu)*(one-two*nu))

    ddsdde(1,1)=one-nu
    ddsdde(2,2)=one-nu
    ddsdde(3,3)=one-nu
    ddsdde(4,4)=(one-two*nu)/two
    ddsdde(5,5)=(one-two*nu)/two
    ddsdde(6,6)=(one-two*nu)/two
    ddsdde(1,2)=nu
    ddsdde(1,3)=nu
    ddsdde(2,1)=nu
    ddsdde(2,3)=nu
    ddsdde(3,1)=nu
    ddsdde(3,2)=nu

    ddsdde=scalar_term*ddsdde

    ! Calculate the predictor strain and stress
    DO i=1,ntens
      DO j=1,ntens
        stress(i)=stress(i)+ddsdde(i,j)*stran(j)
      END DO
    END DO

    stress_dev=stress
    tr=stress(1)+stress(2)+stress(3)
    stress_dev(1)=stress_dev(1)-(one/three)*tr
    stress_dev(2)=stress_dev(2)-(one/three)*tr
    stress_dev(3)=stress_dev(3)-(one/three)*tr

    ! Calculate the equivalent Von Mises stress
    smises_trial=((stress(1)-stress(2))**2)+((stress(2)-stress(3))**2)+ &
     ((stress(3)-stress(1))**2)
    DO i=4,6
      smises_trial=smises_trial+six*(stress(i)**2)
    END DO
    smises_trial=SQRT(smises_trial/two)

    ! Calculate the yield stress
    syield=(stress_inf-syield0)*(one-EXP(-delta*eqplas))+h*eqplas+syield0   

    ! Determine if there is yielding
    IF (smises_trial>(syield+tol)) THEN
 
      ! Calculate the plastic multiplicator through a newton-raphson scheme
      ! Initialize the plastic multiplicator and the residual function
      plastic_mul=zero
      yield_fun=smises_trial-syield

      DO iter=1,max_iter
        ! If convergence is not achieved
        IF (i==max_iter) THEN
          WRITE(*,*) 'Maximum local Newton-Raphson iterations have been reached'
          noncon_flag=.TRUE.
          EXIT
        END IF
        
        ! Calculate the hardening slope
        tangent=(stress_inf-syield0)*delta*EXP(-delta*(eqplas+plastic_mul))+h
      
        ! Calculate the residual derivative
        res_der=-three*shear_mod-tangent

        ! New guess for the plastic multiplicator
        plastic_mul=plastic_mul-(yield_fun/res_der)

        ! Update stresses
        DO i=1,3
          stress(i)=stress_dev(i)*(one-((plastic_mul*three*shear_mod)/smises_trial))+ &
           (one/three)*tr
        END DO

        DO i=4,6
          stress(i)=stress_dev(i)*(one-((plastic_mul*three*shear_mod)/smises_trial))
        END DO
        
        ! Calculate the equivalent Von Mises stress
        smises=((stress(1)-stress(2))**2)+((stress(2)-stress(3))**2)+ &
         ((stress(3)-stress(1))**2)
        DO i=4,6
          smises=smises+six*(stress(i)**2)
        END DO
        smises=SQRT(smises/two)

        ! Check for convergence
        syield=(stress_inf-syield0)*(one-EXP(-delta*(eqplas+plastic_mul)))+ &
         h*(eqplas+plastic_mul)+syield0
        yield_fun=smises-syield

        IF (yield_fun<tol) THEN
          EXIT
        END IF
      END DO

      ! Create the unit flow vector
      sdev_norm=((stress_dev(1)**2)+(stress_dev(2)**2)+(stress_dev(3)**2)+ &
       (two*stress_dev(4)**2)+(two*stress_dev(5)**2)+(two*stress_dev(6)**2))
      sdev_norm=SQRT(sdev_norm)
      flow=stress_dev/sdev_norm

      ! Update the equivalent plastic strain
      eqplas=eqplas+plastic_mul

      ! Update the elastic strains
      DO i=1,3
        stran(i)=stran(i)-plastic_mul*(SQRT(three/two))*flow(i)
      END DO

      DO i=4,6
        stran(i)=stran(i)-two*plastic_mul*(SQRT(three/two))*flow(i)
      END DO

      ! Assemble the tangent operator
      ! Create the IxI tensor
      DO i=1,3
        DO j=1,3
          ixi_tensor(i,j)=one
        END DO
      END DO

      ! Create the unit tensor
      DO i=1,3
        unit_tensor(i,i)=one
      END DO
      unit_tensor(4,4)=half
      unit_tensor(5,5)=half
      unit_tensor(6,6)=half

      DO i=1,6
        DO j=1,6
          ! Create the projection tensor
          proj_tensor(i,j)=unit_tensor(i,j)-(one/three)*ixi_tensor(i,j)

          ! Create the tensor product of the two flow vectors
          flow_flow(i,j)=flow(i)*flow(j)
          
          
          
          ! MAYBE FLOW_FLOW NEEDS A MODIFICATION HERE
          
          

          ! Assemble all the parts into the final tangent operator
          ddsdde(i,j)=ddsdde(i,j)-((plastic_mul*six*(shear_mod**2))/smises_trial)*proj_tensor(i,j)+    &
          six*(shear_mod**2)*((plastic_mul/smises_trial)-(one/(three*shear_mod+tangent)))*flow_flow(i,j)
        END DO
      END DO
    END IF

    ! Update the state variables
    statev(1)=eqplas
    
    RETURN
  END SUBROUTINE umat_vm_nonlinear_hardening
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END PROGRAM xx15_vm_nonlinear_hardening
