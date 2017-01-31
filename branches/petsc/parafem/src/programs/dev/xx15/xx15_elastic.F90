PROGRAM xx15_elastic
!------------------------------------------------------------------------------
!   program xx15:  finite strain elasto-plastic analysis with Newton-Raphson
!------------------------------------------------------------------------------

  USE choose_solvers
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

  INTEGER :: nels, nn, nr, nip, nodof=3, nod, nst=6, loaded_nodes, nn_pp,     &
   nf_start, fmt=1, i, j, k, l, ndim=3, iters, limit, iel, nn_start,          &
   num_load_steps, iload, igauss, dimH, inewton, jump, npes_pp, partitioner  ,&
   argc, iargc, limit_2, fixed_nodes, numfix_pp, fixdim, writetimes=0,        &
   nodes_pp, node_start, node_end, idx1, idx2

  REAL(iwp) :: det, tol, maxdiff, tol2, detF, detFinc, energy, energy1, rn0,  &
   timest(40), detF_mean, initial_guess, lambda, lambda_total, lambda_prev,   &
   energy_prev, energy_prev_prev, min_inc, next_output, temp_inc, dw, e, v,   &
   jacFtransinv(3,3), jacFtrans(3,3), sum_strain(3,3)

  REAL(iwp), PARAMETER :: tol_increment=0.000001_iwp, tol_val=0.00000001_iwp
  REAL(iwp), PARAMETER :: zero=0.0_iwp, one=1._iwp, half=0.5_iwp
  INTEGER, PARAMETER :: statevar_num=10
  
  INTEGER :: iter, prev_iter, max_inc, number_output, counter_output

  CHARACTER(len=15) :: element
  CHARACTER(len=50) :: text, fname_base, fname
 
  LOGICAL :: converged, timewrite=.TRUE., flag=.FALSE., print_output=.FALSE., &
   tol_inc=.FALSE., lambda_inc=.TRUE., noncon_flag=.FALSE.

  CHARACTER(len=choose_solvers_string_length) :: solvers
  LOGICAL                  :: error
  CHARACTER(:),ALLOCATABLE :: message
  CHARACTER,PARAMETER      :: tab = ACHAR(9)
  REAL                     :: memory_use,peak_memory_use
  
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
   reacnodes_pp(:), stiffness_mat_con(:,:,:,:), km(:,:), fixkm_pp(:,:,:)

  INTEGER, ALLOCATABLE  :: num(:), g_num(:,:), g_num_pp(:,:), g_g_pp(:,:),    &
   load_node(:), rest(:,:), nf_pp(:,:), no_pp(:), comp(:,:), fixed_node(:),   &
   fixed_dof(:), fixelem_pp(:), fixdof_pp(:), unload_pp(:,:)
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! 1. Input and initialisation
  !----------------------------------------------------------------------------

  timest = zero

  timest(1) = elap_time()
  timest(2) = elap_time()

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
                     nip,limit,tol,e,v,nod,num_load_steps,jump,tol2,partitioner)

  solvers = get_solvers()
  IF (.NOT. solvers_valid(solvers)) THEN
    CALL SHUTDOWN
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

!------------------------------------------------------------------------------
! 2. Get integration Gauss points and weights in the element
!------------------------------------------------------------------------------
  
  ALLOCATE(points(ndim,nip), weights(nip))
  CALL GET_GAUSS_POINTS(element,points,weights)

!------------------------------------------------------------------------------
! 3. Import and distribute mesh
!------------------------------------------------------------------------------

  CALL CALC_NELS_PP(fname_base,nels,npes,numpe,partitioner,nels_pp)
  
  ALLOCATE(g_num_pp(nod, nels_pp)) 
  
  fname = fname_base(1:INDEX(fname_base," ")-1) // ".d"
  CALL READ_G_NUM_PP(fname_base,iel_start,nn,npes,numpe,g_num_pp)
  ! read_elements_2() does not work for METIS partitions
  CALL READ_ELEMENTS_2(fname,npes,nn,numpe,g_num_pp)
  
  CALL CALC_NN_PP(g_num_pp,nn_pp,nn_start)
  
  ALLOCATE(g_coord_pp(ndim, nn_pp))
  CALL READ_NODES(fname,nn,nn_start,numpe,g_coord_pp)
  
  ALLOCATE(rest(nr,nodof+1))
  fname = fname_base(1:INDEX(fname_base," ")-1) // ".bnd"
  CALL READ_RESTRAINTS(fname,numpe,rest)
  
  timest(30) = timest(30) + elap_time()-timest(2) ! 30 = read
  timest(2) = elap_time()

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
  IF (solvers == parafem_solvers) THEN
    ALLOCATE(storekm_pp(ntot,ntot,nels_pp))
  END IF
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
  ALLOCATE(km(ntot,ntot))

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of 
!    equations to solve. 
!------------------------------------------------------------------------------

  CALL REST_NF(nn_start,rest,nf_pp)

  DO iel = 1,nels_pp
    CALL NUM_TO_G2(g_num_pp(:,iel), nf_pp, g_g_pp(:,iel), nn_start)
  END DO

  CALL CALC_NEQ(nn,rest,neq)
  
!------------------------------------------------------------------------------
! 6. Create interprocessor communications tables
!------------------------------------------------------------------------------  
  
  CALL CALC_NPES_PP(npes,npes_pp)

  CALL CALC_NEQ_PP

  CALL MAKE_GGL(npes_pp,npes,g_g_pp)

  timest(31) = timest(31) + elap_time()-timest(2) ! 31 = setup
  timest(2) = elap_time()
  
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
    
    DEALLOCATE(fixed_node, fixed_dof, fixed_value)

    IF (solvers == petsc_solvers) THEN
      ! storekm_pp is only used by the ParaFEM solvers but would also be needed
      ! to convert the fixed displacements to a force.  Use a (hopefully small)
      ! copy of the elements needed for the fixed displacements.  It would
      ! better to use PETSc matrix-vector multiplies but that would require all
      ! the handling of fixed dofs in ParaFEM to be changed (right down to
      ! changes in gather and scatter I think) to pass them through as zeroes in
      ! the matrix.
      ALLOCATE(fixkm_pp(ntot,ntot,numfix_pp))
    END IF

  END IF
  
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
  
  timest(30) = timest(30) + elap_time()-timest(2) ! 30 = read
  timest(2) = elap_time()

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

  !-----------------------------------------------------------------------------
  ! 9a. Start up PETSc after find_pe_procs (so that MPI has been started) and
  !     after find_g3 so that neq_pp and g_g_pp are set up.
  !-----------------------------------------------------------------------------
  IF (solvers == petsc_solvers) THEN
    CALL p_initialize(fname_base,error)
    IF (error) THEN
      CALL shutdown
    END IF
    ! Set up PETSc.
    CALL p_create(ntot,g_g_pp,error)
    IF (error) THEN
      CALL p_finalize
      CALL shutdown
    END IF
  END IF
  memory_use = p_memory_use()
  peak_memory_use = p_memory_peak()
  IF (numpe == 1) WRITE(*,'(A,2F7.2,A)')                                       &
    "current and peak memory use after setup:        ",                        &
    memory_use,peak_memory_use," GB "

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
  min_inc=0.001_iwp

  ! Establish the number of times the output is printed (not implemented yet)
  number_output=1
  counter_output=1
  
  ! Set the initial guess (normalized)
  initial_guess  = 0.01_iwp
  lambda         = initial_guess
  lambda_total   = zero
  lambda_prev    = zero
  next_output=one/number_output
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!----------------------- Start Load LOOP --------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  timest(31) = timest(31) + elap_time()-timest(2) ! 31 = setup
  timest(2) = elap_time()

  DO iload = 1,max_inc
  
    converged = .FALSE.

    timest(2) = elap_time()

    ! Increase the stage control by 50% if the number of iterations of the two 
    ! previously converged full increments (not affected by the output  
    ! requests) are below 6
    IF ((iter<6).AND.(prev_iter<6).AND.(iload>2)) THEN
      IF (numpe==1) THEN
        WRITE(*,*) 'The load increment is increased by 50%'
      END IF
      lambda=1.5*lambda
    END IF

    timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
    timest(2) = elap_time()

    100 CONTINUE

    timest(2) = elap_time()

    ! If lambda is larger than unity, the simulation is finished
    IF (lambda_total>=(1-tol_increment)) THEN

      timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
      timest(2) = elap_time()
      
      EXIT
    END IF

    ! Update the load increments
    lambda_prev=lambda_total
    lambda_total=lambda_total+lambda
    ! Don't overshoot.
    IF (lambda_total>1._iwp) THEN
      lambda_total = 1._iwp
      lambda = lambda_total - lambda_prev
    END IF

    ! Exit if the time increment is less that the specified minimum
    IF (lambda<min_inc) THEN
      IF(numpe==1) THEN
        WRITE(*,*) 'The load increment is too small'
      END IF

      timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
      timest(2) = elap_time()
      
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

!------------------------------------------------------------------------------
!----------------------- Start Newton-Raphson iterations ----------------------
!------------------------------------------------------------------------------

    deltax_pp_temp = zero
    inewton = 0
    
    timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
    timest(2) = elap_time()

    iterations: DO

      timest(2) = elap_time()

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
      
      timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
      timest(2) = elap_time()
      
      IF (solvers == parafem_solvers) THEN
        storekm_pp = zero
      ELSE IF (solvers == petsc_solvers) THEN
        CALL p_zero_matrix
      END IF
      
      DO iel = 1,nels_pp

        km = zero

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
        
        DO igauss = 1,nip

          ! Initialise the state variables to the same converged value of     
          ! the last time increment
          statev(:)=statevar_con(iel,igauss,:)
          lnstrainelas(:)=lnstrainelas_mat_con(iel,igauss,:)

          ! Calculates the total and incremental deformation gradients
          CALL DEFGRA(igauss,auxm,coord,points,det,detF,beeF,geeF,jacF,ndim,  &
           nod)
          CALL DEFGRAINC(igauss,auxm_inc,upd_coord,points,jacFinc,ndim,nod)
          
          CALL PLASTICITY(deeF,jacF,jacFinc,sigma1C,statev,lnstrainelas,       &
                          sigma,detF,statevar_num,iel,igauss,noncon_flag,      &
                          umat_elastic)
           
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
          km = km + (MATMUL(MATMUL(geeFT,deeF),geeF)*dw)

        END DO

        IF (solvers == parafem_solvers) THEN
          storekm_pp(:,:,iel) = km
        ELSE IF (solvers == petsc_solvers) THEN
          CALL p_add_element(g_g_pp(:,iel),km)
          FORALL (i = 1:numfix_pp, fixelem_pp(i) == iel) fixkm_pp(:,:,i) = km
        END IF
      END DO
      memory_use = p_memory_use()
      peak_memory_use = p_memory_peak()
      IF (numpe == 1) WRITE(*,'(A,2F7.2,A)')                                  &
        "current and peak memory use after add elements: ",                   &
         memory_use,peak_memory_use," GB "

      IF (solvers == petsc_solvers) THEN
        CALL p_assemble
      END IF
      memory_use = p_memory_use()
      peak_memory_use = p_memory_peak()
      IF (numpe == 1) WRITE(*,'(A,2F7.2,A)')                                  &
        "current and peak memory use after assemble:     ",                   &
        memory_use,peak_memory_use," GB "
      
      timest(33) = timest(33) + elap_time()-timest(2) ! 33 = matrix assemble
      timest(2) = elap_time()
      
!------------------------------------------------------------------------------
! 12. Build and invert the preconditioner (ParaFEM only)
!------------------------------------------------------------------------------

      IF (solvers == parafem_solvers) THEN
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

      timest(34) = timest(34) + elap_time()-timest(2) ! 34 = solve
      timest(2) = elap_time()

!------------------------------------------------------------------------------
! 13. Initialize PCG
!------------------------------------------------------------------------------
      
      ! During the first Newton-Raphson iteration, the incremental 
      ! displacements are applied through a linear mapping of these 
      ! displacements to the internal force vector
      IF (inewton==1 .AND. numfix_pp>0) THEN
        DO i = 1,numfix_pp
          DO j = 1,ntot
            IF (g_g_pp(j,fixelem_pp(i))>0) THEN
              IF (solvers == parafem_solvers) THEN
                storefint_pp(j,fixelem_pp(i)) = storefint_pp(j,fixelem_pp(i))  &
                  + fixvalpiece_pp(i)*storekm_pp(j,fixdof_pp(i),fixelem_pp(i))
              ELSE IF (solvers == petsc_solvers) THEN
                storefint_pp(j,fixelem_pp(i)) = storefint_pp(j,fixelem_pp(i))  &
                  + fixvalpiece_pp(i)*fixkm_pp(j,fixdof_pp(i),i)
              END IF
            END IF
          END DO
        END DO
      END IF
      
      fint_pp = zero
      CALL SCATTER(fint_pp(1:),storefint_pp)

      r_pp(1:) = fext_pp(1:) - fint_pp(1:) 
      r_pp(0) = zero

      ! Compute maxdiff of residual 
      maxdiff =  MAXABSVAL_P(r_pp(1:))

      ! Normalise residual vector and stiffness matrix for pcg
      IF (maxdiff == zero) THEN
        IF(numpe==1) THEN
          WRITE(*,*) "maxdiff = zero and now exiting loop"
        END IF

        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()
      
        EXIT
      END IF

!------------------------------------------------------------------------------
!----------------- Solve using preconditioned Krylov solver -------------------
!------------------------------------------------------------------------------
      
      deltax_pp = zero
      res_pp    = r_pp

      timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
      timest(2) = elap_time()

      IF (solvers == parafem_solvers) THEN
        CALL PCG_VER1(inewton,limit,tol,storekm_pp,r_pp(1:),                   &
                      diag_precon_pp(1:),rn0,deltax_pp(1:),iters)
      ELSE IF (solvers == petsc_solvers) THEN
        CALL p_use_solver(1,error)
        IF (error) THEN
          CALL p_shutdown
          CALL shutdown
        END IF
        CALL p_solve(r_pp(1:),deltax_pp(1:))
        CALL p_print_info(1,11)
      END IF

      timest(34) = timest(34) + elap_time()-timest(2) ! 34 = solve
      timest(2) = elap_time()

      IF (numpe==1) THEN
        WRITE(91,*)iload,inewton,iters
        CALL FLUSH(91)
      END IF

      timest(35) = timest(35) + elap_time()-timest(2) ! 35 = write
      timest(2) = elap_time()

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

      timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
      timest(2) = elap_time()

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
        IF (print_output) THEN
          lambda_total=lambda_total-temp_inc
        ELSE
          lambda_total=lambda_total-lambda
        END IF
        lambda=half*lambda
        xnew_pp=xnewprev_pp
        flag=.FALSE.
        
        ! This is not yet implemented
        !IF (print_output) THEN
        !  counter_output=counter_output-1
        !  next_output=FLOAT(counter_output)/number_output
        !  print_output=.FALSE.
        !END IF
        
        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()

        GOTO 100
      END IF
 
      ! After convergence, a last "iteration" is needed to calculate the fully
      ! converged values (some FE algorithms omit this step, but we perform
      ! it because we are interested in a fully precise solution)
      IF (converged) THEN 

        ! If the current increment is not an output request increment, save the  
        ! current and previous increment number of iterations to allow for   
        ! for lambda increment to be increased (not yet implemented)
        !IF (.NOT.(lambda_inc)) THEN
        !  IF (iload>1) THEN
        !    prev_iter=iter
        !  END IF
        !  iter=inewton
        !  lambda_inc=.TRUE.
        !END IF

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

          DO igauss = 1,nip

            statev(:)=statevar_con(iel,igauss,:)
            lnstrainelas(:)=lnstrainelas_mat_con(iel,igauss,:)

            CALL DEFGRA(igauss,auxm,coord,points,det,detF,beeF,geeF,jacF,ndim,&
             nod)

            CALL DEFGRAINC(igauss,auxm_inc,upd_coord,points,jacFinc,ndim,nod)

            CALL PLASTICITY(deeF,jacF,jacFinc,sigma1C,statev,lnstrainelas,     &
                            sigma,detF,statevar_num,iel,igauss,noncon_flag,    &
                            umat_elastic)

            ! Save the variables
            statevar_con(iel,igauss,:)=statev(:)
            lnstrainelas_mat_con(iel,igauss,:)=lnstrainelas(:)
            sigma1C_mat_con(iel,igauss,:)=sigma1C(:)
            stiffness_mat_con(iel,igauss,:,:)=deeF(:,:)
          END DO
        END DO
        
        ! Save the converged displacement
        xnew_pp_previous=xnew_pp

        timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
        timest(2) = elap_time()

        EXIT
      END IF

      timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in load loop
      timest(2) = elap_time()

    END DO iterations

    timest(2) = elap_time()

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
        !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_rea.res"
        !OPEN(26, file=fname, status='replace', action='write')
        !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_fint.res"
        !OPEN(28, file=fname, status='replace', action='write')
                
        ! Homogenized stress and strain
        !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_hom_stress.res"
        !OPEN(29, file=fname, status='replace', action='write')
        !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_hom_strain.res"
        !OPEN(30, file=fname, status='replace', action='write')
        
        ! Load and displacement
        !fname = fname_base(1:INDEX(fname_base, " ")-1) // "_disp_load.res"
        !OPEN(31, file=fname, status='replace', action='write')
      END IF
    END IF

    IF (numpe==1) THEN
      !CALL FLUSH(29)
      !CALL FLUSH(30)
      !CALL FLUSH(31)
    END IF
    
!-----print out displacements, stress, principal stress and reactions -------
    !IF (print_output) THEN
!!$    IF (iload==max_inc) THEN
      
      writetimes = writetimes + 1
      
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

        !coord = TRANSPOSE(g_coord_pp(:,num))
        !DO i = 1,ndim
        !  auxm(:,i) = xnewel_pp(comp(:,i),iel)
        !END DO

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

!      text = "*NODAL REACTIONS"
      !CALL SCATTER_NODES(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp, &
      !        node_start,node_end,storefint_pp,reacnodes_pp,0)
      !CALL WRITE_NODAL_VARIABLE(text,26,iload,nodes_pp,npes,numpe,nodof, &
      !                          reacnodes_pp)
      !DEALLOCATE(reacnodes_pp)

      text = "*ELASTIC STRAIN"
      CALL NODAL_PROJECTION(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,  &
       node_start,node_end,shape_integral_pp,strain_integral_pp,strainnodes_pp)
      CALL WRITE_NODAL_VARIABLE(text,27,iload,nodes_pp,npes,numpe,nst,   &
                                strainnodes_pp)
                                
      DEALLOCATE(strain_integral_pp,strainnodes_pp)
      DEALLOCATE(shape_integral_pp)

      ! Outputting unloading
      !CALL MPI_ALLREDUCE(yield_ip_pp,yield_tot,1,MPI_INT,MPI_SUM,             &
      ! MPI_COMM_WORLD,ier) 

      !CALL MPI_ALLREDUCE(unload_ip_pp,unload_tot,1,MPI_INT,MPI_SUM,           &
      ! MPI_COMM_WORLD,ier) 
      
      !IF (numpe==1) THEN
      !  WRITE(*,*) 'The total number of yielded integration points is ',      &
      !   yield_tot
      !  WRITE(*,*) 'The total number of unloaded integration points is ',     &
      !   unload_tot
      !  WRITE(*,*) 'The percentage of yielded over total is ',                &
      !   FLOAT(yield_tot)/FLOAT(nels*8)
      !  WRITE(*,*) 'The percentage of unloaded over total is ',               &
      !   FLOAT(unload_tot)/FLOAT(nels*8)
      !  WRITE(*,*) 'The percentage of unloaded over yielded is ',             &
      !   FLOAT(unload_tot)/FLOAT(yield_tot*8)
      !END IF
      
      print_output=.false.

!!$    END IF  !printing
    
      timest(35) = timest(35) + elap_time()-timest(2) ! 35 = write
      timest(2) = elap_time()
      
    IF (iload==max_inc) THEN
      EXIT
    END IF
  END DO !iload
  
  timest(2) = elap_time()

  IF(numpe==1) THEN
    CLOSE(24)
    !CLOSE(25)
    !CLOSE(26)
    CLOSE(27)
    !CLOSE(29)
    !CLOSE(30)
    !CLOSE(31)
  END IF

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------- End Load LOOP --------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  peak_memory_use = p_memory_peak()

  IF (numpe==1) THEN
    WRITE(11,'(a,i5,a)') "This job ran on ",npes," processors"
    WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nels," elements and ",&
                           neq," equations"
    WRITE(11,'(A,F10.4)') "Time to read input:       ", timest(30)
    WRITE(11,'(A,F10.4)') "Time for setup:           ", timest(31)
    WRITE(11,'(A,F10.4)') "Time for matrix assemble: ", timest(33)
    WRITE(11,'(A,F10.4)') "Time for linear solve:    ", timest(34)
    WRITE(11,'(A,F10.4)') "Other time in load loop:  ", timest(32)
    WRITE(11,'(A,F10.4)') "Time to write results:    ", timest(35)
    WRITE(11,'(A)')       "                          ----------"
    WRITE(11,'(A,F10.4)') "Total:                    ", SUM(timest(30:35))
    WRITE(11,'(A)')       "                          ----------"
    WRITE(11,'(A,F10.4)') "This analysis took:       ", elap_time()-timest(1)
!   CALL FLUSH(11)
    CLOSE(11)
  END IF

!   Formats
  2000 FORMAT(' Energy  ',i3,1p,i3,1p,e25.15,1p,e25.15) 

  IF (numpe==1) THEN
    WRITE(*,*) 'The simulation is finished'
  END IF

!---------------------------------- shutdown ----------------------------------

  IF (solvers == petsc_solvers) THEN
    CALL p_shutdown
  END IF

  CALL SHUTDOWN()

CONTAINS
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE umat_elastic(stress,statev,ddsdde,stran,dstran,ntens,             &
                          statevar_num,iel,igauss,noncon_flag)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ntens
    INTEGER, INTENT(IN) :: statevar_num, iel, igauss
    REAL(iwp), INTENT(OUT) :: stress(:), ddsdde(:,:)
    REAL(iwp), INTENT(INOUT) :: stran(:), statev(:), dstran(:)
    LOGICAL, INTENT(INOUT)   :: noncon_flag

    REAL(iwp) :: scalar_term, e, nu
    REAL(iwp), PARAMETER :: zero=0._iwp, one=1._iwp, two=2._iwp
    INTEGER :: i, j
	
    ! Assign material properties
    e=12700._iwp
    nu=0.3_iwp
    
    ! Setting the fourth and second-order tensors to zero
    ddsdde=zero
    stress=zero

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
    
    RETURN
  END SUBROUTINE umat_elastic

END PROGRAM xx15_elastic
