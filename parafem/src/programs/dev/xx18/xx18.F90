PROGRAM xx18       
  
  !/****f* dev/xx18
  !*  NAME
  !*    PROGRAM: xx18
  !*  SYNOPSIS
  !*    Usage:   main program
  !*  FUNCTION
  !*    Three dimensional analysis of an elastic solid using 20-node brick
  !*    elements, choice between built in preconditioned conjugate gradient
  !*    solver and PETSc library.  ARCHER eCSE06 project.
  !*    
  !*    Parallel version.  Loaded_nodes only.  See program p121.
  !*    
  !*  AUTHORS
  !*    Lee Margetts, Mark Filipiak
  !*  CREATION DATE
  !*    02.01.2007
  !*  MODIFICATION HISTORY
  !*    Version 2, 19.02.2016, Mark Filipiak
  !*    Version 3, 23.06.2016, Mark Filipiak
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010, University of Edinburgh 2016
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  USE choose_solvers
  USE parafem_petsc
  USE PRECISION; USE global_variables; USE mp_interface; USE input
  USE output; USE loading; USE timing; USE maths; USE gather_scatter
  USE steering; USE new_library; USE large_strain 

  IMPLICIT NONE

  ! neq,ntot are now global variables - must not be declared
  
  INTEGER,PARAMETER   :: nodof=3,ndim=3,nst=6
  INTEGER             :: loaded_nodes,iel,i,j,k,iters,limit,nn,nr,nip,nod
  INTEGER             :: nels,ndof
  INTEGER             :: npes_pp,node_end,node_start,nodes_pp,meshgen
  INTEGER             :: partitioner,nlen
  INTEGER             :: inewton
  REAL(iwp),PARAMETER :: zero=0.0_iwp
  REAL(iwp)           :: rn0
  REAL(iwp)           :: e,v,det,tol,up,alpha,beta,q
  LOGICAL             :: converged=.FALSE.
  CHARACTER(LEN=50)   :: argv
  CHARACTER(LEN=15)   :: element
  CHARACTER(LEN=6)    :: ch 

  CHARACTER(len=choose_solvers_string_length) :: solvers
  LOGICAL                  :: error
  CHARACTER(:),ALLOCATABLE :: message
  CHARACTER,PARAMETER      :: tab = ACHAR(9)
  REAL                     :: peak_memory_use
  
  !-----------------------------------------------------------------------------
  ! 1. Dynamic arrays
  !-----------------------------------------------------------------------------
  
  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:),val(:,:),           &
    disp_pp(:),g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:),bee(:,:),        &
    storkm_pp(:,:,:),eps(:),sigma(:),diag_precon_pp(:),p_pp(:),                &
    r_pp(:),x_pp(:),xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),     &
    timest(:),diag_precon_tmp(:,:),eld_pp(:,:),temp(:),km(:,:)
  INTEGER,ALLOCATABLE :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)
  
  !-----------------------------------------------------------------------------
  ! 2. Input and initialisation
  !-----------------------------------------------------------------------------
  
  ALLOCATE(timest(20))
  timest=zero

  timest(1)=elap_time()

  CALL find_pe_procs(numpe,npes)
  CALL getname(argv,nlen) 
  CALL read_p121(argv,numpe,e,element,limit,loaded_nodes,meshgen,nels,         &
                 nip,nn,nod,nr,partitioner,tol,v)

  solvers = get_solvers(numpe)
  IF (.NOT. solvers_valid(solvers)) THEN
    CALL SHUTDOWN
  END IF

  CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
  ndof=nod*nodof; ntot=ndof
  ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),                &
    rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
  CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
  CALL read_rest(argv,numpe,rest)

  timest(2)=elap_time()

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),der(ndim,nod),         &
    deriv(ndim,nod),bee(nst,ntot),weights(nip),eps(nst),sigma(nst),            &
    km(ntot,ntot),pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),                &
    g_g_pp(ntot,nels_pp))

  IF (solvers == parafem_solvers) THEN
    ALLOCATE(storkm_pp(ntot,ntot,nels_pp))
  END IF

  !-----------------------------------------------------------------------------
  ! 3. Find the steering array and equations per process
  !-----------------------------------------------------------------------------

  CALL rearrange(rest); g_g_pp=0; neq=0
  elements_0: DO iel=1,nels_pp
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
  END DO elements_0
  neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
  ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),             &
    u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp)); diag_precon_pp=zero
  p_pp=zero;  r_pp=zero;  x_pp=zero; xnew_pp=zero; u_pp=zero; d_pp=zero
  
  !-----------------------------------------------------------------------------
  ! 4. Start up PETSc after find_pe_procs (so that MPI has been started)
  !-----------------------------------------------------------------------------
  IF (solvers == petsc_solvers) THEN
    CALL p_initialize(argv,numpe)
    ! Set the approximate number of zeroes per row for the matrix size
    ! pre-allocation.
    CALL p_row_nnz(ndim,nodof,nod,error,message)
    IF (error) THEN
      IF (numpe == 1) THEN
        WRITE(*,'(A)') message
      END IF
      CALL p_finalize
      CALL shutdown
    END IF
    ! Set up PETSc.
    CALL p_setup(neq_pp,ntot,numpe)
  END IF

  IF(numpe==1)THEN
    OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
  END IF

  timest(3) = elap_time()

  !-----------------------------------------------------------------------------
  ! 5. Element stiffness integration and global stiffness matrix creation
  !-----------------------------------------------------------------------------
  
  dee=zero; CALL deemat(dee,e,v); CALL sample(element,points,weights)

  IF (solvers == parafem_solvers) THEN
    storkm_pp=zero
  ELSE IF (solvers == petsc_solvers) THEN
    CALL p_zero_matrix
  END IF

  elements_1: DO iel=1,nels_pp
    km = zero
    gauss_pts_1: DO i=1,nip
      CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
      det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
      CALL beemat(bee,deriv)
      km = km + MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
    END DO gauss_pts_1

    IF (solvers == parafem_solvers) THEN
      storkm_pp(:,:,iel) = km
    ELSE IF (solvers == petsc_solvers) THEN
      CALL p_add_element(g_g_pp(:,iel),km)
    END IF
  END DO elements_1

  IF (solvers == petsc_solvers) THEN
    CALL p_assemble(numpe)
  END IF

  timest(4) = elap_time()
  
  !-----------------------------------------------------------------------------
  ! 6. Build and invert the preconditioner (ParaFEM only)
  !-----------------------------------------------------------------------------

  IF (solvers == parafem_solvers) THEN
    ALLOCATE(diag_precon_tmp(ntot,nels_pp)); diag_precon_tmp=zero
    elements_2: DO iel=1,nels_pp; DO i=1,ndof
      diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel)+storkm_pp(i,i,iel)
    END DO;  END DO elements_2
    CALL scatter(diag_precon_pp,diag_precon_tmp); DEALLOCATE(diag_precon_tmp)
    diag_precon_pp = 1._iwp/diag_precon_pp
  END IF

  timest(5) = elap_time()

  !-----------------------------------------------------------------------------
  ! 7. Get starting r
  !-----------------------------------------------------------------------------
  
  IF(loaded_nodes>0) THEN
    ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes)); node=0; val=zero
    CALL read_loads(argv,numpe,node,val)
    CALL load(g_g_pp,g_num_pp,node,val,r_pp); q=SUM_P(r_pp)
    IF(numpe==1) WRITE(11,'(A,E12.4)') "The total load is: ",q
    DEALLOCATE(node,val)
  END IF
  
  DEALLOCATE(g_g_pp)

  !-----------------------------------------------------------------------------
  ! 8. Preconditioned Krylov solver
  !-----------------------------------------------------------------------------
  
  x_pp      = zero
  
  timest(6) = elap_time()

  IF (solvers == parafem_solvers) THEN
    iters     = 0
    inewton   = 1    ! must be one in this program
    rn0       = zero
    CALL PCG_VER1(inewton,limit,tol,storkm_pp,r_pp,diag_precon_pp,rn0,x_pp,iters)
    xnew_pp = x_pp
    IF(numpe==1)THEN
      WRITE(11,'(A,I6)') "The number of iterations to convergence was ",iters
    END IF
  ELSE IF (solvers == petsc_solvers) THEN
    CALL p_use_solver(1,numpe,error)
    IF (error) THEN
      CALL p_finalize
      CALL shutdown
    END IF
    CALL p_solve(r_pp,xnew_pp)
    CALL p_print_info(numpe,11)
  END IF

  timest(7) = elap_time()
  
  IF(numpe==1)THEN
    WRITE(11,'(A,I7,A)') "This job ran on ",npes," processes"
    WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes", nr,                     &
                             " restrained and ",neq," equations"
    WRITE(11,'(A,E12.4)') "The central nodal displacement is :",xnew_pp(1)
  END IF

  DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,km,pmul_pp)
  IF (solvers == parafem_solvers) THEN
    DEALLOCATE(storkm_pp)
  END IF

  !-----------------------------------------------------------------------------
  ! 9. Recover stresses at centroidal gauss point
  !-----------------------------------------------------------------------------
  
  ALLOCATE(eld_pp(ntot,nels_pp)); eld_pp=zero; points=zero; nip=1; iel=1
  CALL gather(xnew_pp,eld_pp); DEALLOCATE(xnew_pp)
  IF(numpe==1)WRITE(11,'(A)') "The Centroid point stresses for element 1 are"
  gauss_pts_2: DO i=1,nip
    CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
    CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
    eps=MATMUL(bee,eld_pp(:,iel)); sigma=MATMUL(dee,eps)
    IF(numpe==1.AND.i==1) THEN
      WRITE(11,'(A,I5)') "Point ",i ; WRITE(11,'(6E12.4)') sigma
    END IF
  END DO gauss_pts_2; DEALLOCATE(g_coord_pp)
  
  !-----------------------------------------------------------------------------
  ! 10. Write out displacements
  !-----------------------------------------------------------------------------
  
  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  IF(numpe==1) THEN;  WRITE(ch,'(I6.6)') numpe
    OPEN(12,file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',            &
         action='write')
    WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
    WRITE(12,'(A/A/A)') "part", "     1","coordinates"
  END IF
  ALLOCATE(disp_pp(nodes_pp*ndim),temp(nodes_pp)); disp_pp=zero; temp=zero
  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,               &
                     node_start,node_end,eld_pp,disp_pp,1)
  DO i=1,ndim ; temp=zero
    DO j=1,nodes_pp
      k=i+(ndim*(j-1)); temp(j)=disp_pp(k)
    END DO
    CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
  END DO; IF(numpe==1) CLOSE(12)

  timest(8) = elap_time()
  
  peak_memory_use = p_memory_use()

  IF(numpe==1)THEN
    WRITE(11,'(A,F10.2)') "Time to read input:       ",                        &
                          timest(2)-timest(1) + timest(6)-timest(5)
    WRITE(11,'(A,F10.2)') "Time for setup:           ",timest(3)-timest(2)
    WRITE(11,'(A,F10.2)') "Time for matrix assemble: ",timest(4)-timest(3)
    WRITE(11,'(A,F10.2)') "Time to solve equations:  ",                        &
                          timest(5)-timest(4) + timest(7)-timest(6)
    WRITE(11,'(A,F10.2)') "Time to write results:    ",timest(8)-timest(7)
    WRITE(11,'(A)')       "                          ----------"
    WRITE(11,'(A,F10.2)') "This analysis took:       ",elap_time()-timest(1)
    WRITE(11,*)
    WRITE(11,'(A,F10.2,A)') "Peak memory use: ",peak_memory_use," GB"
    WRITE(11,*)
    ! REVISION is substituted with a string like "2108" (including the double
    ! quotes) by the preprocessor, see the makefile.
    WRITE(11,'(2A,2(I0,A),4(F0.2,A),F0.2)')                                    &
      REVISION,tab,                                                            &
      npes,tab,                                                                &
      neq,tab,                                                                 &
      timest(3)-timest(2),tab,                                                 &
      timest(4)-timest(3),tab,                         &
      timest(5)-timest(4) + timest(7)-timest(6),tab,                           &
      timest(3)-timest(2) + timest(4)-timest(3)                                &
      + timest(5)-timest(4) + timest(7)-timest(6),tab,                         &
      peak_memory_use
  END IF
  
  !-----------------------------------------------------------------------------
  ! 11. Shut down PETSc and ParaEFM
  !-----------------------------------------------------------------------------
  
  IF (solvers == petsc_solvers) THEN
    CALL p_shutdown
  END IF

  CALL SHUTDOWN() 
  
END PROGRAM xx18
  
