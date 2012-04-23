PROGRAM xx11         
!------------------------------------------------------------------------------
!      program xx11 three dimensional analysis of Laplace's equation
!      using 8-node brick elements, preconditioned conjugate gradient solver
!      only integrate one element , diagonal preconditioner diag_precon
!      parallel version  ;   central loaded or fixed freedom   ;  box_bc  
!------------------------------------------------------------------------------

  USE precision  ; USE global_variables ; USE mp_interface
  USE input      ; USE output           ; USE loading
  USE timing     ; USE maths            ; USE gather_scatter
  USE partition  ; USE elements         ; USE steering       
  USE geometry   ; USE pcg

  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

 !  neq , ntot  are now global variables - not declared

  INTEGER, PARAMETER    :: ndim=3,nodof=1
  INTEGER               :: nod,nn,nr,nip
  INTEGER               :: i,j,k,iters,limit,iel,num_no,no_index_start
  INTEGER               :: node_end,node_start,nodes_pp
  INTEGER               :: loaded_freedoms,fixed_freedoms,loaded_nodes
  INTEGER               :: fixed_freedoms_pp,fixed_freedoms_start
  INTEGER               :: loaded_freedoms_pp,loaded_freedoms_start
  INTEGER               :: nels,ndof,ielpe,npes_pp
  INTEGER               :: argc,iargc,meshgen,partitioner
  REAL(iwp)             :: kx,ky,kz,det,tol,up,alpha,beta,q
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp,penalty=1.e20_iwp
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: fname,job_name,label
  CHARACTER(LEN=50)     :: program_name='xx11'
  LOGICAL               :: converged 

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),kc(:,:),coord(:,:), weights(:),temp_pp(:)
  REAL(iwp),ALLOCATABLE :: p_g_co_pp(:,:,:), jac(:,:), der(:,:), deriv(:,:)
  REAL(iwp),ALLOCATABLE :: col(:,:),row(:,:),kcx(:,:),kcy(:,:),kcz(:,:)
  REAL(iwp),ALLOCATABLE :: diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:)
  REAL(iwp),ALLOCATABLE :: xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:)
  REAL(iwp),ALLOCATABLE :: d_pp(:),diag_precon_tmp(:,:),eld_pp(:,:),val(:,:),val_f(:)
  REAL(iwp),ALLOCATABLE :: store_pp(:),storkc_pp(:,:,:),eld(:),timest(:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:)
  INTEGER, ALLOCATABLE  :: rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:),no(:)
  INTEGER, ALLOCATABLE  :: no_f(:),no_local_temp(:),no_local_temp_f(:)
  INTEGER, ALLOCATABLE  :: no_local(:),no_pp(:),no_f_pp(:),no_pp_temp(:),no_global(:)
  INTEGER, ALLOCATABLE  :: sense(:),node(:)
 
!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions
!------------------------------------------------------------------------------

  ALLOCATE(timest(25))
  timest    = zero
  timest(1) = elap_time()

  CALL find_pe_procs(numpe,npes)

  argc = iargc()
  IF(argc /= 1) CALL job_name_error(numpe,program_name)
  CALL GETARG(1,job_name)

  CALL read_p123(job_name,numpe,element,fixed_freedoms,kx,ky,kz,limit,        &
                 loaded_nodes,meshgen,nels,nip,nn,nod,nr,partitioner,tol)

  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)

  ndof = nod*nodof 
  ntot = ndof 

  ALLOCATE(g_num_pp(nod,nels_pp))
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  IF (nr>0) ALLOCATE(rest(nr,nodof+1))

  g_num_pp       = 0
  g_coord_pp     = zero
  IF (nr>0) rest = 0

  timest(2) = elap_time()

  CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)
  timest(3) = elap_time()

  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  timest(4) = elap_time()

  CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  timest(5) = elap_time()

  IF (nr>0) CALL read_rest(job_name,numpe,rest)
  timest(6) = elap_time()

  IF(numpe==1) PRINT *, " *** Read input data in: ", timest(6)-timest(1)," s"

!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------
  
  ALLOCATE (points(nip,ndim),coord(nod,ndim),jac(ndim,ndim),                  &
            der(ndim,nod),deriv(ndim,nod),kcx(ntot,ntot),weights(nip),        &
            g(ntot),pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),col(ntot,1), &
            num(nod),g_g_pp(ntot,nels_pp),kcy(ntot,ntot),                     &
            no_local_temp(1),row(1,ntot),                                     &
            eld(ntot),no_f(1),no_local_temp_f(1),kcz(ntot,ntot),              &
            storkc_pp(ntot,ntot,nels_pp))

  IF(numpe==1) PRINT *, " *** Allocated dynamic arrays in: ",                 &
                          elap_time()-timest(6)," s"

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------
      
  IF (nr>0) CALL rearrange_2(rest)  

  g_g_pp = 0

  ! When nr = 0, g_num_pp and g_g_pp are identical
  IF(nr>0) THEN
    elements_1: DO iel = 1, nels_pp
      CALL find_g4(g_num_pp(:,iel),g_g_pp(:,iel),rest)
    END DO elements_1
  ELSE
    g_g_pp = g_num_pp
  END IF
  
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
  CALL make_ggl2(npes_pp,npes,g_g_pp)

  timest(8) = elap_time()

  IF(numpe==1) PRINT *, " *** Created ggl in: ", timest(8)-timest(7), " s"

!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------

  ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),            &
           u_pp(neq_pp),diag_precon_pp(neq_pp),d_pp(neq_pp))

  r_pp           = zero ; p_pp     = zero ; x_pp = zero ; xnew_pp = zero
  diag_precon_pp = zero

  timest(9) = elap_time()

  IF(numpe==1) PRINT *, " *** Allocated arrays dimensioned by neq_pp in: ",   &
                          timest(9)-timest(8), " s"

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------

  CALL sample(element,points,weights)
   
  storkc_pp = zero

  elements_3: DO iel=1,nels_pp
 
    kcx = zero; kcy = zero; kcz = zero

    gauss_pts_1:  DO i=1,nip
      CALL shape_der (der,points,i)
      jac      = MATMUL(der,g_coord_pp(:,:,iel))
      det      = determinant(jac)
      CALL invert(jac)
      deriv    = MATMUL(jac,der)
      row(1,:) = deriv(1,:)
      eld      = deriv(1,:)
      col(:,1) = eld
      kcx      = kcx + MATMUL(col,row)*det*weights(i)
      row(1,:) = deriv(2,:)
      eld      = deriv(2,:)
      col(:,1) = eld
      kcy      = kcy + MATMUL(col,row)*det*weights(i)
      row(1,:) = deriv(3,:)
      eld      = deriv(3,:)
      col(:,1) = eld
      kcz      = kcz + MATMUL(col,row)*det*weights(i)
    END DO gauss_pts_1        

    storkc_pp(:,:,iel) = kcx*kx + kcy*ky + kcz*kz 
    
  END DO elements_3

  timest(10) = elap_time()

  IF(numpe==1) PRINT *, " *** Computed element stiffness matrices in: ",    &
                          timest(10)-timest(9), " s"

!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------

  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  diag_precon_tmp = zero

  elements_4: DO iel = 1,nels_pp
    DO i = 1,ndof
      diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + storkc_pp(i,i,iel)
    END DO
  END DO elements_4

  CALL scatter(diag_precon_pp,diag_precon_tmp)

  DEALLOCATE(diag_precon_tmp)

  timest(11) = elap_time()

  IF(numpe==1) PRINT *, " *** Built diagonal preconditioner in: ",            &
                          timest(11)-timest(10), " s"
  
!------------------------------------------------------------------------------
! 10. Read in fixed nodal displacements and assign to equations
!------------------------------------------------------------------------------

  IF(fixed_freedoms > 0) THEN

    ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),val_f(fixed_freedoms),   &
             no_pp_temp(fixed_freedoms),sense(fixed_freedoms),                &
             no_global(fixed_freedoms))

    node  = 0 ; no = 0 ; no_pp_temp = 0 ; sense = 0 ; no_global = 0
    val_f = zero

    CALL read_fixed(job_name,numpe,node,sense,val_f)
    CALL find_no2(g_g_pp,g_num_pp,node,sense,fixed_freedoms_pp,               &
                  fixed_freedoms_start,no)
    CALL MPI_ALLREDUCE(no,no_global,fixed_freedoms,MPI_INTEGER,MPI_MAX,       &
                       MPI_COMM_WORLD,ier)
    CALL reindex_fixed_nodes(ieq_start,no_global,no_pp_temp,                  &
                             fixed_freedoms_pp,fixed_freedoms_start,neq_pp)

    ALLOCATE(no_f_pp(fixed_freedoms_pp),store_pp(fixed_freedoms_pp))

    no_f_pp    = 0
    store_pp = zero
    no_f_pp    = no_pp_temp(1:fixed_freedoms_pp)

    DEALLOCATE(node,no,sense,no_pp_temp)

  END IF

  IF(fixed_freedoms == 0) fixed_freedoms_pp = 0

  IF (nr>0) DEALLOCATE(rest)
  
  timest(12) = elap_time()
  
  IF(numpe==1) PRINT *, " *** Applied fixed freedoms in:",                             &
                          timest(12)-timest(11), " s"

!------------------------------------------------------------------------------
! 11. Read in loaded nodes and get starting r_pp
!------------------------------------------------------------------------------

  loaded_freedoms = loaded_nodes ! hack

  IF(loaded_freedoms > 0) THEN

    ALLOCATE(node(loaded_freedoms),val(nodof,loaded_freedoms))
    ALLOCATE(no_pp_temp(loaded_freedoms))

    val = zero ; node = 0

    CALL read_loads(job_name,numpe,node,val)
    CALL reindex_fixed_nodes(ieq_start,node,no_pp_temp,                       &
                             loaded_freedoms_pp,loaded_freedoms_start,neq_pp)

    ALLOCATE(no_pp(loaded_freedoms_pp))

    no_pp    = no_pp_temp(1:loaded_freedoms_pp)

    DEALLOCATE(no_pp_temp)

    DO i = 1, loaded_freedoms_pp
      r_pp(no_pp(i) - ieq_start + 1) = val(loaded_freedoms_start + i - 1)
    END DO

!   CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))

    q = SUM_P(r_pp)

    DEALLOCATE(node,val)

  END IF

  DEALLOCATE(g_g_pp)

  timest(12) = elap_time()

  IF(numpe==1) PRINT *, " *** Applied loads in:",                             &
                          timest(12)-timest(11), " s"

!------------------------------------------------------------------------------
! 12. Invert the preconditioner. 
!     If there are fixed freedoms, first apply a penalty
!------------------------------------------------------------------------------

  IF(fixed_freedoms_pp > 0) THEN
    DO i = 1,fixed_freedoms_pp
       j =  no_f_pp(i) - ieq_start + 1
       diag_precon_pp(j) = diag_precon_pp(j) + penalty
       store_pp(i)       = diag_precon_pp(j)
    END DO
  END IF

  diag_precon_pp = 1._iwp/diag_precon_pp

!------------------------------------------------------------------------------
! 13. Initiallize preconditioned conjugate gradient
!------------------------------------------------------------------------------

  IF(fixed_freedoms_pp>0) THEN
    DO i = 1, fixed_freedoms_pp
      j       = no_f_pp(i) - ieq_start + 1
      k       = fixed_freedoms_start + i - 1
      r_pp(j) = store_pp(i) * val_f(k)
    END DO
  END IF
 
  d_pp = diag_precon_pp*r_pp
  p_pp = d_pp
  x_pp = zero

  IF(numpe==1) PRINT *, " *** Initiallized preconditioned conjugate gradients"

!------------------------------------------------------------------------------
! 14. Preconditioned conjugate gradient iterations
!------------------------------------------------------------------------------

    iters = 0

    iterations  :      DO 

      iters    = iters + 1
      u_pp     = zero
      pmul_pp  = zero
      utemp_pp = zero 

      CALL gather(p_pp,pmul_pp) 

      elements_5 : DO iel = 1, nels_pp
        utemp_pp(:,iel) = MATMUL(storkc_pp(:,:,iel),pmul_pp(:,iel)) 
      END DO elements_5  

      CALL scatter(u_pp,utemp_pp)

      IF(fixed_freedoms_pp > 0) THEN
        DO i = 1, fixed_freedoms_pp
          j       = no_f_pp(i) - ieq_start + 1
          u_pp(j) = p_pp(j) * store_pp(i)
        END DO
      END IF

      up      = DOT_PRODUCT_P(r_pp,d_pp)
      alpha   = up/ DOT_PRODUCT_P(p_pp,u_pp)
      xnew_pp = x_pp + p_pp* alpha 
      r_pp    = r_pp - u_pp*alpha
      d_pp    = diag_precon_pp*r_pp 
      beta    = DOT_PRODUCT_P(r_pp,d_pp)/up
      p_pp    = d_pp+p_pp*beta    

      CALL checon_par(xnew_pp,tol,converged,x_pp)    

      IF(converged .OR. iters==limit) EXIT

    END DO iterations

    timest(13) = elap_time()
    
    IF(numpe==1)THEN
        WRITE(11,'(A,I5)')"The number of iterations to convergence was  ",iters 
        WRITE(11,'(A,E12.4)')"The total load is                            ",q
    END IF
    
    IF(nels==1000000)THEN
      IF(numpe==1)THEN
    
        WRITE(11,'(A)')   "The  potentials are   :"
        WRITE(11,'(A)') "   Freedom       Potential"
  
        WRITE(11,'(A,E12.4)') "9901     ", xnew_pp(9901)
        WRITE(11,'(A,E12.4)') "9902     ", xnew_pp(9902)
        WRITE(11,'(A,E12.4)') "9903     ", xnew_pp(9903)
        WRITE(11,'(A,E12.4)') "9904     ", xnew_pp(9904)

      END IF
    END IF
  IF(numpe==1) WRITE(11,*) "This analysis took   ", elap_time( ) - timest(1)

  timest(14) = elap_time()

!------------------------------------------------------------------------------
! 15. Print out results
!------------------------------------------------------------------------------
  
!  IF(numpe==1)THEN
!    WRITE(11,*)
!    WRITE(11,'(/A)')"  Node Nodal Temp"
!    DO k=1,nn
!       WRITE(11,'(I5,1E12.4)')k,xnew_pp(k)
!    END DO
!  END IF
  
  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
    
  IF(numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ")-1)//".tmp"
    OPEN(24, file=fname, status='replace', action='write')
  END IF
  
!------------------------------------------------------------------------------
! 16a. Temperatures
!------------------------------------------------------------------------------

  !IF(numpe==1) PRINT *, "xnew_pp:"
  !IF(numpe==1) PRINT *, xnew_pp
  !IF(numpe==2) PRINT *, xnew_pp
  !IF(numpe==3) PRINT *, xnew_pp
  
  ALLOCATE(eld_pp(ntot,nels_pp))
  eld_pp = zero
  CALL gather(xnew_pp(1:),eld_pp)
  DEALLOCATE(xnew_pp)

  !IF(numpe==1) PRINT *, "eld_pp:"
  !IF(numpe==1) PRINT *, eld_pp
  !IF(numpe==2) PRINT *, eld_pp
  !IF(numpe==3) PRINT *, eld_pp

  ALLOCATE(temp_pp(nodes_pp))
  temp_pp = zero

  label   = "*TEMPERATURE"

  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,              &
                     node_start,node_end,eld_pp,temp_pp,1)
  CALL write_nodal_variable(label,24,1,nodes_pp,npes,numpe,ndim,temp_pp)

  DEALLOCATE(temp_pp)

  IF(numpe==1) CLOSE(24)  

!------------------------------------------------------------------------------
! 16e. Von Mises stress (rho_v)
!      rho_v = sqrt( ( (rho1-rho2)^2 + (rho2-rho3)^2 + (rho1-rho3)^2 ) / 2 )
!------------------------------------------------------------------------------
  
!  label = "*MISES STRESS"
!  
!  DO i = 1,nodes_pp
!    j = ((i-1)*nodof)+1
!    k = j + 1
!    l = j + 2
!    princinodes_pp(j) = SQRT(((princinodes_pp(j)-princinodes_pp(k)) **2 +     &
!                              (princinodes_pp(k)-princinodes_pp(l)) **2 +     &
!                              (princinodes_pp(j)-princinodes_pp(l)) **2)      &
!                              * 0.5_iwp)
!    princinodes_pp(k:l) = zero
!  END DO
!
!  CALL write_nodal_variable(label,27,1,nodes_pp,npes,numpe,nodof,             &
!                            princinodes_pp)
!                            
!  DEALLOCATE(princinodes_pp)
!
!  IF(numpe==1) CLOSE(27)
 
!------------------------------------------------------------------------------
! 17. Output performance data
!------------------------------------------------------------------------------

  CALL WRITE_P123(fixed_freedoms,iters,job_name,loaded_freedoms,neq,nn,npes,  &
                  nr,numpe,timest,q)

  CALL shutdown()

END PROGRAM xx11
