PROGRAM xx12      
!------------------------------------------------------------------------------
!      Program XX.12 Three dimensional anallysis of conduction equation using 
!                    8-node hexahedral elements; pcg version implicit;  
!                    integration in time using 'theta' method parallel version
!
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

! neq,ntot are now global variables - not declared

 INTEGER, PARAMETER  :: ndim=3,nodof=1
 INTEGER             :: nod,nn,nr,nip
 INTEGER             :: i,j,k,l,iters,limit,iel
 INTEGER             :: nxe,nye,nze,neq_temp,nn_temp
 INTEGER             :: nstep,npri,nres,it,is,nlen
 INTEGER             :: node_end,node_start,nodes_pp
 INTEGER             :: loaded_freedoms,fixed_freedoms,loaded_nodes
 INTEGER             :: fixed_freedoms_pp,fixed_freedoms_start
 INTEGER             :: loaded_freedoms_pp,loaded_freedoms_start
 INTEGER             :: nels,ndof,ielpe,npes_pp
 INTEGER             :: argc,iargc,meshgen,partitioner
 INTEGER             :: np_types
 REAL(iwp)           :: aa,bb,cc,kx,ky,kz,det,theta,dtim,real_time
 !REAL(iwp)           :: val0 = 100.0_iwp
 REAL(iwp)           :: tol,alpha,beta,up,big,q
 REAL(iwp)           :: rho,cp,val0
 REAL(iwp),PARAMETER :: zero = 0.0_iwp,penalty=1.e20_iwp
 CHARACTER(LEN=15)   :: element
 CHARACTER(LEN=50)   :: fname,job_name,label
 CHARACTER(LEN=50)   :: program_name='xx12'
 LOGICAL             :: converged = .false.

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------ 

 REAL(iwp),ALLOCATABLE :: loads_pp(:),u_pp(:),p_pp(:),points(:,:),kay(:,:)
 REAL(iwp),ALLOCATABLE :: coord(:,:),fun(:),jac(:,:),der(:,:),deriv(:,:)
 REAL(iwp),ALLOCATABLE :: weights(:),d_pp(:),kc(:,:),pm(:,:),funny(:,:)
 REAL(iwp),ALLOCATABLE :: p_g_co_pp(:,:,:),storka_pp(:,:,:)
 REAL(iwp),ALLOCATABLE :: storkb_pp(:,:,:),x_pp(:),xnew_pp(:),pmul_pp(:,:)
 REAL(iwp),ALLOCATABLE :: utemp_pp(:,:),diag_precon_pp(:),diag_precon_tmp(:,:)
 REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),timest(:)
 REAL(iwp),ALLOCATABLE :: disp_pp(:),eld_pp(:,:)
 REAL(iwp),ALLOCATABLE :: val(:,:),val_f(:),store_pp(:)
 INTEGER,ALLOCATABLE   :: rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:),no(:)
 INTEGER,ALLOCATABLE   :: no_pp(:),no_f_pp(:),no_pp_temp(:),no_global(:)
 INTEGER,ALLOCATABLE   :: sense(:),node(:)
 
!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions
!------------------------------------------------------------------------------ 

  ALLOCATE(timest(25))
  timest    = zero
  timest(1) = elap_time()

  CALL find_pe_procs(numpe,npes)
  
  PRINT *, "FIND_PE_PROCS on processor ", numpe, " of ", npes

  argc = iargc()
  IF(argc /= 1) CALL job_name_error(numpe,program_name)
  CALL GETARG(1,job_name)

  CALL read_xx12(job_name,numpe,dtim,element,fixed_freedoms,kx,ky,kz,limit,   &
                 loaded_nodes,meshgen,nels,nip,nn,nod,npri,nr,nstep,          &
                 partitioner,theta,tol,np_types,rho,cp,val0)

  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)

  ndof = nod*nodof
  ntot = ndof

  ALLOCATE(g_num_pp(nod,nels_pp))
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  IF (nr>0) ALLOCATE(rest(nr,nodof+1))

  g_num_pp   	 = 0
  g_coord_pp 	 = zero
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
  
! nn_temp=0
! neq_temp=0
! nye=nels/nxe/nze
! nr=(nxe+1)*(nye+1)+(nxe+1)*nze+nye*nze
! ielpe=iel_start

!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------

 ALLOCATE (points(nip,ndim),weights(nip),kay(ndim,ndim),                      &
           coord(nod,ndim),fun(nod),jac(ndim,ndim),der(ndim,nod),g(ntot),     &
           deriv(ndim,nod),pm(ntot,ntot),                                     &
           kc(ntot,ntot),funny(1,nod),num(nod),                               &
           g_g_pp(ntot,nels_pp),storka_pp(ntot,ntot,nels_pp),                 &
           utemp_pp(ntot,nels_pp),storkb_pp(ntot,ntot,nels_pp),               &
           pmul_pp(ntot,nels_pp))

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
    DEALLOCATE(rest)
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
  
  IF(numpe==1) PRINT *, " *** Found number of equations in: ", timest(7)-timest(6)," s"

!------------------------------------------------------------------------------
! 6. Create interprocessor communication tables
!------------------------------------------------------------------------------

  CALL calc_neq_pp          
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl2(npes_pp,npes,g_g_pp)

  nres = nxe*(nze-1) + 1

  DO i = 1,neq_pp
    IF(nres==ieq_start+i-1) THEN
      it = numpe; is = i
      IF(numpe==1) PRINT *, " *** it = ", it, " is = ", i
    END IF
  END DO

  timest(8) = elap_time()

  IF(numpe==1) PRINT *, " *** Created ggl in: ", timest(8)-timest(7), " s"

!------------------------------------------------------------------------------
! 7. Allocate arrays dimensioned by neq_pp 
!------------------------------------------------------------------------------

  ALLOCATE(loads_pp(neq_pp),diag_precon_pp(neq_pp),u_pp(neq_pp),              &
           d_pp(neq_pp),p_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp))

  loads_pp  = zero ; diag_precon_pp = zero ; u_pp = zero
  d_pp      = zero ; p_pp           = zero ; x_pp = zero ; xnew_pp = zero

  timest(9) = elap_time()

  IF(numpe==1) PRINT *, " *** Allocated arrays dimensioned by neq_pp in: ",   &
                          timest(9)-timest(8), " s"  

  
! kay=0.0_iwp
! kay(1,1)=kx
! kay(2,2)=ky
! kay(3,3)=kz
! CALL box_bc8(nxe,nye,nze,rest)

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------

 CALL sample(element,points,weights)

 storka_pp = zero 
 storkb_pp = zero
 kay       = zero
 kay(1,1)  = kx
 kay(2,2)  = ky
 kay(3,3)  = kz
 
 elements_3: DO iel=1,nels_pp

   kc = zero ; pm = zero

   gauss_pts: DO i=1,nip
     CALL shape_der(der,points,i)
     CALL shape_fun(fun,points,i)
     funny(1,:) = fun(:)
     jac        = MATMUL(der,g_coord_pp(:,:,iel))
     det        = determinant(jac)
     CALL invert(jac)
     deriv      = MATMUL(jac,der)
     kc         = kc + MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
     pm         = pm + MATMUL(TRANSPOSE(funny),funny)*det*weights(i)*rho*cp
   END DO gauss_pts

! storka & storkb to be equivalent to (3.94) S&G_4ed pp.94
! ([Mm] + theta*dt[Kc]){phi}_1 = ([Mm]-(1-theta)dt[Kc]){phi}_0 + theta*dt{Q}_1+(1-theta)dt{Q}_0
! storka = ([Mm] + theta*dt[Kc])
! storkb = ([Mm]-(1-theta)dt[Kc])
! Heat Equation must be rearranged into same format as (3.90)
! -k*phi + cp*rho*dphi/dt = q
! [Kc}{phi} + [Mm]{dphi/dt} = {q}	(3.90)
!
! Therefore
! Kc multiplied by -kay	!---- However, applying this below causes solution to explode
! Mm multiplied by cp*rho
! 
   storka_pp(:,:,iel)=pm+kc*theta*dtim
   storkb_pp(:,:,iel)=pm-kc*(1._iwp-theta)*dtim
   
 END DO elements_3

 timest(10) = elap_time()

 IF(numpe==1) PRINT *, " *** 8. Element stiffness integration and storage in: ",   &
                           timest(10)-timest(9), " s"  

 
!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------ 
  
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))
  diag_precon_tmp = zero

  elements_4: DO iel = 1,nels_pp 
    DO k=1,ntot
      diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)+storka_pp(k,k,iel)
    END DO
  END DO elements_4
 
 CALL scatter(diag_precon_pp,diag_precon_tmp)

! DEALLOCATE(diag_precon_tmp)

! diag_precon_pp=1._iwp/diag_precon_pp ! needs moving - MOVED to sec.11d
 
 timest(11) = elap_time()

 IF(numpe==1) PRINT *, " *** 9. Build the diagonal preconditioner in: ",   &
                          timest(11)-timest(10), " s" 

!------------------------------------------------------------------------------
! 10. Allocate disp_pp array and open file to write temperature output
!------------------------------------------------------------------------------

 !IF(numpe==it)THEN
 IF(numpe==1)THEN
   fname = job_name(1:INDEX(job_name, " ")-1) // ".res"
   OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
 END IF

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  ALLOCATE(disp_pp(nodes_pp))
  ALLOCATE(eld_pp(ntot,nels_pp))
  
  IF(numpe==1) THEN
    fname   = job_name(1:INDEX(job_name, " ")-1)//".ttr"
    OPEN(24, file=fname, status='replace', action='write')
    label   = "*TEMPERATURE"  
  END IF
  
!------------------------------------------------------------------------------
! 11a. Read in the initial conditions and assign to equations
!------------------------------------------------------------------------------
 
  j=1
 
! Started time-stepping loop here to include setting of fixed_freedoms
! and loaded_freedoms at start of each step
  timesteps: DO j=1,nstep

! Reassign the diagonal precondtioner at start of each time step (changes in sec.11e)
! Done using diag_precon_tmp - reason for commenting out deallocation of variable in sec.9
  IF(j>1)THEN
     CALL scatter(diag_precon_pp,diag_precon_tmp)
!     diag_precon_pp=1._iwp/diag_precon_pp
  END IF
 
!  val0 = 100.0_iwp
!  loads_pp = val0    ! needs to be read in from file
  pmul_pp  = .0_iwp
  utemp_pp = zero
  
!------------------------------------------------------------------------------
! 11b. Read in fixed nodal displacements and assign to equations
!------------------------------------------------------------------------------

! In sec.11b & sec.11c several variables are allocated
! To avoid reallocation whilst looping they are deallocated at end of section
! These could be taken outside time-stepping loop

  IF(fixed_freedoms > 0) THEN

    ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),                         &
             no_pp_temp(fixed_freedoms),sense(fixed_freedoms))
    ALLOCATE(val_f(fixed_freedoms),no_global(fixed_freedoms))

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
    DEALLOCATE(no_global)

  END IF

  IF(fixed_freedoms == 0) fixed_freedoms_pp = 0
                           
!------------------------------------------------------------------------------
! 11c. Read in loaded nodes and get starting r_pp
!------------------------------------------------------------------------------

! timesteps: DO j=1,nstep	! - Gives same results as if started in sec.11a

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
      loads_pp(no_pp(i) - ieq_start + 1) = val(loaded_freedoms_start + i - 1)
    END DO

!   CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))

    q = SUM_P(loads_pp) ! Not really needed - sums loads

    DEALLOCATE(node,val)
    DEALLOCATE(no_pp)

  END IF

  !DEALLOCATE(g_g_pp)

  timest(12) = elap_time()
!  IF(numpe==1) PRINT *, " *** 10. Read in the initial conditions and assign to equations in: ",   &
!                           timest(12)-timest(11), " s"

!------------------------------------------------------------------------------
! 11d. Invert the preconditioner. 
!     If there are fixed freedoms, first apply a penalty
!------------------------------------------------------------------------------

! timesteps: DO j=1,nstep	! - starting here allows solution to diffuse to zero

  IF(fixed_freedoms_pp > 0) THEN
    DO i = 1,fixed_freedoms_pp
       l =  no_f_pp(i) - ieq_start + 1
       diag_precon_pp(l) = diag_precon_pp(l) + penalty
       store_pp(i)       = diag_precon_pp(l)
    END DO
  END IF

  diag_precon_pp = 1._iwp/diag_precon_pp

!------------------------------------------------------------------------------
! 11e. Initiallize preconditioned conjugate gradient
!------------------------------------------------------------------------------

  IF(fixed_freedoms_pp>0) THEN
    DO i = 1, fixed_freedoms_pp
      l       = no_f_pp(i) - ieq_start + 1
      k       = fixed_freedoms_start + i - 1
      loads_pp(l) = store_pp(i) * val_f(k) ! As above (section 10c)
    END DO
    
!    DEALLOCATE(no_f_pp,store_pp,val_f)
  END IF

!------------------------------------------------------------------------------
! 12. Time-stepping loop
!------------------------------------------------------------------------------
 
! timesteps: DO j=1,nstep ! - starting here allows solution to diffuse to zero
 
   real_time = j*dtim
   u_pp      = zero

   CALL gather(loads_pp,pmul_pp)
   elements_5: DO iel=1,nels_pp
     utemp_pp(:,iel)=MATMUL(storkb_pp(:,:,iel),pmul_pp(:,iel))
   END DO elements_5
   CALL scatter(u_pp,utemp_pp)

   loads_pp=u_pp

  timest(13) = elap_time()

!  IF(numpe==1) PRINT *, " *** 11. Time-stepping loop in: ",   &
!                           timest(13)-timest(12), " s"
 
!------------------------------------------------------------------------------
! 13. Solve simultaneous equations by pcg
!------------------------------------------------------------------------------

   d_pp = diag_precon_pp*loads_pp
   p_pp = d_pp
   x_pp = zero

   iters = 0

   iterations: DO 

     iters   = iters+1
     u_pp    = zero
     pmul_pp = zero
     utemp_pp = zero

     CALL gather(p_pp,pmul_pp)  
     elements_6: DO iel=1,nels_pp
       utemp_pp(:,iel)=MATMUL(storka_pp(:,:,iel),pmul_pp(:,iel)) 
     END DO elements_6
     CALL scatter(u_pp,utemp_pp)

  timest(14) = elap_time()

!  IF(numpe==1) PRINT *, " *** 12. Solve simultaneous equations by pcg in: ",   &
!                           timest(14)-timest(13), " s"



!------------------------------------------------------------------------------
! 14. PCG equation solution
!------------------------------------------------------------------------------

! Copied directly from xx11, not sure if needed
     IF(fixed_freedoms_pp > 0) THEN
       DO i = 1, fixed_freedoms_pp
         l       = no_f_pp(i) - ieq_start + 1
         u_pp(l) = p_pp(l) * store_pp(i)
       END DO
     END IF

     up       = DOT_PRODUCT_P(loads_pp,d_pp)
     alpha    = up/DOT_PRODUCT_P(p_pp,u_pp)
     xnew_pp  = x_pp+p_pp*alpha
     loads_pp = loads_pp-u_pp*alpha
     d_pp     = diag_precon_pp*loads_pp
     beta     = DOT_PRODUCT_P(loads_pp,d_pp)/up
     p_pp     = d_pp+p_pp*beta
     u_pp     = xnew_pp

     CALL checon_par(xnew_pp,tol,converged,x_pp)
     IF(converged.OR.iters==limit)EXIT

   END DO iterations

   IF(fixed_freedoms_pp > 0) THEN
     DEALLOCATE(no_f_pp,store_pp,val_f)
   END IF

   loads_pp=xnew_pp
   
!   utemp_pp = zero
   eld_pp   = zero
   disp_pp  = zero
!   CALL gather(loads_pp(1:),utemp_pp)
   CALL gather(xnew_pp(1:),eld_pp)
   
   CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,              &
!                      node_start,node_end,utemp_pp,disp_pp,1)
                      node_start,node_end,eld_pp,disp_pp,1)
   IF(j/npri*npri==j.AND.numpe==1)THEN
     CALL write_nodal_variable(label,24,j,nodes_pp,npes,numpe,nodof,disp_pp)
  
     IF(numpe==1)THEN
       WRITE(11,'(E12.4,8E19.8)')real_time,disp_pp
     END IF
   END IF

 END DO timesteps
 
 timest(15) = elap_time()

 IF(numpe==1) PRINT *, " *** 13. PCG equation solution in: ",   &
                          timest(15)-timest(14), " s"
 
 IF(numpe==1) PRINT *, "This analysis took  :",elap_time()-timest(1)
 
 IF(numpe==1)THEN
   CLOSE(11)
   CLOSE(24)
 END IF

 CALL shutdown() 

END PROGRAM xx12
