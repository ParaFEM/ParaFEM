PROGRAM p123         
!------------------------------------------------------------------------------
!      program 12.3 three dimensional analysis of Laplace's equation
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

  INTEGER, PARAMETER    :: ndim=3,nodof=1,nod=8,partitioner=1
  INTEGER               :: nn,nr,nip
  INTEGER               :: nxe,nye,nze,nres,is,it
  INTEGER               :: nn_temp,neq_temp
  INTEGER               :: i,j,k,iters,limit,iel,num_no,no_index_start
  INTEGER               :: loaded_freedoms,fixed_freedoms
  INTEGER               :: nels,ndof,ielpe,npes_pp
  REAL(iwp)             :: aa,bb,cc 
  REAL(iwp)             :: kx,ky,kz,det,tol,up,alpha,beta,q
  REAL(iwp),PARAMETER   :: zero = 0.0_iwp,penalty=1.e20_iwp
  CHARACTER(LEN=15)     :: element= 'hexahedron'
  CHARACTER(LEN=50)     :: job_name
  LOGICAL               :: converged 

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),kc(:,:),coord(:,:), weights(:)
  REAL(iwp),ALLOCATABLE :: p_g_co_pp(:,:,:), jac(:,:), der(:,:), deriv(:,:)
  REAL(iwp),ALLOCATABLE :: col(:,:),row(:,:),kcx(:,:),kcy(:,:),kcz(:,:)
  REAL(iwp),ALLOCATABLE :: diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:)
  REAL(iwp),ALLOCATABLE :: xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:)
  REAL(iwp),ALLOCATABLE :: d_pp(:),diag_precon_tmp(:,:),val(:),val_f(:)
  REAL(iwp),ALLOCATABLE :: store_pp(:),eld(:),timest(:)
  INTEGER, ALLOCATABLE  :: rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:),no(:)
  INTEGER, ALLOCATABLE  :: no_f(:),no_local_temp(:),no_local_temp_f(:)
  INTEGER, ALLOCATABLE  :: no_local(:)
 
!------------------------------------------------------------------------------
! 3. Read job_name from the command line. 
!    Read control data, mesh data, boundary and loading conditions
!------------------------------------------------------------------------------

  ALLOCATE(timest(25))
  timest = zero
  timest(1) = elap_time()

  CALL find_pe_procs(numpe,npes)

! argc = iargc()
! IF(argc /= 1) CALL job_name_error(numpe,program_name)
! CALL GETARG(1,job_name)

! CALL read_p123(job_name,numpe,fixed_freedoms,kx,ky,kz,limit,loaded_nodes,   &
!                meshgen,nels,nip,nn,nod,nr,partitioner,tol,v)

!-- remove

  IF (numpe==npes) THEN
   OPEN (10,FILE='p123.dat',STATUS=    'OLD',ACTION='READ')
   READ (10,*) nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz, tol,limit ,             &
               loaded_freedoms,fixed_freedoms
  END IF
  CALL bcast_inputdata_p123(numpe,npes,nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz, &
                             tol,limit,loaded_freedoms,fixed_freedoms)
!-- end remove

  CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)

  ndof = nod*nodof 
  ntot = ndof 

!-- remove

  nye      = nels/nxe/nze  
  neq_temp = 0 
  nn_temp  = 0 
  nr       = (nxe+1)*(nye+1)+(nxe+1)*nze+nye*nze  

!-- end remove

  ALLOCATE(g_num_pp(nod,nels_pp))
! ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  ALLOCATE(p_g_co_pp(nod,ndim,nels_pp))
  ALLOCATE(rest(nr,nodof+1))

  g_num_pp   = 0
! g_coord_pp = zero
  p_g_co_pp  = zero
  rest       = 0

  timest(2) = elap_time()

! CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)
  timest(3) = elap_time()

! IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  timest(4) = elap_time()

! CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
  timest(5) = elap_time()

! CALL read_rest(job_name,numpe,rest)
  timest(6) = elap_time()

  IF(numpe==1) PRINT *, " *** Read input data in: ", timest(6)-timest(1)," s"

!------------------------------------------------------------------------------
! 4. Allocate dynamic arrays used in main program
!------------------------------------------------------------------------------
  
  ALLOCATE (points(nip,ndim),coord(nod,ndim),jac(ndim,ndim),kc(ntot,ntot),    &
            der(ndim,nod),deriv(ndim,nod),kcx(ntot,ntot),weights(nip),        &
            g(ntot),pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),col(ntot,1), &
            num(nod),g_g_pp(ntot,nels_pp),no(1),kcy(ntot,ntot),val(1),        &
            no_local_temp(1),row(1,ntot),diag_precon_tmp(ntot,nels_pp),       &
            val_f(1),eld(ntot),no_f(1),no_local_temp_f(1),kcz(ntot,ntot))

  IF(numpe==1) PRINT *, " *** Allocated dynamic arrays in: ",                 &
                          elap_time()-timest(6)," s"


  CALL box_bc8(rest,nxe,nye,nze)
  ielpe = iel_start
  nres  = nxe*(nze-1)+1
  IF(loaded_freedoms>0) THEN ; no  = nres ;   val = 10.  ; END IF
  IF(fixed_freedoms>0)  THEN ; no_f= nres ; val_f = 100. ; END IF 

  CALL sample(element,points,weights)

!------------------------------------------------------------------------------
! 5. Loop the elements to find the steering array and the number of equations
!    to solve.
!------------------------------------------------------------------------------
      
  CALL rearrange_2(rest)  

  g_g_pp = 0

! elements_1: DO iel = 1, nels_pp
!   CALL find_g4(g_num_pp(:,iel),g_g_pp(:,iel),rest)
! END DO elements_1
!   
! neq = 0

!-- remove

  elements_0: DO iel = 1 , nels_pp
    CALL geometry_8bxz(ielpe,nxe,nze,aa,bb,cc,coord,num)
    CALL find_g4(num,g,rest) ; g_num_pp(:,iel) = num
    p_g_co_pp(:,:,iel) = coord; g_g_pp(:,iel) = g; ielpe = ielpe+1 
    i = MAXVAL(g); j = MAXVAL(num)
    IF(i>neq_temp)neq_temp = i; IF(j>nn_temp)nn_temp = j
  END DO elements_0
  neq = reduce(neq_temp)  ;  nn = reduce(nn_temp)

!-- remove
   
! elements_2: DO iel = 1, nels_pp
!   i = MAXVAL(g_g_pp(:,iel))
!   IF(i > neq) neq = i
! END DO elements_2
!
! neq = MAX_INTEGER_P(neq)

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
           u_pp(neq_pp),diag_precon_pp(neq_pp),d_pp(neq_pp),store_pp(neq_pp))

  r_pp           = zero ; p_pp     = zero ; x_pp = zero ; xnew_pp = zero
  diag_precon_pp = zero ; store_pp = zero

  timest(9) = elap_time()

  IF(numpe==1) PRINT *, " *** Allocated arrays dimensioned by neq_pp in: ",   &
                          timest(9)-timest(8), " s"

!------------------------------------------------------------------------------
! 8. Element stiffness integration and storage
!------------------------------------------------------------------------------

!-- remove

   diag_precon_tmp = zero

   DO i=1,neq_pp; IF(nres==ieq_start+i-1) THEN; it = numpe;is = i; END IF
   END DO

   iel=1;CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,num)
   kcx = .0_iwp; kcy = .0_iwp; kcz = .0_iwp
   gauss_pts_1:  DO i=1,nip
     CALL shape_der (der,points,i) ; jac = MATMUL(der,coord)
     det=determinant(jac);CALL invert(jac);deriv = MATMUL(jac,der)
     row(1,:) = deriv(1,:); eld=deriv(1,:); col(:,1) = eld
     kcx = kcx + MATMUL(col,row)*det*weights(i)
     row(1,:) = deriv(2,:); eld=deriv(2,:); col(:,1) = eld
     kcy = kcy + MATMUL(col,row)*det*weights(i)
     row(1,:) = deriv(3,:); eld=deriv(3,:); col(:,1) = eld
     kcz = kcz + MATMUL(col,row)*det*weights(i)
    END DO gauss_pts_1        ; kc = kcx*kx + kcy*ky + kcz*kz 
    elements_1: DO iel = 1,nels_pp  
    DO k=1,ntot;diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)+kc(k,k);END DO    
    END DO elements_1
    CALL scatter(diag_precon_pp,diag_precon_tmp);DEALLOCATE(diag_precon_tmp)
     IF(numpe==it)THEN
      OPEN (11,FILE='p123.res',STATUS='REPLACE',ACTION='WRITE')              
      WRITE(11,'(A,I5,A)') "This job ran on ", npes , "  processors" 
      WRITE(11,'(A)') "Global coordinates and node numbers "
      DO i = 1 , nels_pp,nels_pp-1                                         
                 WRITE(11,'(A,I8)')"Element ",i ; num = g_num_pp(:,i)
      DO k = 1,nod;WRITE(11,'(A,I8,3E12.4)')                               &
                 "  Node",num(k),p_g_co_pp(k,:,i); END DO
      END DO
      WRITE(11,'(A,3(I8,A))') "There are ",nn," nodes",nr," restrained and",&
                           neq," equations"
      WRITE(11,*) "Time after setup is   :", elap_time( ) - timest(1)
     END IF

!-- end remove
 
! CALL sample(element,points,weights)
   
! kcx = zero; kcy = zero; kcz = zero

! elements_3: DO i=1,nels_pp
!   gauss_pts_1:  DO i=1,nip
!     CALL shape_der (der,points,i)
!     jac      = MATMUL(der,g_coord_pp(:,:,iel))
!     det      = determinant(jac)
!     CALL invert(jac)
!     deriv    = MATMUL(jac,der)
!     row(1,:) = deriv(1,:)
!     eld      = deriv(1,:)
!     col(:,1) = eld
!     kcx      = kcx + MATMUL(col,row)*det*weights(i)
!     row(1,:) = deriv(2,:)
!     eld      = deriv(2,:)
!     col(:,1) = eld
!     kcy      = kcy + MATMUL(col,row)*det*weights(i)
!     row(1,:) = deriv(3,:)
!     eld      = deriv(3,:)
!     col(:,1) = eld
!     kcz      = kcz + MATMUL(col,row)*det*weights(i)
!   END DO gauss_pts_1        
!
!   kc = kcx*kx + kcy*ky + kcz*kz 
!    
!   storkc_pp here?
!
!   timest(10) = elap_time()
!
!   IF(numpe==1) PRINT *, " *** Computed element stiffness matrices in: ",    &
!                           timest(10)-timest(9), " s"
!
!------------------------------------------------------------------------------
! 9. Build the diagonal preconditioner
!------------------------------------------------------------------------------

! ALLOCATE(diag_precon_tmp(ntot,nels_pp))
! diag_precon_tmp = zero
!
! elements_4: DO iel = 1,nels_pp
!   DO i = 1,ndof
!     diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel) + storkc_pp(i,i,iel)
!   END DO
! END DO elements_4
!
! CALL scatter(diag_precon_pp,diag_precon_tmp)
!
! DEALLOCATE(diag_precon_tmp)
!
! timest(11) = elap_time()
!
! IF(numpe==1) PRINT *, " *** Built diagonal preconditioner in: ",            &
!                         timest(11)-timest(10), " s"
!  
!------------------------------------------------------------------------------
! 10. Read in fixed nodal displacements and assign to equations
!------------------------------------------------------------------------------
!
! IF(fixed_freedoms > 0) THEN
!
!   ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),valf(fixed_freedoms),    &
!            no_pp_temp(fixed_freedoms),sense(fixed_freedoms),                &
!            no_global(fixed_freedoms))
!
!   node = 0 ; no = 0 ; no_pp_temp = 0 ; sense = 0 ; no_global = 0
!   valf = zero
!
!   CALL read_fixed(job_name,numpe,node,sense,valf)
!   find_no2(g_g_pp,g_num_pp,node,sense,fixed_freedoms_pp,                    &
!            fixed_freedoms_start,no)
!   CALL MPI_ALLREDUCE(no,no_global,fixed_freedoms,MPI_INTEGER,MPI_MAX,       &
!                      MPI_COMM_WORLD,ier)
!   CALL reindex_fixed_nodes(ieq_start,no_global,no_pp_temp,                  &
!                            fixed_freedoms_pp,fixed_freedoms_start,neq_pp)
!
!   ALLOCATE(no_pp(fixed_freedoms_pp),store_pp(fixed_freedoms_pp))
!
!   no_pp    = 0
!   store_pp = zero
!   no_pp    = no_pp_temp(1:fixed_freedoms_pp)
!
!   DEALLOCATE(node,no,sense,no_pp_temp)
!
! END IF
!
! IF(fixed_freedoms == 0) fixed_freedoms_pp = 0
!
! DEALLOCATE(rest)
!
!------------------------------------------------------------------------------
! 11. Read in loaded nodes and get starting r_pp
!------------------------------------------------------------------------------
!
! IF(loaded_nodes > 0) THEN
!
!   ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes))
!
!   val = zero ; node = 0
!
!   CALL read_loads(job_name,numpe,node,val)
!   CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))
!
!   tload = SUM_P(r_pp(1:))
!
!   DEALLOCATE(node,val)
!
! ELSE
!
!   tload = zero
!
! END IF
!
! DEALLOCATE(g_g_pp)
!
! timest(12) = elap_time()
!
! IF(numpe==1) PRINT *, " *** Applied loads in:",                             &
!                         timest(12)-timest(11), " s"
!
!------------------------------------------------------------------------------
! 12. Invert the preconditioner. 
!     If there are fixed freedoms, first apply a penalty
!------------------------------------------------------------------------------
!
! IF(fixed_freedoms_pp > 0) THEN
!   DO i = 1,fixed_freedoms_pp
!      j =  no_pp(i) - ieq_start + 1
!      diag_precon_pp(j) = diag_precon_pp(j) + penalty
!      store_pp(i)       = diag_precon_pp(j)
!   END DO
! END IF
!
! diag_precon_pp = 1._iwp/diag_precon_pp
!
!-------------------- get starting r-----------------------------------------
    IF(loaded_freedoms>0) THEN
     CALL reindex_fixed_nodes(ieq_start,no,no_local_temp,num_no,              &
                              no_index_start,neq_pp)
     ALLOCATE(no_local(1:num_no)) ; no_local = no_local_temp(1:num_no)
     DEALLOCATE(no_local_temp)
        DO i = 1 , num_no
         r_pp(no_local(i)-ieq_start+1) = val(no_index_start + i - 1)
        END DO
    END IF            ;     q = SUM_P(r_pp)
    IF(numpe==it)THEN  ;  WRITE(11,'(A,E12.4)') "The total load is ", q
    END IF
    IF(fixed_freedoms>0) THEN
     CALL reindex_fixed_nodes(ieq_start,no_f,no_local_temp_f,        &
                              num_no,no_index_start,neq_pp)
     ALLOCATE(no_local(1:num_no)) ; no_local = no_local_temp_f(1:num_no)
     DEALLOCATE(no_local_temp_f)
        DO i = 1 , num_no        ; j=no_local(i) - ieq_start + 1
         diag_precon_pp(j)=diag_precon_pp(j) + penalty 
         r_pp(j) = diag_precon_pp(j) *  val_f(no_index_start + i - 1)
         store_pp(j) =  diag_precon_pp(j)
        END DO
    END IF
    diag_precon_pp=1._iwp/diag_precon_pp;d_pp=diag_precon_pp*r_pp; p_pp=d_pp
!--------------------preconditioned c. g. iterations---------------------------
       iters = 0
     iterations  :      DO 
             iters = iters + 1     ;    u_pp = 0._iwp  ; pmul_pp = .0_iwp
       CALL gather(p_pp,pmul_pp) 
       elements_2 : DO iel = 1, nels_pp
                    utemp_pp(:,iel) = MATMUL(kc,pmul_pp(:,iel)) 
       END DO elements_2  ;       CALL scatter(u_pp,utemp_pp)
    IF(fixed_freedoms>0) THEN
        DO i = 1 , num_no;  j = no_local(i)-ieq_start+1
         u_pp(j)=p_pp(j) * store_pp(j)
        END DO
    END IF
!--------------------------pcg equation solution-------------------------------
           up=DOT_PRODUCT_P(r_pp,d_pp); alpha= up/ DOT_PRODUCT_P(p_pp,u_pp)
           xnew_pp = x_pp + p_pp* alpha ; r_pp=r_pp - u_pp*alpha
           d_pp = diag_precon_pp*r_pp ;  beta=DOT_PRODUCT_P(r_pp,d_pp)/up
           p_pp=d_pp+p_pp*beta    
           CALL checon_par(xnew_pp,tol,converged,x_pp)    
           IF(converged .OR. iters==limit) EXIT
     END DO iterations
     IF(numpe==it)THEN
       WRITE(11,'(A,I5)')"The number of iterations to convergence was  ",iters 
       WRITE(11,'(A)')   "The  potentials are   :"
       WRITE(11,'(A)') "   Freedom       Potential"
       DO i = 1 , 4
        WRITE(11,'(I5,A,E12.4)') nres+i-1, "     ", xnew_pp(is+i-1)
       END DO
     END IF
  IF(numpe==it) WRITE(11,*) "This analysis took   ", elap_time( ) - timest(1)
  CALL MPI_FINALIZE(ier)
 END PROGRAM p123
