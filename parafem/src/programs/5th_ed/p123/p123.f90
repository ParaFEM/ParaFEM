PROGRAM p123         
!-------------------------------------------------------------------------
!      program p12.3 three dimensional analysis of Laplace's equation
!      using 8-node bricks, preconditioned conjugate gradient solver
!      diagonal preconditioner; parallel; externally generated model 
!-------------------------------------------------------------------------
!USE mpi_wrapper  !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE new_library; IMPLICIT NONE
!neq,ntot  are now global variables - not declared
 INTEGER, PARAMETER::ndim=3,nodof=1
 INTEGER::nod,nn,nr,nip,i,j,k,iters,limit,iel,partitioner,meshgen,       &
   node_end,node_start,nodes_pp,loaded_freedoms,fixed_freedoms,          &
   fixed_freedoms_pp,fixed_freedoms_start,nlen,nres,is,it,               &
   loaded_freedoms_pp,loaded_freedoms_start,nels,ndof,npes_pp   
 REAL(iwp),PARAMETER::zero=0.0_iwp,penalty=1.e20_iwp
 REAL(iwp)::kx,ky,kz,det,tol,up,alpha,beta,q; CHARACTER(LEN=6)::ch
 CHARACTER(LEN=15)::element; CHARACTER(LEN=50)::argv 
 LOGICAL::converged=.false.
 REAL(iwp),ALLOCATABLE::points(:,:),weights(:),eld_pp(:,:),kay(:,:),     &
   fun(:),jac(:,:),der(:,:),deriv(:,:),col(:,:),row(:,:),                &
   kcx(:,:),kcy(:,:),kcz(:,:),diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:), &
   xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),val(:,:),       &
   diag_precon_tmp(:,:),store_pp(:),storkc_pp(:,:,:),eld(:),timest(:),   &
   val_f(:),g_coord_pp(:,:,:),ptl_pp(:)
 INTEGER,ALLOCATABLE::rest(:,:),g_num_pp(:,:),g_g_pp(:,:),no(:),         &
   no_pp(:),no_f_pp(:),no_pp_temp(:),sense(:),node(:)
!--------------------------input and initialisation-----------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)
 CALL read_p123(argv,numpe,element,fixed_freedoms,kx,ky,kz,limit,        &
   loaded_freedoms,meshgen,nels,nip,nn,nod,nr,nres,partitioner,tol)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof 
 ALLOCATE(g_num_pp(nod,nels_pp),g_coord_pp(nod,ndim,nels_pp)) 
 g_num_pp=0; g_coord_pp=zero
 CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 IF (nr>0) THEN; ALLOCATE(rest(nr,nodof+1)); rest=0
 CALL read_rest(argv,numpe,rest); END IF
 ALLOCATE (points(nip,ndim),jac(ndim,ndim),storkc_pp(ntot,ntot,nels_pp), &
   deriv(ndim,nod),kcx(ntot,ntot),weights(nip),der(ndim,nod),            &
   pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),col(ntot,1),eld(ntot),   &
   g_g_pp(ntot,nels_pp),kcy(ntot,ntot),row(1,ntot),kcz(ntot,ntot))
!----------  find the steering array and equations per process -----------
 timest(2)=elap_time(); g_g_pp=0; neq=0
 IF(nr>0) THEN; CALL rearrange_2(rest)
   elements_0: DO iel=1, nels_pp
     CALL find_g4(g_num_pp(:,iel),g_g_pp(:,iel),rest)
   END DO elements_0
 ELSE
   g_g_pp=g_num_pp  !When nr = 0, g_num_pp and g_g_pp are identical
 END IF
 neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
 DO i=1,neq_pp;IF(nres==ieq_start+i-1)THEN;it=numpe;is=i;END IF;END DO
 IF(numpe==it)THEN
   OPEN(11,FILE=argv(1:nlen)//'.res',STATUS='REPLACE',ACTION='WRITE') 
   WRITE(11,'(A,I5,A)')"This job ran on ", npes,"  processes"
   WRITE(11,'(A,3(I7,A))')"There are ",nn," nodes",nr,                   &
     " restrained and   ",neq," equations"
   WRITE(11,'(A,F10.4)')"Time after setup is ",elap_time()-timest(1)
 END IF
 ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),        &
   u_pp(neq_pp),diag_precon_pp(neq_pp),d_pp(neq_pp))
 r_pp=zero; p_pp=zero; x_pp=zero; xnew_pp=zero; diag_precon_pp=zero
!-------------- element stiffness integration and storage ----------------
 CALL sample(element,points,weights); storkc_pp=zero
 elements_1: DO iel=1,nels_pp
   kcx=zero; kcy=zero; kcz=zero
   gauss_pts_1:  DO i=1,nip
     CALL shape_der (der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
     det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
     row(1,:)=deriv(1,:); eld=deriv(1,:); col(:,1)=eld
     kcx=kcx+MATMUL(col,row)*det*weights(i); row(1,:)=deriv(2,:)
     eld=deriv(2,:); col(:,1)=eld
     kcy=kcy+MATMUL(col,row)*det*weights(i); row(1,:)=deriv(3,:)
     eld=deriv(3,:); col(:,1)=eld
     kcz=kcz+MATMUL(col,row)*det*weights(i)
   END DO gauss_pts_1        
   storkc_pp(:,:,iel)=kcx*kx+kcy*ky+kcz*kz 
 END DO elements_1
!------------------ build the diagonal preconditioner --------------------
 ALLOCATE(diag_precon_tmp(ntot,nels_pp)); diag_precon_tmp=zero
 elements_1a: DO iel=1,nels_pp
   DO i=1,ndof
     diag_precon_tmp(i,iel)=diag_precon_tmp(i,iel)+storkc_pp(i,i,iel)
   END DO
 END DO elements_1a; CALL scatter(diag_precon_pp,diag_precon_tmp)
 DEALLOCATE(diag_precon_tmp)
!------------- read in fixed freedoms and assign to equations ------------
 IF(fixed_freedoms>0) THEN
   ALLOCATE(node(fixed_freedoms),no_pp_temp(fixed_freedoms),             &
     val_f(fixed_freedoms),no(fixed_freedoms),sense(fixed_freedoms))
   node=0; no=0; no_pp_temp=0; sense=0; val_f=zero
   CALL read_fixed(argv,numpe,node,sense,val_f)
   CALL find_no2(g_g_pp,g_num_pp,node,sense,no)
   CALL reindex(ieq_start,no,no_pp_temp,fixed_freedoms_pp,               &
     fixed_freedoms_start,neq_pp)
   ALLOCATE(no_f_pp(fixed_freedoms_pp),store_pp(fixed_freedoms_pp))
   no_f_pp=0; store_pp=zero; no_f_pp=no_pp_temp(1:fixed_freedoms_pp)
   DEALLOCATE(node,no,sense,no_pp_temp)
 END IF
 IF(fixed_freedoms==0) fixed_freedoms_pp=0; IF(nr>0) DEALLOCATE(rest)
!------------ read in loaded freedoms and get starting r_pp --------------
 IF(loaded_freedoms>0) THEN
   ALLOCATE(node(loaded_freedoms),val(nodof,loaded_freedoms),            &
     no_pp_temp(loaded_freedoms)); val=zero; node=0
   CALL read_loads(argv,numpe,node,val)
   CALL reindex(ieq_start,node,no_pp_temp,loaded_freedoms_pp,            &
     loaded_freedoms_start,neq_pp); ALLOCATE(no_pp(loaded_freedoms_pp))
   no_pp=no_pp_temp(1:loaded_freedoms_pp); DEALLOCATE(no_pp_temp)
   DO i = 1, loaded_freedoms_pp
     r_pp(no_pp(i)-ieq_start+1) = val(1,loaded_freedoms_start+i-1)
   END DO; q=SUM_P(r_pp); DEALLOCATE(node,val)
 END IF
!-------------------------- invert preconditioner ------------------------
 IF(fixed_freedoms_pp>0) THEN
   DO i=1,fixed_freedoms_pp; j=no_f_pp(i)-ieq_start+1
     diag_precon_pp(j)=diag_precon_pp(j)+penalty
     store_pp(i)=diag_precon_pp(j)
   END DO
 END IF;  diag_precon_pp = 1._iwp/diag_precon_pp
!---------------- initialise preconditioned conjugate gradient -----------
 IF(fixed_freedoms_pp>0) THEN
   DO i=1,fixed_freedoms_pp
     j=no_f_pp(i)-ieq_start+1; k=fixed_freedoms_start+i-1
     r_pp(j)=store_pp(i)*val_f(k)
   END DO
 END IF;  d_pp=diag_precon_pp*r_pp; p_pp=d_pp; x_pp=zero
!--------------- preconditioned conjugate gradient iterations ------------
 iters=0; timest(3)=elap_time()
 iterations: DO 
   iters=iters+1; u_pp=zero; pmul_pp=zero; utemp_pp=zero 
   CALL gather(p_pp,pmul_pp) 
   elements_2 : DO iel = 1, nels_pp
     utemp_pp(:,iel) = MATMUL(storkc_pp(:,:,iel),pmul_pp(:,iel)) 
   END DO elements_2 ; CALL scatter(u_pp,utemp_pp)
   IF(fixed_freedoms_pp>0) THEN
     DO i=1,fixed_freedoms_pp
       j=no_f_pp(i)-ieq_start+1; u_pp(j)=p_pp(j)*store_pp(i)
     END DO
   END IF
   up=DOT_PRODUCT_P(r_pp,d_pp); alpha=up/DOT_PRODUCT_P(p_pp,u_pp)
   xnew_pp=x_pp+p_pp*alpha; r_pp=r_pp-u_pp*alpha
   d_pp=diag_precon_pp*r_pp; beta=DOT_PRODUCT_P(r_pp,d_pp)/up
   p_pp=d_pp+p_pp*beta; CALL checon_par(xnew_pp,tol,converged,x_pp)    
   IF(converged .OR. iters==limit) EXIT
 END DO iterations
 timest(4)=elap_time()
 IF(numpe==it)THEN
   WRITE(11,'(A,I5)')"The number of iterations to convergence was ",iters 
   WRITE(11,'(A,E12.4)')"The total load is ",q
   WRITE(11,'(A)')   "The potentials are:"
   WRITE(11,'(A)') " Freedom       Potential"
   DO i=1,4
     WRITE(11,'(I8,A,E12.4)') nres+i-1, "     ", xnew_pp(is+i-1)
   END DO
   WRITE(11,'(A,F10.4)')"Integration, preconditioning and loading took ",&
     timest(3)-timest(2)
   WRITE(11,'(A,F10.4)') "Time spent in the solver was ",                &
     timest(4)-timest(3)
 END IF
!------------------- output potentials for ParaView ----------------------
 IF(numpe==1)THEN; WRITE(ch,'(I6.6)') numpe
   OPEN(12,file=argv(1:nlen)//".ensi.NDPTL-"//ch,status='replace',       &
     action='write')
   WRITE(12,'(A)') "Alya Ensight Gold --- Scalar per-node variable file"
   WRITE(12,'(A/A/A)') "part", "    1","coordinates"
 END IF
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 ALLOCATE(ptl_pp(nodes_pp*ndim))
 ptl_pp=zero; utemp_pp=zero; CALL gather(xnew_pp(1:),utemp_pp)
 CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,         &
   node_start,node_end,utemp_pp,ptl_pp,1)
 CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,ptl_pp)
 IF(numpe==it)                                                           &
   WRITE(11,'(A,F10.4)') "This analysis took ",elap_time()-timest(1)
 IF(numpe==it) CLOSE(11); IF(numpe==1) CLOSE(12); CALL SHUTDOWN()
END PROGRAM p123
