PROGRAM p123         
!------------------------------------------------------------------------------
!      program p12.3 three dimensional analysis of steady state heat equation
!      using 8-node brick elements, preconditioned conjugate gradient solver
!      diagonal preconditioner ; parallel version ; externally generated model 
!------------------------------------------------------------------------------
 USE precision; USE global_variables; USE mp_interface; USE input; USE output
 USE loading; USE timing; USE maths; USE gather_scatter; USE new_libary
 IMPLICIT NONE
 !  neq , ntot  are now global variables - not declared
 INTEGER, PARAMETER::ndim=3,nodof=1
 INTEGER::nod,nn,nr,nip,i,j,k,iters,limit,iel,num_no,no_index_start,l,   &
   idx1,node_end,node_start,nodes_pp,loaded_freedoms,fixed_freedoms,     &
   loaded_nodes,fixed_freedoms_pp,fixed_freedoms_start,                  &
   loaded_freedoms_pp,loaded_freedoms_start,nels,ndof,ielpe,npes_pp,     &
   argc,iargc,meshgen,partitioner
 REAL(iwp),PARAMETER::zero=0.0_iwp,penalty=1.e20_iwp
 REAL(iwp)::kx,ky,kz,det,tol,up,alpha,beta,q
  CHARACTER(LEN=15)     :: element
  CHARACTER(LEN=50)     :: fname,job_name,label
  CHARACTER(LEN=50)     :: program_name='xx11'
  LOGICAL               :: converged 

!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: points(:,:),kc(:,:),coord(:,:), weights(:)
  REAL(iwp),ALLOCATABLE :: p_g_co_pp(:,:,:), jac(:,:), der(:,:), deriv(:,:)
  REAL(iwp),ALLOCATABLE :: col(:,:),row(:,:),kcx(:,:),kcy(:,:),kcz(:,:)
  REAL(iwp),ALLOCATABLE :: diag_precon_pp(:),p_pp(:),r_pp(:),x_pp(:)
  REAL(iwp),ALLOCATABLE :: xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:)
  REAL(iwp),ALLOCATABLE :: d_pp(:),diag_precon_tmp(:,:),val(:,:),val_f(:)
  REAL(iwp),ALLOCATABLE :: store_pp(:),storkc_pp(:,:,:),eld(:),timest(:)
  REAL(iwp),ALLOCATABLE :: g_coord_pp(:,:,:),disp_pp(:),eld_pp(:,:)
  REAL(iwp),ALLOCATABLE :: flux_pp(:,:),storeflux_pp(:),kay(:,:),fun(:)
  REAL(iwp),ALLOCATABLE :: shape_integral_pp(:,:),flux_integral_pp(:,:),fluxnodes_pp(:)
  REAL(iwp),ALLOCATABLE :: shape_integral2_pp(:,:),flux_integral2_pp(:,:)
  INTEGER, ALLOCATABLE  :: rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:),no(:)
  INTEGER, ALLOCATABLE  :: no_f(:),no_local_temp(:),no_local_temp_f(:)
  INTEGER, ALLOCATABLE  :: no_local(:),no_pp(:),no_f_pp(:),no_pp_temp(:),no_global(:)
  INTEGER, ALLOCATABLE  :: sense(:),node(:)
 
!--------------------------input and initialisation-----------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)
 CALL read_p123(argv,numpe,element,fixed_freedoms,kx,ky,kz,limit,        &
   loaded_nodes,meshgen,nels,nip,nn,nod,nr,partitioner,tol)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof 
 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp))           &
 g_num_pp=0; g_coord_pp=zero
 CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
 IF (nr>0) THEN; ALLOCATE(rest(nr,nodof+1)); rest=0
 CALL read_rest(job_name,numpe,rest); END IF
 ALLOCATE (points(nip,ndim),coord(nod,ndim),jac(ndim,ndim),              &
   deriv(ndim,nod),kcx(ntot,ntot),weights(nip),g(ntot),der(ndim,nod),    &
   pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),col(ntot,1),num(nod),    &
   g_g_pp(ntot,nels_pp),kcy(ntot,ntot),no_local_temp(1),row(1,ntot),     &
   eld(ntot),no_f(1),no_local_temp_f(1),kcz(ntot,ntot),                  &
   storkc_pp(ntot,ntot,nels_pp))
!----------  find the steering array and equations per process -----------
 IF(nr>0) CALL rearrange_2(rest);  g_g_pp=0; neq=0
 IF(nr>0) THEN
   elements_1: DO iel = 1, nels_pp
     CALL find_g4(g_num_pp(:,iel),g_g_pp(:,iel),rest)
   END DO elements_1
 ELSE
   g_g_pp=g_num_pp  !When nr = 0, g_num_pp and g_g_pp are identical
 END IF
 neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
 IF(numpe==1)THEN
   OPEN(11,FILE=argv(1:nlen)//'.res',STATUS='REPLACE',ACTION='WRITE')              
   WRITE(11,'(A,I5,A)')"This job ran on ", npes,"  processes"
   WRITE(11,'(A,3(I7,A))')"There are ",nn," nodes",nr,                   &
     " restrained and   ",neq," equations"
   WRITE(11,*)"Time after setup  is  : ",elap_time()-timest(1)
 END IF 
 ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),        &
   u_pp(neq_pp),diag_precon_pp(neq_pp),d_pp(neq_pp))
 r_pp=zero; p_pp=zero; x_pp=zero; xnew_pp=zero; diag_precon_pp=zero
!-------------- element stiffness integration and storage ----------------
 CALL sample(element,points,weights); storkc_pp=zero
 elements_3: DO iel=1,nels_pp
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
 END DO elements_3
!------------------ build the diagonal preconditioner --------------------
 ALLOCATE(diag_precon_tmp(ntot,nels_pp)); diag_precon_tmp=zero
 elements_4: DO iel=1,nels_pp
   DO i=1,ndof
     diag_precon_tmp(i,iel)=diag_precon_tmp(i,iel)+storkc_pp(i,i,iel)
   END DO
 END DO elements_4; CALL scatter(diag_precon_pp,diag_precon_tmp)
 DEALLOCATE(diag_precon_tmp)
!------------- read in fixed freedoms and assign to equations ------------
 IF(fixed_freedoms>0) THEN
   ALLOCATE(node(fixed_freedoms),no(fixed_freedoms),
     val_f(fixed_freedoms),no_pp_temp(fixed_freedoms),
     sense(fixed_freedoms),no_global(fixed_freedoms))
   node=0; no=0; no_pp_temp=0; sense=0; no_global=0; val_f=zero
   CALL read_fixed(job_name,numpe,node,sense,val_f)
   CALL find_no2(g_g_pp,g_num_pp,node,sense,fixed_freedoms_pp,           &
     fixed_freedoms_start,no)
   CALL MPI_ALLREDUCE(no,no_global,fixed_freedoms,MPI_INTEGER,MPI_MAX,   &
     MPI_COMM_WORLD,ier)
   CALL reindex(ieq_start,no_global,no_pp_temp,fixed_freedoms_pp,        &
     fixed_freedoms_start,neq_pp)
   ALLOCATE(no_f_pp(fixed_freedoms_pp),store_pp(fixed_freedoms_pp))
   no_f_pp=0; store_pp=zero
   no_f_pp=no_pp_temp(1:fixed_freedoms_pp)
   DEALLOCATE(node,no,sense,no_pp_temp)
 END IF
 IF(fixed_freedoms == 0) fixed_freedoms_pp = 0
 IF (nr>0) DEALLOCATE(rest)
!--------------- read in loaded nodes and get starting r_pp --------------
 loaded_freedoms=loaded_nodes ! hack
 IF(loaded_freedoms>0) THEN
   ALLOCATE(node(loaded_freedoms),val(nodof,loaded_freedoms),            &
     no_pp_temp(loaded_freedoms)); val=zero; node=0
   CALL read_loads(job_name,numpe,node,val)
   CALL reindex(ieq_start,node,no_pp_temp,loaded_freedoms_pp,            &
     loaded_freedoms_start,neq_pp); ALLOCATE(no_pp(loaded_freedoms_pp))
   no_pp=no_pp_temp(1:loaded_freedoms_pp); DEALLOCATE(no_pp_temp)
   DO i = 1, loaded_freedoms_pp
     r_pp(no_pp(i)-ieq_start+1) = val(loaded_freedoms_start+i-1)
   END DO; q=SUM_P(r_pp); DEALLOCATE(node,val)
 END IF
!-------------------------- invert preconditioner ------------------------
 IF(fixed_freedoms_pp>0) THEN
   DO i=1,fixed_freedoms_pp
     j=no_f_pp(i)-ieq_start+1
     diag_precon_pp(j)=diag_precon_pp(j)+penalty
     store_pp(i)=diag_precon_pp(j)
   END DO
 END IF
 diag_precon_pp = 1._iwp/diag_precon_pp
!---------------- initialise preconditioned conjugate gradient -----------
 IF(fixed_freedoms_pp>0) THEN
   DO i=1,fixed_freedoms_pp
     j=no_f_pp(i)-ieq_start+1; k=fixed_freedoms_start+i-1
     r_pp(j)=store_pp(i)*val_f(k)
   END DO
 END IF
 d_pp=diag_precon_pp*r_pp; p_pp=d_pp; x_pp=zero
!--------------- preconditioned conjugate gradient iterations ------------
 iters=0
 iterations: DO 
   iters=iters+1; u_pp=zero; pmul_pp=zero; utemp_pp=zero 
   CALL gather(p_pp,pmul_pp) 
   elements_5 : DO iel = 1, nels_pp
     utemp_pp(:,iel) = MATMUL(storkc_pp(:,:,iel),pmul_pp(:,iel)) 
   END DO elements_5  
   CALL scatter(u_pp,utemp_pp)
   IF(fixed_freedoms_pp>0) THEN
     DO i=1,fixed_freedoms_pp
       j=no_f_pp(i)-ieq_start+1
       u_pp(j)=p_pp(j)*store_pp(i)
     END DO
   END IF
   up=DOT_PRODUCT_P(r_pp,d_pp); alpha=up/DOT_PRODUCT_P(p_pp,u_pp)
   xnew_pp=x_pp+p_pp*alpha; r_pp=r_pp-u_pp*alpha
   d_pp=diag_precon_pp*r_pp; beta=DOT_PRODUCT_P(r_pp,d_pp)/up
   p_pp=d_pp+p_pp*beta
   CALL checon_par(xnew_pp,tol,converged,x_pp)    
   IF(converged .OR. iters==limit) EXIT
 END DO iterations
 IF(numpe==1)THEN
   WRITE(11,'(A,I5)')"The number of iterations to convergence was  ",iters 
   WRITE(11,'(A,E12.4)')"The total load is                            ",q
   WRITE(11,'(A)')   "The  potentials are   :"
   WRITE(11,'(A)') "   Freedom       Potential"
   WRITE(11,'(A,E12.4)') "9901     ", xnew_pp(9901)
   WRITE(11,'(A,E12.4)') "9902     ", xnew_pp(9902)
   WRITE(11,'(A,E12.4)') "9903     ", xnew_pp(9903)
   WRITE(11,'(A,E12.4)') "9904     ", xnew_pp(9904)
 END IF
!----------------------- output nodal temperatures -----------------------
 IF(numpe==1)THEN; WRITE(ch,'(I6.6)') iy
   OPEN(12,file=argv(1:nlen)//".ensi.TTR-"//ch,status='replace',         &
     action='write')
   WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
   WRITE(12,'(A/A/A)') "part", "    1","coordinates"
 END IF
 disp_pp=zero; utemp_pp=zero; CALL gather(totd_pp(1:),utemp_pp)
 CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,          &
   node_start,node_end,utemp_pp,disp_pp,1)
 DO i=1,ndim; temp=zero
   DO j=1,nodes_pp; k=i+(ndim*(j-1)); temp(j)=disp_pp(k); END DO
   CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
 END DO; IF(numpe==1) CLOSE(12)
 IF(numpe==1)WRITE(11,*)"This analysis took : ",elap_time()-timest(1)
 CALL shutdown()
END PROGRAM p123
