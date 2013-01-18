PROGRAM p129       
!-------------------------------------------------------------------------
!      Program 12.9 forced vibration of a 3d elastic solid
!      Lumped or consistent mass
!      Implicit integration by theta method : parallel version
!-------------------------------------------------------------------------
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE partition; USE elements; USE steering; USE pcg; IMPLICIT NONE
! neq,ntot are now global variables - not declared 
 INTEGER,PARAMETER::nodof=3,nst=6
 INTEGER::nn,nr,nip,nod,i,j,k,iel,ndim=3,nstep,npri,iters,limit,ndof,    &
   nels,npes_pp,node_end,node_start,nodes_pp,loaded_nodes,               &
   meshgen,partitioner
 REAL(iwp),PARAMETER::zero=0.0_iwp    
 REAL(iwp)::e,v,det,rho,alpha1,beta1,omega,theta,period,pi,dtim,volume,  &
   c1,c2,c3,c4,real_time,tol,big,up,alpha,beta,tload     
 CHARACTER(LEN=15)::element,io_type; CHARACTER(LEN=50)::argv
 CHARACTER(LEN=6)::step; LOGICAL::consistent=.TRUE.,converged
!------------------------- dynamic arrays --------------------------------
 REAL(iwp),ALLOCATABLE::loads_pp(:),points(:,:),dee(:,:),fun(:),jac(:,:),&
   der(:,:),deriv(:,:),weights(:),bee(:,:),km(:,:),g_coord_pp(:,:,:),    &
   x1_pp(:),d1x1_pp(:),d2x1_pp(:),emm(:,:),ecm(:,:),x0_pp(:),d1x0_pp(:), &
   d2x0_pp(:),store_km_pp(:,:,:),vu_pp(:),store_mm_pp(:,:,:),u_pp(:),    &
   p_pp(:),d_pp(:),x_pp(:),xnew_pp(:),pmul_pp(:,:),utemp_pp(:,:),        &
   diag_precon_pp(:),diag_precon_tmp(:,:),temp_pp(:,:,:),disp_pp(:),     &
   loadNodeValue(:,:),fext_pp(:),val(:,:),timest(:)
 INTEGER,ALLOCATABLE::rest(:,:),g(:),g_num_pp(:,:),g_g_pp(:,:),          &
   loadNodeNum(:),node(:)       
!----------------------- input and initialisation ------------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)
 CALL read_p129(argv,numpe,alpha1,beta1,e,element,limit,loaded_nodes,&
  meshgen,nels,nip,nn,nod,npri,nr,nstep,omega,partitioner,rho,theta,tol,v)
 CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof
 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),           &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nels,nn,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest); timest(2)=elap_time()
 ALLOCATE(points(nip,ndim),g(ntot),fun(nod),dee(nst,nst),jac(ndim,ndim), &
   weights(nip),der(ndim,nod),deriv(ndim,nod),bee(nst,ntot),             &
   store_km_pp(ntot,ntot,nels_pp),utemp_pp(ntot,nels_pp),ecm(ntot,ntot), &
   pmul_pp(ntot,nels_pp),store_mm_pp(ntot,ntot,nels_pp),emm(ntot,ntot),  &
   temp_pp(ntot,ntot,nels_pp),diag_precon_tmp(ntot,nels_pp),             &
   g_g_pp(ntot,nels_pp))
!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_1: DO iel = 1, nels_pp
   CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_1
 elements_2: DO iel=1,nels_pp  
   i=MAXVAL(g_g_pp(:,iel)); IF(i>neq) neq=i
 END DO elements_2  
 neq=MAX_INTEGER_P(neq); CALL calc_neq_pp; CALL calc_npes_pp(npes,npes_pp)
 CALL make_ggl(npes_pp,npes,g_g_pp)
 DO i=1,neq_pp; IF(nres==ieq_start+i-1)THEN;it=numpe;is=i;END IF;END DO
 IF(numpe==it) THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A,I6,A)') "This job ran on ",npes," processes"
   WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes ",nr,                &
     " restrained and ", neq," equations"
   WRITE(11,'(A,F10.4)') "Time to read input was:",timest(2)-timest(1)
   WRITE(11,'(A,F10.4)') "Time after setup was:",elap_time()-timest(1)
 END IF 
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 ALLOCATE(disp_pp(nodes_pp*ndim),temp(nodes_pp),eld_pp(ntot,nels_pp))
 disp_pp=zero; eld_pp=zero, temp=0
 ALLOCATE(x0_pp(neq_pp),d1x0_pp(neq_pp),x1_pp(neq_pp),vu_pp(neq_pp),     &
   diag_precon_pp(neq_pp),u_pp(neq_pp),d2x0_pp(neq_pp),loads_pp(neq_pp), &
   d1x1_pp(neq_pp),d2x1_pp(neq_pp),d_pp(neq_pp),p_pp(neq_pp),            &
   x_pp(neq_pp),xnew_pp(neq_pp),fext_pp(neq_pp))
 x0_pp=zero; d1x0_pp=zero; x1_pp=zero; vu_pp=zero; diag_precon_pp=zero
 u_pp=zero; d2x0_pp=zero; loads_pp=zero; d1x1_pp=zero; d2x1_pp=zero
 d_pp=zero; p_pp=zero; x_pp=zero; xnew_pp=zero; fext_pp=zero 
!--- element stiffness and mass integration, storage and preconditioner --
 CALL deemat(e,v,dee); CALL sample(element,points,weights)
 store_km_pp=zero; store_mm_pp=zero; diag_precon_tmp=zero
 pi=ACOS(-1._iwp); period=2._iwp*pi/omega; dtim=period/20._iwp
 c1=(1._iwp-theta)*dtim; c2=beta1-c1; c3=alpha1+1._iwp/(theta * dtim)
 c4=beta1+theta*dtim
 elements_3: DO iel=1,nels_pp
   volume=zero; emm=zero; ecm=zero
   gauss_points_1: DO i=1,nip     
     CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
     det=determinant(jac); CALL invert(jac); deriv=matmul(jac,der)
     CALL beemat(deriv,bee)
     store_km_pp(:,:,iel)=store_km_pp(:,:,iel) +                         & 
                          MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) *       &
                          det*weights(i)
     volume=volume+det*weights(i); CALL shape_fun(fun,points,i)
     IF(consistent)THEN; CALL ecmat(ecm,fun,ntot,nodof)
       ecm=ecm*det*weights(i)*rho; emm=emm+ecm
     END IF
   END DO gauss_points_1   
   IF(.NOT.consistent)THEN
     DO i=1,ntot; emm(i,i)=volume*rho/13._iwp; END DO
     DO i=1,19,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=2,20,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=3,21,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=37,55,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=38,56,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=39,57,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
   END IF
   store_mm_pp(:,:,iel)=emm
 END DO elements_3
 diag_precon_tmp=zero
 elements_4: DO iel=1,nels_pp
   DO k=1,ntot
     diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)    +                  &
                           store_mm_pp(k,k,iel)*c3+store_km_pp(k,k,iel)*c4
   END DO
 END DO elements_4; CALL scatter(diag_precon_pp,diag_precon_tmp)
 diag_precon_pp=1._iwp/diag_precon_pp; DEALLOCATE(diag_precon_tmp)
!---------------------------------- loads --------------------------------
 IF(loaded_nodes>0) THEN
   ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes)); val=zero; node=0
   CALL read_loads(job_name,numpe,node,val)
   CALL load(g_g_pp,g_num_pp,node,val,fext_pp(1:))
   tload=SUM_P(fext_pp(1:)); DEALLOCATE(node,val)
 END IF
!---------------------------- initial conditions -------------------------
 x0_pp=zero; d1x0_pp=zero; d2x0_pp=zero; real_time=zero
!---------------------------- time stepping loop ------------------------- 
 timesteps: DO j=1,nstep
   real_time=real_time+dtim; loads_pp=zero; u_pp=zero; vu_pp=zero
   elements_5: DO iel=1,nels_pp    ! gather for rhs multiply
     temp_pp(:,:,iel)=store_km_pp(:,:,iel)*c2+store_mm_pp(:,:,iel)*c3
   END DO elements_5; CALL gather(x0_pp,pmul_pp)
   DO iel=1,nels_pp
     utemp_pp(:,iel)=MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
   END DO; CALL scatter(u_pp,utemp_pp)
!---------------------------- velocity part ------------------------------
   temp_pp=store_mm_pp/theta; CALL gather(d1x0_pp,pmul_pp)
   DO iel=1,nels_pp
     utemp_pp(:,iel)=MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
   END DO; CALL scatter(vu_pp,utemp_pp)   ! doesn't add to last u_pp
   loads_pp=fext_pp*(theta*dtim*cos(omega*real_time)+c1*                 &
                           cos(omega*(real_time-dtim))) 
   loads_pp=u_pp+vu_pp+loads_pp
!----------- solve simultaneous equations by PCG -------------------------
   d_pp=diag_precon_pp*loads_pp; p_pp=d_pp; x_pp=zero; iters=0
   iterations: DO 
     iters=iters+1; u_pp=zero; vu_pp=zero
     temp_pp=store_mm_pp*c3+store_km_pp*c4
     CALL gather(p_pp,pmul_pp)
     elements_6: DO iel=1,nels_pp
       utemp_pp(:,iel)=MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
     END DO elements_6; CALL scatter(u_pp,utemp_pp)
     up=DOT_PRODUCT_P(loads_pp,d_pp); alpha=up/DOT_PRODUCT_P(p_pp,u_pp)
     xnew_pp=x_pp+p_pp*alpha; loads_pp=loads_pp-u_pp*alpha
     d_pp=diag_precon_pp*loads_pp; beta=DOT_PRODUCT_P(loads_pp,d_pp)/up
     p_pp=d_pp+p_pp*beta; u_pp=xnew_pp
     CALL checon_par(xnew_pp,tol,converged,x_pp)
     IF(converged.OR.iters==limit)EXIT
   END DO iterations
   x1_pp=xnew_pp 
   d1x1_pp=(x1_pp-x0_pp)/(theta*dtim)-d1x0_pp*(1._iwp-theta)/theta
   d2x1_pp=(d1x1_pp-d1x0_pp)/(theta*dtim)-d2x0_pp*(1._iwp-theta)/theta
   x0_pp=x1_pp; d1x0_pp=d1x1_pp; d2x0_pp=d2x1_pp; utemp_pp=zero
   IF(numpe==1) THEN;  WRITE(ch,'(I6.6)') j
     OPEN(12,file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',     &
          action='write')
     WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
     WRITE(12,'(A/A/A)') "part", "     1","coordinates"
   END IF
   CALL gather(x1_pp(1:),eld_pp); disp_pp=zero
   CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,        &
                       node_start,node_end,eld_pp,disp_pp,1)
   DO i=1,ndim ; temp=zero
     DO j=1,nodes_pp; k=i+(ndim*(j-1)); temp(j)=disp_pp(k); END DO
       CALL dismsh_ensi_p(12,j,nodes_pp,npes,numpe,1,temp)
   END DO ; IF(numpe==1) CLOSE(12)
 END DO timesteps
 IF(numpe==it) THEN
   WRITE(11,'(A,F10.4)')"This analysis took  :",elap_time()-timest(1)
 END IF; CALL shutdown() 
END PROGRAM p129
