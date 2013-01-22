PROGRAM p126     
!-------------------------------------------------------------------------
!      Program 12.6 steady state  3-d Navier-Stokes equation
!      using 20-node velocity hexahedral elements  
!      coupled to 8-node pressure hexahedral elements : u-p-v-w order
!      element by element solution using BiCGSTAB(L) : parallel version
!-------------------------------------------------------------------------
!USE mpi_wrapper  !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE input; 
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE partition; USE elements; USE steering; USE bicg; USE fluid
 USE pcg ! This is to access CHECON_PAR, perhaps should move out of module
 IMPLICIT NONE
! neq,ntot are now global variables - not declared
 INTEGER::nn,nip,nodof=4,nod=20,nodf=8,ndim=3,cj_tot,i,j,k,l,iel,ell,    &
   limit,fixed_equations,iters,cjiters,cjits,nr,n_t,num_no,              &
   no_index_start,nres,is,it,nlen,nels,ndof,npes_pp,argc,iargc,meshgen,  &
   partitioner,node_end,node_start,nodes_pp  
 REAL(iwp):: visc,rho,rho1,det,ubar,vbar,wbar,tol,cjtol,alpha,beta,      &
   penalty,x0,pp,kappa,gama,omega,norm_r,r0_norm,error
 REAL(iwp),PARAMETER::zero=0.0_iwp,one=1.0_iwp
 LOGICAL::converged,cj_converged
 CHARACTER(LEN=15)::element='hexahedron'
 CHARACTER(LEN=50)::argv
!--------------------------- dynamic arrays ------------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),coord(:,:),derivf(:,:),fun(:),       &
   jac(:,:),kay(:,:),der(:,:),deriv(:,:),weights(:),derf(:,:),funf(:),   &
   coordf(:,:),g_coord_pp(:,:,:),c11(:,:),c21(:,:),c12(:,:),val(:),      &
   wvel(:),ke(:,:),c23(:,:),c32(:,:),x_pp(:),b_pp(:),r_pp(:,:),temp(:),  &
   funny(:,:),row1(:,:),row2(:,:),uvel(:),vvel(:),funnyf(:,:),rowf(:,:), &
   storke_pp(:,:,:),diag_pp(:),utemp_pp(:,:),xold_pp(:),c24(:,:),        &
   c42(:,:),row3(:,:),u_pp(:,:),rt_pp(:),y_pp(:),y1_pp(:),s(:),Gamma(:), &
   GG(:,:),diag_tmp(:,:),store_pp(:),pmul_pp(:,:),timest(:),disp_pp(:)
 INTEGER,ALLOCATABLE::rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:),   &
   no(:),g_t(:),no_local(:),no_local_temp(:)
!---------------------- input and initialisation -------------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)
 CALL read_p126(argv,numpe,cjits,cjtol,ell,fixed_equations,kappa,limit,  &
   meshgen,nels,nip,nn,nr,nres,partitioner,penalty,rho,tol,x0,visc) 
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ntot=nod+nodf+nod+nod; n_t=nod*nodof
 ALLOCATE(g_num_pp(nod,nels_pp),g_coord_pp(nod,ndim,nels_pp),            &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nels,nn,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest);timest(2)=elap_time()
 ALLOCATE(points(nip,ndim),coord(nod,ndim),derivf(ndim,nodf),fun(nod),   &
   jac(ndim,ndim),kay(ndim,ndim),der(ndim,nod),deriv(ndim,nod),          &
   derf(ndim,nodf),funf(nodf),coordf(nodf,ndim),funny(nod,1),            &
   g_g_pp(ntot,nels_pp),c11(nod,nod),c12(nod,nodf),c21(nodf,nod),        &
   ke(ntot,ntot),c24(nodf,nod),c42(nod,nodf),num(nod),c32(nod,nodf),     &
   c23(nodf,nod),uvel(nod),vvel(nod),row1(1,nod),funnyf(nodf,1),         &
   rowf(1,nodf),no_local_temp(fixed_equations),wvel(nod),row3(1,nod),    &
   storke_pp(ntot,ntot,nels_pp),g_t(n_t),GG(ell+1,ell+1),g(ntot),        &
   Gamma(ell+1),no(fixed_equations),val(fixed_equations),                &
   diag_tmp(ntot,nels_pp),utemp_pp(ntot,nels_pp),row2(1,nod),s(ell+1),   &
   pmul_pp(ntot,nels_pp),weights(nip))
!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_0: DO iel=1,nels_pp
   CALL find_g3(g_num_pp(:,iel),g_t,rest)
   CALL g_t_g_ns(nod,g_t,g_g_pp(:,iel))
 END DO elements_0
 neq=MAXVAL(g_g_pp); neq=MAX_INTEGER_P(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp)
 CALL make_ggl(npes_pp,npes,g_g_pp); DEALLOCATE(g_g_pp)
 DO i=1,neq_pp; IF(nres==ieq_start+i-1)THEN;it=numpe;is=i;END IF;END DO
 IF(numpe==it) THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A,I6,A)') "This job ran on ",npes," processes"
   WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes ",nr,                &
     " restrained and ", neq," equations"
   WRITE(11,'(A,F10.4)') "Time to read input was:",timest(2)-timest(1)
   WRITE(11,'(A,F10.4)') "Time after setup was:",elap_time()-timest(1)
 END IF
 ALLOCATE(x_pp(neq_pp),rt_pp(neq_pp),r_pp(neq_pp,ell+1),b_pp(neq_pp),    &
   u_pp(neq_pp,ell+1),diag_pp(neq_pp),xold_pp(neq_pp),y_pp(neq_pp),      &
   y1_pp(neq_pp),store_pp(neq_pp))
 x_pp=zero; rt_pp=zero; r_pp=zero; u_pp=zero; b_pp=zero; diag_pp=zero
 xold_pp=zero; y_pp=zero; y1_pp=zero; store_pp=zero
!-------------------------- organise fixed equations ---------------------
 CALL read_loads_ns(argv,numpe,no,val)
 CALL reindex_fixed_nodes(ieq_start,no,no_local_temp,num_no,             &
   no_index_start,neq_pp); ALLOCATE(no_local(1:num_no))
 no_local = no_local_temp(1:num_no); DEALLOCATE(no_local_temp)
!------------------------- main iteration loop ---------------------------
 CALL sample(element,points,weights); uvel=zero; vvel=zero; wvel=zero
 kay=zero; iters=0; cj_tot=0; kay(1,1)=visc/rho; kay(2,2)=visc/rho
 kay(3,3)=visc/rho; timest(3)=elap_time()
 iterations: DO
   iters=iters+1; ke=zero; storke_pp=zero; diag_pp=zero; utemp_pp=zero
   b_pp=zero; pmul_pp=zero; CALL gather(x_pp,utemp_pp)
   CALL gather(xold_pp,pmul_pp)
!-------------------- element stiffness integration ----------------------
   elements_1: DO iel=1,nels_pp
     uvel=(utemp_pp(1:nod,iel)+pmul_pp(1:nod,iel))*.5_iwp
     DO i=nod+nodf+1,nod+nodf+nod
       vvel(i-nod-nodf)=(utemp_pp(i,iel)+pmul_pp(i,iel))*.5_iwp
     END DO
     DO i=nod+nodf+nod+1,ntot
       wvel(i-nod-nodf-nod)=(utemp_pp(i,iel)+pmul_pp(i,iel))*.5_iwp
     END DO                                                        
     c11=zero; c12=zero; c21=zero; c23=zero; c32=zero; c24=zero; c42=zero
     gauss_points_1: DO i=1,nip
!------------------------ velocity contribution --------------------------
       CALL shape_fun(funny(:,1),points,i)
       ubar=DOT_PRODUCT(funny(:,1),uvel);vbar=DOT_PRODUCT(funny(:,1),vvel)
       wbar=DOT_PRODUCT(funny(:,1),wvel)
       IF(iters==1)THEN; ubar=one; vbar=zero; wbar=zero; END IF
       CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
       det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
       row1(1,:)=deriv(1,:); row2(1,:)=deriv(2,:); row3(1,:)=deriv(3,:)
       c11=c11+                                                          &
           MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i) +   &
           MATMUL(funny,row1)*det*weights(i)*ubar +                      &
           MATMUL(funny,row2)*det*weights(i)*vbar +                      &
           MATMUL(funny,row3)*det*weights(i)*wbar
!------------------------ pressure contribution --------------------------
       CALL shape_fun(funnyf(:,1),points,i); CALL shape_der(derf,points,i)
       coordf(1:4,:)=g_coord_pp(1:7:2,:,iel)
       coordf(5:8,:)=g_coord_pp(13:19:2,:,iel); jac=MATMUL(derf,coordf)
       det=determinant(jac); CALL invert(jac); derivf=MATMUL(jac,derf)
       rowf(1,:)=derivf(1,:)
       c12=c12+MATMUL(funny,rowf)*det*weights(i)/rho
       rowf(1,:)=derivf(2,:)
       c32=c32+MATMUL(funny,rowf)*det*weights(i)/rho
       rowf(1,:)=derivf(3,:)
       c42=c42+MATMUL(funny,rowf)*det*weights(i)/rho
       c21=c21+MATMUL(funnyf,row1)*det*weights(i)
       c23=c23+MATMUL(funnyf,row2)*det*weights(i) 
       c24=c24+MATMUL(funnyf,row3)*det*weights(i)
     END DO gauss_points_1
     CALL formupvw(storke_pp,iel,c11,c12,c21,c23,c32,c24,c42)
   END DO elements_1
!----------------------- build the preconditioner ------------------------
 diag_tmp=zero
 elements_2: DO iel=1,nels_pp; DO k=1,ntot
   diag_tmp(k,iel)=diag_tmp(k,iel)+storke_pp(k,k,iel); END DO
 END DO elements_2; CALL scatter(diag_pp,diag_tmp)
!------------------- prescribed values of velocity and pressure ----------
 DO i=1,num_no; k=no_local(i)-ieq_start+1
   diag_pp(k)=diag_pp(k)+penalty
   b_pp(k)=diag_pp(k)*val(no_index_start+i-1); store_pp(k)=diag_pp(k)
 END DO   
!---------- solve the equations element-by-element using BiCGSTAB --------
!-------------------------- initialisation phase -------------------------
 IF(iters==1) x_pp=x0; pmul_pp=zero; y1_pp=zero; y_pp=x_pp
   CALL gather(y_pp,pmul_pp)
   elements_3: DO iel=1,nels_pp  
     utemp_pp(:,iel)=MATMUL(storke_pp(:,:,iel),pmul_pp(:,iel))
   END DO elements_3; CALL scatter(y1_pp,utemp_pp)
   DO i=1,num_no; k=no_local(i)-ieq_start+1
     y1_pp(k)=y_pp(k)*store_pp(k)
   END DO; y_pp=y1_pp; rt_pp=b_pp-y_pp; r_pp=zero; r_pp(:,1)=rt_pp
   u_pp=zero; gama=one; omega=one; k=0; norm_r=norm_p(rt_pp)
   r0_norm=norm_r; error=one; cjiters=0
!----------------------- BiCGSTAB(ell) iterations ------------------------
   bicg_iterations: DO; cjiters=cjiters+1   
     cj_converged=error<cjtol; IF(cjiters==cjits.OR.cj_converged) EXIT
     gama=-omega*gama; y_pp=r_pp(:,1)
     DO j=1,ell
       rho1=DOT_PRODUCT_P(rt_pp,y_pp); beta=rho1/gama
       u_pp(:,1:j)=r_pp(:,1:j)-beta*u_pp(:,1:j)
       pmul_pp=zero; y_pp=u_pp(:,j); y1_pp=zero; CALL gather(y_pp,pmul_pp)
       elements_4: DO iel=1,nels_pp
         utemp_pp(:,iel)=MATMUL(storke_pp(:,:,iel),pmul_pp(:,iel))
       END DO elements_4; CALL scatter(y1_pp,utemp_pp)
       DO i=1,num_no; l=no_local(i)-ieq_start+1
         y1_pp(l)=y_pp(l)*store_pp(l)
       END DO; y_pp=y1_pp; u_pp(:,j+1)=y_pp
       gama=DOT_PRODUCT_P(rt_pp,y_pp); alpha=rho1/gama
       x_pp=x_pp+alpha*u_pp(:,1)
       r_pp(:,1:j)=r_pp(:,1:j)-alpha*u_pp(:,2:j+1)
       pmul_pp=zero; y_pp=r_pp(:,j); y1_pp=zero; CALL gather(y_pp,pmul_pp)
       elements_5: DO iel=1,nels_pp
         utemp_pp(:,iel)=MATMUL(storke_pp(:,:,iel),pmul_pp(:,iel))
       END DO elements_5; CALL scatter(y1_pp,utemp_pp)
       DO i=1,num_no; l=no_local(i)-ieq_start+1
         y1_pp(l)=y_pp(l)*store_pp(l)
       END DO; y_pp=y1_pp; r_pp(:,j+1)=y_pp
     END DO
     DO i=1,ell+1; DO j=1,ell+1
       GG(i,j)=DOT_PRODUCT_P(r_pp(:,i),r_pp(:,j))
     END DO; END DO
     CALL form_s(GG,ell,kappa,omega,Gamma,s) 
     x_pp=x_pp-MATMUL(r_pp,s); r_pp(:,1)=MATMUL(r_pp,Gamma)
     u_pp(:,1)=MATMUL(u_pp,Gamma); norm_r=norm_p(r_pp(:,1))
     error=norm_r/r0_norm; k=k+1
   END DO bicg_iterations            
   b_pp=x_pp-xold_pp; pp=norm_p(b_pp); cj_tot=cj_tot+cjiters
   IF(numpe==it) THEN
     WRITE(11,'(A,E12.4)') "Norm of the error is:", pp
     WRITE(11,'(A,I6,A)') "It took BiCGSTAB(L) ", cjiters,               &
       " iterations to converge"; END IF
   CALL checon_par(x_pp,tol,converged,xold_pp)
   IF(converged.OR.iters==limit)EXIT
 END DO iterations; timest(4)=elap_time()
 DEALLOCATE(rt_pp,r_pp,u_pp,b_pp,diag_pp,xold_pp,y_pp,y1_pp,store_pp)
 DEALLOCATE(storke_pp,pmul_pp)
!------------------------- output results --------------------------------
 IF(numpe==it) THEN
   WRITE(11,'(A)') "The pressure at the corner of the box is: "
   WRITE(11,'(A)') "Freedom  Pressure "
   WRITE(11,'(I6,E12.4)') nres, x_pp(is)
   WRITE(11,'(A,I6)')"The total number of BiCGSTAB iterations was:",cj_tot
   WRITE(11,'(A,I5,A)')"The solution took",iters," iterations to converge"
   WRITE(11,'(A,F10.4)')"Time spent in solver was:",timest(4)-timest(3)
 END IF
 IF(numpe==1) THEN
   OPEN(12,file=argv(1:nlen)//".ensi.VEL",status='replace',              &
     action='write')
   WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
   WRITE(12,'(A/A/A)') "part", "     1","coordinates"
 END IF
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 ALLOCATE(disp_pp(nodes_pp*nodof),temp(nodes_pp)); disp_pp=zero          
 temp=zero; utemp_pp=zero; CALL gather(x_pp(1:),utemp_pp)
 CALL scatter_nodes_ns(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,      &
   node_start,node_end,utemp_pp,disp_pp,1)
 DO i=1,nodof ; temp=zero
   IF(i/=2) THEN
     DO j=1,nodes_pp; k=i+(nodof*(j-1)); temp(j)=disp_pp(k); END DO
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp); END IF
 END DO ; IF(numpe==1) CLOSE(12)
 IF(numpe==it) THEN 
   WRITE(11,'(A,F10.4)') "This analysis took :", elap_time()-timest(1)
   CLOSE(11); END IF; CALL shutdown()    
END PROGRAM p126
