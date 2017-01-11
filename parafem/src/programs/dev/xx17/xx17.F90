PROGRAM xx17    
!-------------------------------------------------------------------------
!      Program xx17 steady state  3-d Navier-Stokes equation
!      using 20-node velocity hexahedral elements  
!      coupled to 8-node pressure hexahedral elements : u-p-v-w order
!      element by element solution using BiCGSTAB(L) : parallel version
!-------------------------------------------------------------------------
!USE mpi_wrapper  !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE input; 
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE steering; USE fluid; USE new_library; USE bicg
 IMPLICIT NONE
!neq,ntot are now global variables - not declared
 INTEGER,PARAMETER::nodof=4,nod=20,nodf=8,ndim=3
 INTEGER::nn,nip,cj_tot,i,j,k,l,iel,ell,limit,fixed_freedoms,iters,      &
   cjits,nr,n_t,fixed_freedoms_pp,nres,is,it,nlen,nels,ndof,     &
   npes_pp,meshgen,partitioner,node_end,node_start,nodes_pp,             &
   fixed_freedoms_start,cjiters=0  
 REAL(iwp):: visc,rho,det,ubar,vbar,wbar,tol,cjtol,alpha,      &
   penalty,x0,pp,kappa
 REAL(iwp),PARAMETER::zero=0.0_iwp,one=1.0_iwp
 LOGICAL::converged
 CHARACTER(LEN=15)::element='hexahedron'; CHARACTER(LEN=50)::argv
!--------------------------- dynamic arrays ------------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),derivf(:,:),fun(:),store_pp(:),      &
   jac(:,:),kay(:,:),der(:,:),deriv(:,:),weights(:),derf(:,:),funf(:),   &
   coordf(:,:),g_coord_pp(:,:,:),c11(:,:),c21(:,:),c12(:,:),val(:),      &
   wvel(:),c23(:,:),c32(:,:),x_pp(:),b_pp(:),temp(:),                    &
   funny(:,:),row1(:,:),row2(:,:),uvel(:),vvel(:),funnyf(:,:),rowf(:,:), &
   storke_pp(:,:,:),diag_pp(:),utemp_pp(:,:),xold_pp(:),c24(:,:),        &
   c42(:,:),row3(:,:),                                                   &
   diag_tmp(:,:),pmul_pp(:,:),timest(:),upvw_pp(:)
 INTEGER,ALLOCATABLE::rest(:,:),g_num_pp(:,:),g_g_pp(:,:),no(:),g_t(:),  &
   no_pp(:),no_pp_temp(:)
!---------------------- input and initialisation -------------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)
 CALL read_p126(argv,numpe,cjits,cjtol,ell,fixed_freedoms,kappa,limit,   &
   meshgen,nels,nip,nn,nr,nres,partitioner,penalty,rho,tol,x0,visc) 
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ntot=nod+nodf+nod+nod; n_t=nod*nodof
 ALLOCATE(g_num_pp(nod,nels_pp),g_coord_pp(nod,ndim,nels_pp),            &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest); timest(2)=elap_time()
 ALLOCATE(points(nip,ndim),derivf(ndim,nodf),pmul_pp(ntot,nels_pp),      &
   jac(ndim,ndim),kay(ndim,ndim),der(ndim,nod),deriv(ndim,nod),          &
   derf(ndim,nodf),funf(nodf),coordf(nodf,ndim),funny(nod,1),            &
   g_g_pp(ntot,nels_pp),c11(nod,nod),c12(nod,nodf),c21(nodf,nod),        &
   c24(nodf,nod),c42(nod,nodf),c32(nod,nodf),fun(nod),row2(1,nod),       &
   c23(nodf,nod),uvel(nod),vvel(nod),row1(1,nod),funnyf(nodf,1),         &
   rowf(1,nodf),no_pp_temp(fixed_freedoms),wvel(nod),row3(1,nod),        &
   storke_pp(ntot,ntot,nels_pp),g_t(n_t),                                &
   no(fixed_freedoms),val(fixed_freedoms),weights(nip),                  &
   diag_tmp(ntot,nels_pp),utemp_pp(ntot,nels_pp))
!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_1: DO iel=1,nels_pp
   CALL find_g3(g_num_pp(:,iel),g_t,rest)
   CALL g_t_g_ns(nod,g_t,g_g_pp(:,iel))
 END DO elements_1
 neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
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
 ALLOCATE(x_pp(neq_pp),b_pp(neq_pp),                                     &
   diag_pp(neq_pp),xold_pp(neq_pp),                                      &
   store_pp(neq_pp))
 x_pp=zero; b_pp=zero; diag_pp=zero
 xold_pp=zero; store_pp=zero
!-------------------------- organise fixed equations ---------------------
 CALL read_loads_ns(argv,numpe,no,val)
 CALL reindex(ieq_start,no,no_pp_temp,fixed_freedoms_pp,                &
   fixed_freedoms_start,neq_pp); ALLOCATE(no_pp(1:fixed_freedoms_pp))
 no_pp = no_pp_temp(1:fixed_freedoms_pp); DEALLOCATE(no_pp_temp)
!------------------------- main iteration loop ---------------------------
 CALL sample(element,points,weights); uvel=zero; vvel=zero; wvel=zero
 kay=zero; iters=0; cj_tot=0; kay(1,1)=visc/rho; kay(2,2)=visc/rho
 kay(3,3)=visc/rho; timest(3)=elap_time()
 iterations: DO
   iters=iters+1; storke_pp=zero; diag_pp=zero; utemp_pp=zero
   b_pp=zero; pmul_pp=zero; CALL gather(x_pp,utemp_pp)
   CALL gather(xold_pp,pmul_pp)
!-------------------- element stiffness integration ----------------------
   elements_2: DO iel=1,nels_pp
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
   END DO elements_2
!----------------------- build the preconditioner ------------------------
 diag_tmp=zero
 elements_2a: DO iel=1,nels_pp; DO k=1,ntot
   diag_tmp(k,iel)=diag_tmp(k,iel)+storke_pp(k,k,iel); END DO
 END DO elements_2a; CALL scatter(diag_pp,diag_tmp)
!------------------- prescribed values of velocity and pressure ----------
 DO i=1,fixed_freedoms_pp; k=no_pp(i)-ieq_start+1
   diag_pp(k)=diag_pp(k)+penalty
   b_pp(k)=diag_pp(k)*val(fixed_freedoms_start+i-1)
   store_pp(k)=diag_pp(k)
 END DO   
!---------- solve the equations element-by-element using BiCGSTAB --------
 IF(iters==1) x_pp = x0

 CALL bicgstabl_p(b_pp,cjiters,cjits,cjtol,ell,kappa,ieq_start,no_pp,          &
                  pmul_pp,store_pp,storke_pp,x_pp)

 b_pp=x_pp-xold_pp; pp=norm_p(b_pp); cj_tot=cj_tot+cjiters

 IF(numpe==it) THEN
   WRITE(11,'(A,E12.4)') "Norm of the error is:", pp
   WRITE(11,'(A,I6,A)') "It took BiCGSTAB(L) ", cjiters,               &
      " iterations to converge"; END IF
 CALL checon_par(x_pp,tol,converged,xold_pp)
 IF(converged.OR.iters==limit)EXIT

 END DO iterations

 timest(4)=elap_time()
 DEALLOCATE(b_pp,diag_pp,xold_pp,store_pp)
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
 ALLOCATE(upvw_pp(nodes_pp*nodof),temp(nodes_pp)); upvw_pp=zero          
 temp=zero; utemp_pp=zero; CALL gather(x_pp(1:),utemp_pp)
 CALL scatter_nodes_ns(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,      &
   node_start,node_end,utemp_pp,upvw_pp,1)
 DO i=1,nodof ; temp=zero
   IF(i/=2) THEN
     DO j=1,nodes_pp; k=i+(nodof*(j-1)); temp(j)=upvw_pp(k); END DO
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp); END IF
 END DO ; IF(numpe==1) CLOSE(12)
 IF(numpe==it) THEN 
   WRITE(11,'(A,F10.4)') "This analysis took :", elap_time()-timest(1)
   CLOSE(11); END IF; CALL SHUTDOWN() 

END PROGRAM xx17
