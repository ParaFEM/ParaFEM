PROGRAM p127
!-------------------------------------------------------------------------
!      Program 12.7 3-D consolidation of a cuboidal Biot elastic
!      solid using 20-node solid hexahedral elements  coupled to 8-node
!      fluid elements-parallel pcg version - biot_cube
!-------------------------------------------------------------------------
!USE mpi_wrapper  !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE loading
 USE timing; USE maths; USE gather_scatter; USE new_library; USE geometry
 USE input; IMPLICIT NONE
! neq,ntot are now global variables - not declared
 INTEGER::nxe,nye,nze,nn,nr,nip,nodof=4,nod=20,nodf=8,nst=6,ndim=3,i,j,k,&
   l,iel,ns,nstep,cjiters,cjits,loaded_freedoms,loaded_freedoms_pp,n_t,  &
   neq_temp,nn_temp,nle,nlen,nels,partitioner=1,ndof,ielpe,npes_pp,      &
   loaded_freedoms_start,nlfp,nls
 REAL(iwp)::kx,ky,kz,e,v,det,dtim,theta,real_time,up,alpha,beta,cjtol,aa,&
   bb,cc,q
 REAL(iwp),PARAMETER::zero=0._iwp
 LOGICAL::cj_converged; CHARACTER(LEN=15)::element='hexahedron'
 CHARACTER(LEN=50)::argv
!---------------------------- dynamic arrays------------------------------
 REAL(iwp),ALLOCATABLE::dee(:,:),points(:,:),coord(:,:),derivf(:,:),     &
   jac(:,:),kay(:,:),der(:,:),deriv(:,:),weights(:),derf(:,:),funf(:),   &
   coordf(:,:),bee(:,:),km(:,:),eld(:),sigma(:),kc(:,:),ke(:,:),         &
   g_coord_pp(:,:,:),kd(:,:),fun(:),c(:,:),loads_pp(:),pmul_pp(:,:),     &
   storke_pp(:,:,:),ans_pp(:),volf(:,:),p_pp(:),x_pp(:),xnew_pp(:),      &
   u_pp(:),eld_pp(:,:),diag_precon_pp(:),diag_precon_tmp(:,:),d_pp(:),   &
   utemp_pp(:,:),storkc_pp(:,:,:),val(:),vol(:),timest(:),al(:),lf(:,:)
 INTEGER,ALLOCATABLE::rest(:,:),g(:),num(:),g_g_pp(:,:),g_num_pp(:,:),   &
   g_t(:),no(:),no_pp_temp(:),no_pp(:)
!-------------------------input and initialisation------------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)
 CALL read_p127(argv,numpe,nels,nxe,nze,aa,bb,cc,nip,kx,ky,kz,e,v,dtim,  &
   nstep,theta,cjits,cjtol)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*ndim; ntot=ndof+nodf; n_t=nod*nodof
 nye=nels/nxe/nze; nle=nxe/5; loaded_freedoms=3*nle*nle+4*nle+1
 nr=3*nxe*nye*nze+4*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+2
 ALLOCATE(points(nip,ndim),coord(nod,ndim),derivf(ndim,nodf),            &
   jac(ndim,ndim),kay(ndim,ndim),der(ndim,nod),deriv(ndim,nod),          &
   derf(ndim,nodf),funf(nodf),coordf(nodf,ndim),bee(nst,ndof),           &
   km(ndof,ndof),eld(ndof),sigma(nst),kc(nodf,nodf),weights(nip),        &
   g_g_pp(ntot,nels_pp),diag_precon_tmp(ntot,nels_pp),ke(ntot,ntot),     &
   kd(ntot,ntot),fun(nod),c(ndof,nodf),g_t(n_t),vol(ndof),               &
   rest(nr,nodof+1),g(ntot),volf(ndof,nodf),g_coord_pp(nod,ndim,nels_pp),&
   g_num_pp(nod,nels_pp),num(nod),storke_pp(ntot,ntot,nels_pp),          &
   storkc_pp(nodf,nodf,nels_pp),pmul_pp(ntot,nels_pp),dee(nst,nst),      &
   utemp_pp(ntot,nels_pp),eld_pp(ntot,nels_pp),no(loaded_freedoms),      &
   val(loaded_freedoms),no_pp_temp(loaded_freedoms))
 kay=0.0_iwp; kay(1,1)=kx; kay(2,2)=ky; kay(3,3)=kz
 CALL biot_cube_bc20(nxe,nye,nze,rest); CALL rearrange(rest)
 CALL biot_loading(nxe,nze,nle,no,val); val=-val*aa*bb/12._iwp
 CALL sample(element,points,weights); CALL deemat(dee,e,v)
 ielpe=iel_start
 nlfp=3; ALLOCATE(lf(2,nlfp)); lf=zero
 lf(1,1)=0.0_iwp; lf(2,1)=0.0_iwp
 lf(1,2)=10.0_iwp; lf(2,2)=1.0_iwp
 lf(1,3)=20.0_iwp; lf(2,3)=1.0_iwp
 nls=FLOOR(lf(1,nlfp)/dtim); IF(nstep>nls)nstep=nls
 ALLOCATE(al(nstep)); al=zero; CALL load_function(lf,dtim,al)
 IF(numpe==1) PRINT *, al
!----------------- loop the elements to set up global arrays -------------
 neq_temp=0; nn_temp=0
 elements_1: DO iel=1,nels_pp
   CALL geometry_20bxz(ielpe,nxe,nze,aa,bb,cc,coord,num)
   CALL find_g3(num,g_t,rest); CALL g_t_g(nod,g_t,g)
   g_coord_pp(:,:,iel)=coord; g_g_pp(:,iel)=g; ielpe=ielpe+1
   i=MAXVAL(g); j=MAXVAL(num); g_num_pp(:,iel)=num
   IF(i>neq_temp)neq_temp=i; IF(j>nn_temp)nn_temp=j
 END DO elements_1
 neq=max_p(neq_temp); nn=max_p(nn_temp); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
 IF(numpe==1)THEN
   OPEN(11,FILE=argv(1:nlen)//'.res',STATUS='REPLACE',ACTION='WRITE') 
   WRITE(11,'(A,I5,A)')"This job ran on ",npes,"  processes"
   WRITE(11,'(A)')"Global coordinates and node numbers"
   DO i=1,nels_pp,nels_pp-1
     WRITE(11,'(A,I8)')"Element ",i
     num=g_num_pp(:,i)
     DO k=1,nod
       WRITE(11,'(A,I8,3E12.4)')"   Node",num(k),g_coord_pp(k,:,i)
     END DO
   END DO
   WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nr," restrained and ",&
     neq,"  equations"
   WRITE(11,'(A,F10.4)')"Time after setup is:",elap_time()-timest(1)
 END IF
 ALLOCATE(loads_pp(neq_pp),ans_pp(neq_pp),p_pp(neq_pp),x_pp(neq_pp),     &
   xnew_pp(neq_pp),u_pp(neq_pp),diag_precon_pp(neq_pp),d_pp(neq_pp))
 loads_pp=.0_iwp; p_pp=.0_iwp; xnew_pp=.0_iwp; diag_precon_pp=.0_iwp
 diag_precon_tmp=.0_iwp
!-------- element stiffness integration, storage and preconditioner ------ 
 elements_2: DO iel=1,nels_pp 
   coord=g_coord_pp(:,:,iel); coordf(1:4,:)=coord(1:7:2,:)
   coordf(5:8,:)=coord(13:20:2,:); km=zero; c=zero; kc=zero
   gauss_points_1: DO i=1,nip
     CALL shape_der(der,points,i); jac=MATMUL(der,coord)
     det=determinant(jac); CALL invert(jac)
     deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
     vol(:)=bee(1,:)+bee(2,:)+bee(3,:)
     km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
!------------------------- now the fluid contribution --------------------
     CALL shape_fun(funf,points,i); CALL shape_der(derf,points,i)
     derivf=MATMUL(jac,derf)
     kc=kc+                                                              &
        MATMUL(MATMUL(TRANSPOSE(derivf),kay),derivf)*det*weights(i)*dtim
     DO l=1,nodf; volf(:,l)=vol(:)*funf(l); END DO
     c=c+volf*det*weights(i)               
   END DO gauss_points_1
   storkc_pp(:,:,iel)=kc; CALL formke(km,kc,c,ke,theta)
   storke_pp(:,:,iel)=ke 
   DO k=1,ndof
     diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)+theta*km(k,k)
   END DO
   DO k=1,nodf
     diag_precon_tmp(ndof+k,iel)=diag_precon_tmp(ndof+k,iel)             &
       -theta*theta*kc(k,k)
   END DO
 END DO elements_2
 CALL scatter(diag_precon_pp,diag_precon_tmp)
 diag_precon_pp=1._iwp/diag_precon_pp
 DEALLOCATE(diag_precon_tmp)
!----------------------------loaded freedoms -----------------------------
 no_pp_temp=0
 CALL reindex(ieq_start,no,no_pp_temp,loaded_freedoms_pp,                &
   loaded_freedoms_start,neq_pp)
 ALLOCATE(no_pp(1:loaded_freedoms_pp)); no_pp=0
 no_pp=no_pp_temp(1:loaded_freedoms_pp); DEALLOCATE(no_pp_temp)
!----------------------- enter the time-stepping loop --------------------
 real_time=zero     
 time_steps: DO ns=1,nstep
   real_time=real_time+dtim
   IF(numpe==1) WRITE(11,'(A,E12.4)')"The time is", real_time
   pmul_pp=zero; utemp_pp=zero; ans_pp=zero
   CALL gather(loads_pp,pmul_pp)
   elements_3: DO iel=1,nels_pp
     utemp_pp(ndof+1:,iel)=MATMUL(storkc_pp(:,:,iel),pmul_pp(ndof+1:,iel))
   END DO elements_3; CALL scatter(ans_pp,utemp_pp)
!---------------------------- ramp loading -------------------------------
   IF(loaded_freedoms_pp>0) THEN
     DO i=1,loaded_freedoms_pp; j=no_pp(i)-ieq_start+1
       ans_pp(j)=val(loaded_freedoms_start+i-1)*al(ns)
     END DO
   END IF
   d_pp=diag_precon_pp*ans_pp; p_pp=d_pp
   x_pp=zero  ! depends on starting x = zero
!----------------- solve the simultaneous equations by pcg ---------------
   cjiters=0
   conjugate_gradients: DO 
     cjiters=cjiters+1; u_pp=zero; pmul_pp=zero; utemp_pp=zero
     CALL gather(p_pp,pmul_pp)
     elements_4: DO iel=1,nels_pp
       utemp_pp(:,iel)=MATMUL(storke_pp(:,:,iel),pmul_pp(:,iel))
     END DO elements_4; CALL scatter(u_pp,utemp_pp)
!--------------------------- pcg process ---------------------------------
     up=DOT_PRODUCT_P(ans_pp,d_pp); alpha=up/DOT_PRODUCT_P(p_pp,u_pp)
     xnew_pp=x_pp+p_pp*alpha; ans_pp=ans_pp-u_pp*alpha
     d_pp=diag_precon_pp*ans_pp; beta=DOT_PRODUCT_P(ans_pp,d_pp)/up
     p_pp=d_pp+p_pp*beta
     CALL checon_par(xnew_pp,cjtol,cj_converged,x_pp)
     IF(cj_converged.OR.cjiters==cjits)EXIT
   END DO conjugate_gradients
   ans_pp=xnew_pp; loads_pp=loads_pp+ans_pp
   IF(numpe==1)THEN
     WRITE(11,'(A,I5,A)')                                                &
       "Conjugate gradients took ",cjiters,"  iterations to converge"
     WRITE(11,'(A)')" The nodal displacements and pore pressures are    :"
     WRITE(11,'(4E12.4)')loads_pp(1:4)
   END IF
!------------------ recover stresses at gauss-points ---------------------
   eld_pp=zero; CALL gather(loads_pp,eld_pp)
   iel=1; coord=g_coord_pp(:,:,iel); eld=eld_pp(:,iel)
   IF(numpe==1)WRITE(11,'(A,I5,A)')                                      &
     "The Gauss Point effective stresses for element",iel,"  are"
   gauss_pts_2: DO i=1,nip
     CALL shape_der(der,points,i); jac=MATMUL(der,coord)
     CALL invert(jac); deriv=MATMUL(jac,der)
     CALL beemat(bee,deriv); sigma=MATMUL(dee,MATMUL(bee,eld))
     IF(numpe==1.AND.i==1)THEN
       WRITE(11,'(A,I5)')"Point  ",i
       WRITE(11,'(6E12.4)')sigma
     END IF
   END DO gauss_pts_2 
 END DO time_steps
 IF(numpe==1) THEN; WRITE(11,'(A,F10.4)')"This analysis took:",          &
   elap_time()-timest(1); CLOSE(11); END IF; CALL SHUTDOWN()
END PROGRAM p127
