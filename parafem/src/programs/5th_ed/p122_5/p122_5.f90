PROGRAM p122      
!-----------------------------------------------------------------------
!    Program 6.11 3-d strain of an elastic-plastic(Mohr-Coulomb) solid
!    using 20-node brick elements: viscoplastic strain method,pcg parallel
!    cube_bc  : pick up current x not x = .0  : load control
!------------------------------------------------------------------------
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE elements; USE steering; USE plasticity; IMPLICIT NONE
!--------------- ndof,nels,neq,ntot are global - must not be declared-----
 INTEGER::nxe,nye,nze,nn,nr,nip,nodof=3,nod=20,nst=6,i,j,k,iel,plasiters, &
   plasits,cjiters,cjits,cjtot,incs,iy,ndim=3,loaded_freedoms,num_no,     &
   no_index_start,neq_temp,nn_temp,nle,nlen
 LOGICAL::plastic_converged,cj_converged
 CHARACTER(LEN=15)::element='hexahedron',argv;  CHARACTER(LEN=5)::ch
 REAL(iwp)::e,v,det,phi,c,psi,dt,f,dsbar,dq1,dq2,dq3,lode_theta,presc,    &
   sigm,pi,snph,cons,aa,bb,cc,plastol,cjtol,up,alpha,beta,big
 REAL(iwp)::zero=0.0_iwp
!---------------------------- dynamic arrays-----------------------------
 REAL(iwp),ALLOCATABLE::loads_pp(:),points(:,:),bdylds_pp(:),             &
   evpt_pp(:,:,:),pmul_pp(:,:),dee(:,:),coord(:,:),jac(:,:),weights(:),   &
   oldis_pp(:),der(:,:),deriv(:,:),bee(:,:),km(:,:),eld(:),eps(:),        &
   sigma(:),bload(:),eload(:),erate(:),p_g_co_pp(:,:,:),evp(:),devp(:),   &
   m1(:,:),m2(:,:),m3(:,:),flow(:,:),storkm_pp(:,:,:),r_pp(:),val(:),     &
   tensor_pp(:,:,:),stress(:),totd_pp(:),qinc(:),p_pp(:),x_pp(:),         &
   xnew_pp(:),u_pp(:),utemp_pp(:,:),diag_precon_pp(:),d_pp(:),            &
   diag_precon_tmp(:,:),timest(:)
 INTEGER,ALLOCATABLE::rest(:,:),g(:),no(:),num(:),g_num_pp(:,:),          &
   g_g_pp(:,:),no_local(:),no_local_temp(:)
!--------------------------input and initialisation----------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen) 
!IF(numpe==npes)THEN
!  CALL getname(argv,nlen)
!  OPEN(10,FILE=argv(1:nlen)//'.dat',STATUS='OLD',ACTION='READ')
!  READ(10,*)phi,c,psi,e,v,cons, nels,nxe,nze,nip,aa,bb,cc,incs,plasits,  &
!  cjits,plastol,cjtol
!END IF
 CALL read_p122(argv,numpe,phi,c,psi,e,v,cons,nels,nip,incs,plasits,      &
   cjits,plastol,cjtol)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof
 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),           &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp2(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest); timest(2)=elap_time()
 ALLOCATE(,points(nip,ndim),weights(nip),m1(nst,nst),num(nod),            &
   dee(nst,nst),evpt_pp(nst,nip,nels_pp),m2(nst,nst),                     &
   tensor_pp(nst,nip,nels_pp),coord(nod,ndim),g_g_pp(ntot,nels_pp),       &
   stress(nst),storkm_pp(ntot,ntot,nels_pp),jac(ndim,ndim),der(ndim,nod), &
   deriv(ndim,nod),bee(nst,ntot),km(ntot,ntot),eld(ntot),eps(nst),        &
   sigma(nst),bload(ntot),eload(ntot),erate(nst),evp(nst),devp(nst),      &
   g(ntot),m3(nst,nst),flow(nst,nst),pmul_pp(ntot,nels_pp),               &
   utemp_pp(ntot,nels_pp),no(loaded_freedoms),                            &
   no_local_temp(loaded_freedoms),val(loaded_freedoms),qinc(incs))
 IF(numpe==1)READ(10,*)qinc
 CALL MPI_BCAST(qinc,incs,MPI_REAL8,npes-1,MPI_COMM_WORLD,ier)
 CALL rearrange(rest)
!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_1: DO iel=1,nels_pp
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
!   CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_1
 elements_2: DO iel=1,nels_pp  
    i=MAXVAL(g_g_pp(:,iel)); IF(i>neq) neq=i
 END DO elements_2  
 neq=MAX_INTEGER_P(neq); CALL calc_neq_pp; CALL calc_npes_pp(npes,npes_pp)
 CALL make_ggl2(npes_pp,npes,g_g_pp)
 IF(numpe==1)THEN
   OPEN(11,FILE=argv(1:nlen)//'.res',STATUS='REPLACE',ACTION='WRITE')              
   WRITE(11,'(A,I5,A)')"This job ran on ", npes,"  processes"
   WRITE(11,'(A,3(I5,A))')"There are ",nn," nodes",nr,                    &
     " restrained and   ",neq," equations"
   WRITE(11,*)"Time after setup  is  : ",elap_time()-timest(1)
 END IF 
 ALLOCATE(loads_pp(neq_pp),bdylds_pp(neq_pp),oldis_pp(neq_pp),            &
   r_pp(neq_pp),totd_pp(neq_pp),p_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),&        
   u_pp(neq_pp),diag_precon_pp(neq_pp),d_pp(neq_pp));  oldis_pp=zero
 totd_pp=zero; tensor_pp=zero; p_pp=zero; xnew_pp=zero; diag_precon_pp=zero
 pi=ACOS(-1._iwp )
 snph=SIN(phi*pi/180._iwp)
 dt=4._iwp*(1._iwp+v)*(1._iwp-2._iwp*v)/(e*(1._iwp-2._iwp*v+snph*snph))
 IF(numpe==1)WRITE(11,'(A,E12.4)')"The critical timestep is   ",dt
!---- element stiffness integration, preconditioner & set initial stress-- 
 CALL deemat(e,v,dee); CALL sample(element,points,weights)
 storkm_pp=zero
 elements_2: DO iel=1,nels_pp
   gauss_pts_1: DO i=1,nip    
     tensor_pp(1:3,i,iel)=cons; CALL shape_der(der,points,i)
     jac=MATMUL(der,g_coord_pp(:,:,iel); det=determinant(jac)
     CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
     storkm_pp(:,:,iel)=storkm_pp(:,:,iel)+                              &
                    MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
   END DO gauss_pts_1
 END DO elements_2
 ALLOCATE(diag_precon_tmp(ntot,nels_pp)); diag_precon_tmp=zero
 elements_4: DO iel=1,nels_pp ; DO i=1,ndof
   diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel)+storkm_pp(i,i,iel)
 END DO;  END DO elements_4
 CALL scatter(diag_precon_pp,diag_precon_tmp); DEALLOCATE(diag_precon_tmp)
!-------------------- invert preconditioner -----------------------------
 CALL reindex_fixed_nodes                                                 &
   (ieq_start,no,no_local_temp,num_no,no_index_start)
 ALLOCATE(no_local(1:num_no)); no_local=no_local_temp(1:num_no)
 DEALLOCATE(no_local_temp); diag_precon_pp=1._iwp/diag_precon_pp
!------------------------  load  increment loop --------------------------
 load_increments: do iy=1,incs
   plasiters=0; bdylds_pp=zero; evpt_pp=zero; cjtot=0
   IF(numpe==npes)WRITE(11,'(/,A,I5)')"Load Increment   ",iy
!--------------------------   plastic iteration loop   -------------------
   plastic_iterations: DO
     plasiters=plasiters+1
     loads_pp=zero
     DO i=1,num_no
       j= no_local(i)-ieq_start+1
       loads_pp(j)=val(no_index_start+i-1)*qinc(iy) 
     END DO
     loads_pp=loads_pp+bdylds_pp
!------ if x=.0 p and r are just loads but in general p=r=loads-A*x ------
!-----------------------so form r = A * x --------------------------------
     r_pp=zero
     CALL gather(x_pp,pmul_pp)
     elements_2a: DO iel=1,nels_pp
       utemp_pp(:,iel)=MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
     END DO elements_2a
     CALL scatter(r_pp,utemp_pp)
!------------------------now precondition r and p ------------------------
     r_pp=loads_pp-r_pp; d_pp=diag_precon_pp*r_pp; p_pp=d_pp
!-------------------   solve the simultaneous equations by pcg -----------
     cjiters=0
     conjugate_gradients: DO
       cjiters=cjiters+1
       u_pp=zero; pmul_pp=zero
       CALL gather(p_pp,pmul_pp)
       elements_3: DO iel=1,nels_pp
         utemp_pp(:,iel)=MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
       END DO elements_3
       CALL scatter(u_pp,utemp_pp)
!-------------------------------pcg process ------------------------------
       up=DOT_PRODUCT_P(r_pp,d_pp); alpha=up/DOT_PRODUCT_P(p_pp,u_pp)
       xnew_pp=x_pp+p_pp*alpha; r_pp=r_pp-u_pp*alpha
       d_pp=diag_precon_pp*r_pp; beta=DOT_PRODUCT_P(r_pp,d_pp)/up
       p_pp=d_pp+p_pp*beta; cj_converged=.TRUE.
       CALL checon_par(xnew_pp,x_pp,cjtol,cj_converged,neq_pp)
       IF(cj_converged.or.cjiters==cjits)EXIT
     END DO conjugate_gradients
     cjtot=cjtot+cjiters  
!---------------------------- end of pcg process -------------------------
     loads_pp=xnew_pp; pmul_pp=zero
!-----------------------   check plastic convergence  --------------------
     CALL checon_par(loads_pp,oldis_pp,plastol,plastic_converged,neq_pp)
     IF(plasiters==1)plastic_converged=.FALSE. 
     IF(plastic_converged.OR.plasiters==plasits)bdylds_pp=zero
     CALL gather(loads_pp,pmul_pp); utemp_pp=zero
!------------------------ go round the Gauss Points ----------------------
     elements_4: DO iel=1,nels_pp
       bload=zero; eld=pmul_pp(:,iel)
       gauss_points_2: DO i=1,nip
         CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
         det=determinant(jac); CALL invert(jac)
         deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
         eps=MATMUL(bee,eld); eps=eps-evpt_pp(:,i,iel)
         sigma=MATMUL(dee,eps); stress=sigma+tensor_pp(:,i,iel)
         CALL invar(stress,sigm,dsbar,lode_theta)                            
!----------------------  check whether yield is violated  ----------------
         CALL mocouf(phi,c,sigm,dsbar,lode_theta,f)
         IF(plastic_converged.OR.plasiters==plasits)THEN
           devp=stress 
         ELSE
           IF(f>=zero)THEN
             CALL mocouq(psi,dsbar,lode_theta,dq1,dq2,dq3)
             CALL formm(stress,m1,m2,m3); flow=f*(m1*dq1+m2*dq2+m3*dq3)
             erate=MATMUL(flow,stress); evp=erate*dt
             evpt_pp(:,i,iel)=evpt_pp(:,i,iel)+evp
             devp=MATMUL(dee,evp)
           END IF
         END IF
         IF(f>=zero)THEN
           eload=MATMUL(TRANSPOSE(bee),devp)
           bload=bload+eload*det*weights(i)
         END IF
         IF(plastic_converged.OR.plasiters==plasits)THEN                          
!---------------------- update the Gauss Point stresses  -----------------
           tensor_pp(:,i,iel)=stress          
         END IF
       END DO gauss_points_2
!-------------compute the total bodyloads vector -------------------------
       utemp_pp(:,iel)=utemp_pp(:,iel)+bload
     END DO elements_4
     CALL scatter(bdylds_pp,utemp_pp)
     IF(plastic_converged.OR.plasiters==plasits)EXIT
   END DO plastic_iterations
   totd_pp=totd_pp+loads_pp
   IF(numpe==1)THEN
     WRITE(11,'(A,E12.4)')"The displacement is  ",totd_pp(1)
     WRITE(11,'(A)')"  sigma z    sigma x     sigma y"
     WRITE(11,'(3E12.4)')tensor_pp(3,1,1),tensor_pp(1,1,1),tensor_pp(2,1,1)
     WRITE(11,'(A,I12)')"The total number of cj iterations was  ",cjtot 
     WRITE(11,'(A,I12)')"The number of plastic iterations was   ",plasiters 
     WRITE(11,'(A,F11.2)')"cj iterations per plastic iteration were ",    &
       REAL(cjtot)/REAL(plasiters)
   END IF
   IF(plasiters==plasits)EXIT
 END DO load_increments
 IF(numpe==1)WRITE(11,*)"This analysis took : ",elap_time()-timest(1)
 CALL shutdown()
END PROGRAM p122
