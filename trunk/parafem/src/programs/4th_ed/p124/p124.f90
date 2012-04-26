  PROGRAM p124      
!-------------------------------------------------------------------------
!      Program 8.3 conduction equation on 3-d box shaped 
!      volume using 8-node hexahedral elements : pcg version implicit
!      integration in time using 'theta' method : parallel version : box_bc
!-------------------------------------------------------------------------
 USE new_library
 USE geometry_lib
 USE precision
 USE utility
 USE mp_module
 USE timing
 USE global_variables1
 USE gather_scatter6
 IMPLICIT NONE
!----- ndof,nels,neq,ntot are now global variables - not declared---------
 INTEGER::nxe,nye,nze,nn,nr,nip,nodof=1,nod=8,ndim=3,neq_temp,nn_temp,i,j,&
   k,iel,nstep,npri,nres,iters,limit,it,is,nlen
 REAL(iwp)::aa,bb,cc,kx,ky,kz,det,theta,dtim,val0,real_time,tol,alpha,    &
   beta,up,big
 LOGICAL::converged
 CHARACTER(LEN=15)::element='hexahedron',argv
!------------------------- dynamic arrays----------------------------------
 REAL(iwp),ALLOCATABLE::loads_pp(:),u_pp(:),p_pp(:),points(:,:),kay(:,:), &
   coord(:,:),fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),d_pp(:),     &
   kc(:,:),pm(:,:),funny(:,:),p_g_co_pp(:,:,:),storka_pp(:,:,:),          &
   storkb_pp(:,:,:),x_pp(:),xnew_pp(:),pmul_pp(:,:),utemp_pp(:,:),        &
   diag_precon_pp(:),diag_precon_tmp(:,:)
 INTEGER,ALLOCATABLE::rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:)
!----------------------input and initialisation---------------------------
 timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes)
 IF(numpe==npes)THEN
   CALL getname(argv,nlen)
   OPEN(10,FILE=argv(1:nlen)//'.dat',STATUS='OLD',ACTION='READ')
   READ(10,*)nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz,dtim,nstep,theta,npri,tol,&
     limit,val0
 END IF
 CALL bcast_inputdata_p124(numpe,npes,nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz, &
   dtim,nstep,theta,npri,tol,limit,val0)
 CALL calc_nels_pp
 ndof=nod*nodof
 ntot=ndof
 nn_temp=0
 neq_temp=0
 nye=nels/nxe/nze
 nr=(nxe+1)*(nye+1)+(nxe+1)*nze+nye*nze
 ielpe=iel_start
 ALLOCATE(rest(nr,nodof+1),points(nip,ndim),weights(nip),kay(ndim,ndim),  &
   coord(nod,ndim),fun(nod),jac(ndim,ndim),der(ndim,nod),g(ntot),         &
   p_g_co_pp(nod,ndim,nels_pp),deriv(ndim,nod),pm(ntot,ntot),             &
   g_num_pp(nod,nels_pp),kc(ntot,ntot),funny(1,nod),num(nod),             &
   g_g_pp(ntot,nels_pp),storka_pp(ntot,ntot,nels_pp),                     &
   utemp_pp(ntot,nels_pp),storkb_pp(ntot,ntot,nels_pp),                   &
   pmul_pp(ntot,nels_pp),diag_precon_tmp(ntot,nels_pp))
 kay=0.0_iwp
 kay(1,1)=kx
 kay(2,2)=ky
 kay(3,3)=kz
 CALL sample(element,points,weights)
 CALL box_bc8(nxe,nye,nze,rest)
 CALL rearrange_2(rest)
!-------------loop the elements to  set up global arrays -----------------
 elements_1: DO iel=1,nels_pp
   CALL geometry_8bxz(ielpe,nxe,nze,aa,bb,cc,coord,num)
   CALL find_g4(num,g,rest)
   g_num_pp(:,iel)=num
   p_g_co_pp(:,:,iel)=coord
   g_g_pp(:,iel)=g
   ielpe=ielpe+1
   i=MAXVAL(g)
   j=MAXVAL(num)
   IF(i>neq_temp)neq_temp=i
   IF(j>nn_temp)nn_temp=j
 END DO elements_1
 neq=reduce(neq_temp)
 nn=reduce(nn_temp)
 nres=nxe*(nze-1)+1
 CALL calc_neq_pp
 CALL make_ggl(g_g_pp)
 DO i=1,neq_pp
   IF(nres==ieq_start+i-1)THEN
     it=numpe
     is=i
   END IF
 END DO
 IF(numpe==it)THEN
   OPEN(11,FILE=argv(1:nlen)//'.res',STATUS='REPLACE',ACTION='WRITE')              
   WRITE(11,'(A,I5,A)')"This job ran on  ",npes,"  processors"  
   WRITE(11,'(A)')"Global coordinates and node numbers"
   DO i=1,nels_pp,nels_pp - 1                                        
     WRITE(11,'(A,I8)')"Element ",i
     num=g_num_pp(:,i)
     DO k=1,nod
       WRITE(11,'(A,I8,3E12.4)')"  Node",num(k),p_g_co_pp(k,:,i)
     END DO
   END DO
   WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nr," restrained and",  &
     neq," equations"
   WRITE(11,*)"Time after setup  is   :",elap_time()-timest(1)
 END IF
 ALLOCATE(loads_pp(neq_pp),diag_precon_pp(neq_pp),u_pp(neq_pp),           &
   d_pp(neq_pp),p_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp))
 storka_pp=.0_iwp
 storkb_pp=.0_iwp
 diag_precon_tmp=.0_iwp
 p_pp=.0_iwp
 diag_precon_pp=.0_iwp
 xnew_pp=.0_iwp
!----------- element integration ,storage and build preconditioner --------
 elements_2: DO iel=1,nels_pp
   num=g_num_pp(:,iel)
   g=g_g_pp(:,iel)
   coord=p_g_co_pp(:,:,iel)
   kc=0.0_iwp
   pm=0.0_iwp
   gauss_pts: DO i=1,nip
     CALL shape_der (der,points,i)
     CALL shape_fun(fun,points,i)
     funny(1,:)=fun(:)
     jac=MATMUL(der,coord)
     det=determinant(jac)
     CALL invert(jac)
     deriv=MATMUL(jac,der)
     kc=kc+MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
     pm=pm+MATMUL(TRANSPOSE(funny),funny)*det*weights(i) 
   END DO gauss_pts
   storka_pp(:,:,iel)=pm+kc*theta*dtim
   storkb_pp(:,:,iel)=pm-kc*(1._iwp-theta)*dtim
   DO k=1,ntot
     diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)+storka_pp(k,k,iel)
   END DO
 END DO elements_2
 CALL scatter(diag_precon_pp,diag_precon_tmp)
 DEALLOCATE(diag_precon_tmp)
 diag_precon_pp=1._iwp/diag_precon_pp 
!---------------------------initial conditions ---------------------------
 loads_pp=val0
 pmul_pp=.0_iwp
!----------------------time stepping recursion ---------------------------
 IF(numpe==it)THEN
   WRITE(11,'(A)')"    Time     Pressure  Iterations"
 END IF
 timesteps: DO j=1,nstep
   real_time=j*dtim
   u_pp=.0_iwp
   CALL gather(loads_pp,pmul_pp)
   elements_3: DO iel=1,nels_pp   
     utemp_pp(:,iel)=MATMUL(storkb_pp(:,:,iel),pmul_pp(:,iel))
   END DO elements_3
   CALL scatter(u_pp,utemp_pp)
   loads_pp=u_pp
!------------------- solve simultaneous equations by pcg -----------------
   d_pp=diag_precon_pp*loads_pp
   p_pp=d_pp
   x_pp=.0_iwp
   iters=0
   iterations: DO 
     iters=iters+1
     u_pp=0._iwp
     pmul_pp=.0_iwp
     CALL gather(p_pp,pmul_pp)  
     elements_4: DO iel=1,nels_pp
       utemp_pp(:,iel)=MATMUL(storka_pp(:,:,iel),pmul_pp(:,iel)) 
     END DO elements_4
     CALL scatter(u_pp,utemp_pp)
!--------------------------pcg equation solution--------------------------
     up=DOT_PRODUCT_P(loads_pp,d_pp)
     alpha= up/DOT_PRODUCT_P(p_pp,u_pp)
     xnew_pp=x_pp+p_pp*alpha
     loads_pp=loads_pp-u_pp*alpha
     d_pp=diag_precon_pp*loads_pp
     beta=DOT_PRODUCT_P(loads_pp,d_pp)/up
     p_pp=d_pp+p_pp*beta
     u_pp=xnew_pp
     CALL checon_par(xnew_pp,x_pp,tol,converged,neq_pp)
     IF(converged.OR.iters==limit)EXIT
   END DO iterations
   loads_pp=xnew_pp
   IF(j/npri*npri==j.AND.numpe==it)WRITE(11,'(2E12.4,I5)')real_time,      &
     loads_pp(is),iters
 END DO timesteps
 IF(numpe==it)WRITE(11,*)"This analysis took  :",elap_time()-timest(1)
 CALL shutdown() 
END PROGRAM p124
