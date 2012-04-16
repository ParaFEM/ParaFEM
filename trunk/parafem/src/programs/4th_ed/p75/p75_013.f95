!     Last change:  DV   19 Oct 2004    7:26 pm
PROGRAM p75
!-------------------------------------------------------------------------
! Program 7.5 General two- (plane) or three-dimensional analysis of steady
!             seepage. No global conductivity matrix assembly.
!             Diagonally preconditioned conjugate gradient solver.
!-------------------------------------------------------------------------
 USE main
 USE geom
 IMPLICIT NONE    
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER::cg_iters,cg_limit,fixed_freedoms,i,iel,j,k,loaded_nodes,nci,ndim, &
   nels,neq,nip,nod,nn,np_types,nlen
   !--------------- Added j
 REAL(iwp)::alpha,beta,cg_tol,det,one=1.0_iwp,penalty=1.0e20_iwp,up,        &
   zero=0.0_iwp
 CHARACTER(LEN=15)::element,argv
 LOGICAL::cg_converged
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g_num(:,:),node(:),num(:)
 REAL(iwp),ALLOCATABLE::coord(:,:),d(:),der(:,:),deriv(:,:),              &
   diag_precon(:),flux(:,:),g_coord(:,:),jac(:,:),kay(:,:),kc(:,:),        &
   loads(:),p(:),points(:,:),prop(:,:),store(:),storkc(:,:,:),u(:),       &
   value(:),weights(:),x(:),xnew(:),storeflux(:)
   !--------------- Added storeflux and changed flux dimensions
!-----------------------input and initialisation--------------------------
 CALL getname(argv,nlen)
 OPEN(10,FILE=argv(1:nlen)//'.dat')
 OPEN(11,FILE=argv(1:nlen)//'.res')
 READ(10,*)element,nod
 READ(10,*)nels,nn,nip,ndim,cg_tol,cg_limit
 READ(10,*)np_types
 neq=nn
 WRITE(11,'(A,I5,A)')" There are",neq," equations"
 ALLOCATE(points(nip,ndim),g_coord(ndim,nn),coord(nod,ndim),etype(nels),  &
   jac(ndim,ndim),weights(nip),num(nod),g_num(nod,nels),der(ndim,nod),    &
   deriv(ndim,nod),kc(nod,nod),kay(ndim,ndim),prop(ndim,np_types),        &
   p(0:neq),loads(0:neq),x(0:neq),xnew(0:neq),u(0:neq),diag_precon(0:neq),&
   d(0:neq),flux(0:nels*nip,ndim),storkc(nod,nod,nels),storeflux(ndim))
   !--------------- Added storeflux and changed flux(0:neq) for x, y, z and size of nels*nip
 READ(10,*)prop
 etype=1
 IF(np_types>1)READ(10,*)etype
 READ(10,*)g_coord
 READ(10,*)g_num
 IF(ndim==2)CALL mesh(g_coord,g_num,argv,nlen,12)
 diag_precon=zero
 CALL sample(element,points,weights)   
!----------element conductivity integration, storage and preconditioner--- 
 elements_1: DO iel=1,nels
   kay=zero
   DO i=1,ndim
     ! Give each element its property value (i.e. conductivity)
     kay(i,i)=prop(i,etype(iel))
   END DO
   num=g_num(:,iel)
   coord=TRANSPOSE(g_coord(:,num))
   kc=zero
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i)
     jac=MATMUL(der,coord)
     det=determinant(jac)
     CALL invert(jac)
     deriv=MATMUL(jac,der)
     kc=kc+MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
   END DO gauss_pts_1
   storkc(:,:,iel)=kc
   DO k=1,nod
     diag_precon(num(k))=diag_precon(num(k))+kc(k,k)
   END DO
 END DO elements_1
!-----------------------invert the preconditioner and get starting loads--
 loads=zero
 READ(10,*)loaded_nodes
 READ(10,*)(k,loads(k),i=1,loaded_nodes)
 READ(10,*)fixed_freedoms
 IF(fixed_freedoms/=0)THEN
   ALLOCATE(node(fixed_freedoms),value(fixed_freedoms),                   &
     store(fixed_freedoms))
   READ(10,*)(node(i),value(i),i=1,fixed_freedoms)
   diag_precon(node)=diag_precon(node)+penalty
   loads(node)=diag_precon(node)*value
   store=diag_precon(node)
 END IF
 diag_precon(1:)=one/diag_precon(1:)
 diag_precon(0)=zero
 d=diag_precon*loads
 p=d
!-----------------------pcg equation solution-----------------------------
 x=zero
 cg_iters=0
 pcg: DO
   cg_iters=cg_iters+1
   u=zero
   elements_2: DO iel=1,nels
     num=g_num(:,iel)
     kc=storkc(:,:,iel)
     u(num)=u(num)+MATMUL(kc,p(num)) 
   END DO elements_2
   IF(fixed_freedoms/=0)u(node)=p(node)*store
   up=DOT_PRODUCT(loads,d)
   alpha=up/DOT_PRODUCT(p,u)
   xnew=x+p*alpha
   loads=loads-u*alpha
   d=diag_precon*loads
   beta=DOT_PRODUCT(loads,d)/up
   p=d+p*beta
   CALL checon(xnew,x,cg_tol,cg_converged)
   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
 END DO pcg
 WRITE(11,'(A,I5)')" Number of cg iterations to convergence was",cg_iters
!-----------------------retrieve nodal net flow rates---------------------
 loads=xnew
 flux=zero
 storeflux=zero
 j=1
  
 elements_3: DO iel=1,nels
   
   !---Properties and node coords to be recalculated for each element
   kay=zero
   DO i=1,ndim
     kay(i,i)=prop(i,etype(iel))
   END DO
   num=g_num(:,iel)
   coord=TRANSPOSE(g_coord(:,num))    
 
   gauss_pts_2: DO i=1,nip
     !---Heat flux components calculated from equation HFL = -k*(dT/dx)
     !---dt/dx = derivative of shape function * nodal temp
     CALL shape_der(der,points,i)
     jac=MATMUL(der,coord)
     det=determinant(jac)
     CALL invert(jac)
     deriv=MATMUL(jac,der)
     
     storeflux=-MATMUL(kay,MATMUL(deriv,loads(num)))
   
     DO k=1,ndim
       flux(j,k)=storeflux(k)
     END DO     
     
     j=j+1
   END DO gauss_pts_2   
 END DO elements_3
 
 !----------------------- Write results to file ---------------------
 
 WRITE(11,'(/A)')"  Node Nodal Temp"
 DO k=1,nn
   WRITE(11,'(I5,1E12.4)')k,loads(k)
 END DO

 IF(ndim==2.AND.nod==4)THEN
   READ(10,*)nci
   CALL contour(loads,g_coord,g_num,nci,argv,nlen,13)
 END IF
 
 !---At each element Abaqus outputs the integration point heat fluxes
 !---in a specific order, therefore the calculated heat flux values are 
 !---reordered when written to file for direct comparison
 WRITE(11,'(/A)')"  IP   Flux_x      Flux_y      Flux_z"
 DO k=1,nels 
   WRITE(11,'(I5,3E12.4)')k,flux(8*k-0,1),flux(8*k-0,2),flux(8*k-0,3)
   WRITE(11,'(I5,3E12.4)')k,flux(8*k-2,1),flux(8*k-2,2),flux(8*k-2,3)
   WRITE(11,'(I5,3E12.4)')k,flux(8*k-4,1),flux(8*k-4,2),flux(8*k-4,3)
   WRITE(11,'(I5,3E12.4)')k,flux(8*k-5,1),flux(8*k-5,2),flux(8*k-5,3)
   WRITE(11,'(I5,3E12.4)')k,flux(8*k-1,1),flux(8*k-1,2),flux(8*k-1,3)
   WRITE(11,'(I5,3E12.4)')k,flux(8*k-3,1),flux(8*k-3,2),flux(8*k-3,3)
   WRITE(11,'(I5,3E12.4)')k,flux(8*k-6,1),flux(8*k-6,2),flux(8*k-6,3)
   WRITE(11,'(I5,3E12.4)')k,flux(8*k-7,1),flux(8*k-7,2),flux(8*k-7,3)
 END DO

 WRITE(11,'(/A)')"       HFL11       HFL22       HFL33"
 WRITE(11,'(5X,3E12.4)')                                    &
 SUM(flux(:,1),MASK=flux(:,1)>zero),SUM(flux(:,2),MASK=flux(:,2)>zero),SUM(flux(:,3),MASK=flux(:,3)>zero)
 WRITE(11,'(5X,3E12.4)')                                    &
 SUM(flux(:,1),MASK=flux(:,1)<zero),SUM(flux(:,2),MASK=flux(:,2)<zero),SUM(flux(:,3),MASK=flux(:,3)<zero)

STOP
END PROGRAM p75
