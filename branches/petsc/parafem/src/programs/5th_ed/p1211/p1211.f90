PROGRAM p1211         
!-------------------------------------------------------------------------
! Program 12.11 Three-dimensional strain of an elastic solid using
!               8-, 14- or 20-node brick hexahedra. Mesh numbered in x-y
!               planes then in the z-direction. No global stiffness matrix
!               assembly. Diagonally preconditioned conjugate gradient 
!               solver. GPU version.
!-------------------------------------------------------------------------
 USE new_library; USE geometry; USE precision; USE timing
 USE maths; USE global_variables; USE input; IMPLICIT NONE
!INCLUDE 'altcublas.inc' ! Include HMPPALT cublas proxy interface
!neq is global, must not be declared 
 INTEGER::err ! For HMPPALT error management
 REAL(iwp)::coeffa=1_iwp, coeffb=0_iwp ! Coefficents for DGEMM 
 INTEGER::cg_iters,cg_limit,fixed_freedoms,i,iel,k,loaded_nodes,ndim=3,   &
   ndof,nels,nip,nn,nprops=2,np_types,nod,nodof=3,nr,nst=6,nxe,nye,nze,nlen
 REAL(iwp)::alpha,beta,big,cg_tol,det,one=1.0_iwp,penalty=1.0e20_iwp,up,  &
   zero=0.0_iwp
 CHARACTER(LEN=15)::element='hexahedron',argv
 LOGICAL::cg_converged
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),no(:),    &
   node(:),num(:),sense(:)
 REAL(iwp),ALLOCATABLE::bee(:,:),coord(:,:),d(:),dee(:,:),der(:,:),       &
   deriv(:,:),diag_precon(:),eld(:),fun(:),gc(:),g_coord(:,:),g_pmul(:,:),&
   g_utemp(:,:),jac(:,:),km(:,:),loads(:),p(:),points(:,:),prop(:,:),     &
   sigma(:),store(:),storkm(:,:,:),u(:),value(:),weights(:),x(:),xnew(:), &
   x_coords(:),y_coords(:),z_coords(:),timest(:),oldlds(:)
 !-----------------------external subroutines------------------------------
 EXTERNAL::DDOT ; REAL(iwp)::DDOT ! For MKL
!-----------------------input and initialisation--------------------------
 ALLOCATE(timest(2)); timest=zero; timest(1) = elap_time()
 CALL getname(argv,nlen)
 OPEN(10,FILE=argv(1:nlen)//'.dat'); OPEN(11,FILE=argv(1:nlen)//'.res')
 READ(10,*)nod,nxe,nye,nze,nip,cg_tol,cg_limit,np_types
 CALL mesh_size(element,nod,nels,nn,nxe,nye,nze); ndof=nod*nodof
 ALLOCATE(nf(nodof,nn),points(nip,ndim),dee(nst,nst),coord(nod,ndim),     &
   jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),fun(nod),gc(ndim),        &
   bee(nst,ndof),km(ndof,ndof),eld(ndof),sigma(nst),g_coord(ndim,nn),     &
   g_num(nod,nels),weights(nip),num(nod),g_g(ndof,nels),x_coords(nxe+1),  &
   g(ndof),y_coords(nye+1),z_coords(nze+1),etype(nels),g_pmul(ndof,nels), &
   prop(nprops,np_types),storkm(ndof,ndof,nels),g_utemp(ndof,nels))
 READ(10,*)prop; etype=1; IF(np_types>1)READ(10,*)etype
 READ(10,*)x_coords,y_coords,z_coords
 nf=1;  READ(10,*)nr,(k,nf(:,k),i=1,nr); CALL formnf(nf); neq=MAXVAL(nf)
 ALLOCATE(p(0:neq),loads(0:neq),x(0:neq),xnew(0:neq),u(0:neq),            &
   diag_precon(0:neq),d(0:neq),oldlds(0:neq))
 diag_precon=zero; CALL sample(element,points,weights)
 timest(2)=elap_time()
!----------element stiffness integration, storage and preconditioner------
 elements_1: DO iel=1,nels
   CALL hexahedron_xz(iel,x_coords,y_coords,z_coords,coord,num)
   CALL num_to_g(num,nf,g); g_num(:,iel)=num
   g_coord(:,num)=TRANSPOSE(coord); g_g(:,iel)=g
   CALL deemat(dee,prop(1,etype(iel)),prop(2,etype(iel)));num=g_num(:,iel)
   coord=TRANSPOSE(g_coord(:,num)); g=g_g(:,iel); km=zero
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i); jac=MATMUL(der,coord)
     det=determinant(jac); CALL invert(jac)
     deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
     km=km+MATMUL(matmul(transpose(bee),dee),bee)*det*weights(i)
   END DO gauss_pts_1
   storkm(:,:,iel) = km
   DO k=1,ndof; diag_precon(g(k))=diag_precon(g(k))+km(k,k); END DO
 END DO elements_1  
!-----------------------invert the preconditioner and get starting loads--
 loads=zero; READ(10,*)loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)
 oldlds=loads; READ(10,*)fixed_freedoms
 IF(fixed_freedoms/=0)THEN
   ALLOCATE(node(fixed_freedoms),sense(fixed_freedoms),                   &
     value(fixed_freedoms),no(fixed_freedoms),store(fixed_freedoms))
   READ(10,*)(node(i),sense(i),value(i),i=1,fixed_freedoms)
   DO  i=1,fixed_freedoms; no(i)=nf(sense(i),node(i)); END DO
   diag_precon(no)=diag_precon(no)+penalty;loads(no)=diag_precon(no)*value
   store=diag_precon(no)
 END IF
 diag_precon(1:)=one/diag_precon(1:); diag_precon(0)=zero
 d=diag_precon*loads; p=d
!-----------------------pcg equation solution-----------------------------
 !$hmpp htspt1 acquire
 !$hmpp htspt1 region, target=CUDA
 !$hmpp htspt1 endregion
 x=zero; cg_iters=0
 !$hmpp htspt1 allocate data["km"], size={ndof,ndof}
 !$hmpp htspt1 allocate data["g_pmul"], size={ndof,nels}
 !$hmpp htspt1 allocate data["g_utemp"], size={ndof,nels}
 !$hmpp htspt1 advancedload data["km"], asynchronous
 pcg: DO
   cg_iters=cg_iters+1; u=zero; g_utemp=zero
   elements_2: DO iel=1,nels
     g_pmul(:,iel)=p(g_g(:,iel))
   END DO elements_2
   !$hmpp htspt1 advancedload data["g_pmul"]
   !$hmpp htspt1 waitload data["km"]
   !$hmppalt cublas call, name="dgemm", error="err"
   CALL DGEMM('N','N',ndof,nels,ndof,coeffa,km,ndof,g_pmul,ndof,coeffb,  &
     g_utemp,ndof)
   !$hmpp htspt1 delegatedstore data["g_utemp"]
   elements_2a: DO iel=1,nels
     u(g_g(:,iel)) = u(g_g(:,iel))+g_utemp(:,iel)
   END DO elements_2a
   IF(fixed_freedoms/=0) u(no)=p(no)*store
   up=DDOT(neq,loads,1,d,1); alpha=up/DDOT(neq,p,1,u,1)
   xnew=x+p*alpha; loads=loads-u*alpha; d=diag_precon*loads
   beta=DDOT(neq,loads,1,d,1)/up; p=d+p*beta
   CALL checon(xnew(1:),x(1:),cg_tol,cg_converged)
   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
 END DO pcg
 !$hmpp htspt1 release
 WRITE(11,'(A)')                                                         &
  "This analysis ran on 1 process and was accelerated using 1 GPU"
 WRITE(11,'(A,I8,A)') "There are",neq," equations"
 WRITE(11,'(A,I5)') "Number of iterations to convergence was ",cg_iters
 WRITE(11,'(/A)')"       Node   x-disp      y-disp      z-disp"
 DO k=1,nn,nn-1; WRITE(11,'(I10,3E12.4)')k,xnew(nf(:,k)); END DO
!----------------- recover stresses at nip integrating point ---------------
 nip       = 1; DEALLOCATE(points,weights)
 ALLOCATE(points(nip,ndim),weights(nip))
 CALL sample(element,points,weights)
 WRITE(11,'(/A,I2,A)')"The integration point (nip=",nip,") stresses are:"
 WRITE(11,'(A,/,A)')"    Element     x-coord     y-coord     z-coord",    &
  "    sig_x       sig_y       sig_z       tau_xy      tau_yz      tau_zx" 
 xnew(0)=zero
 elements_3: DO iel=1,nels
   CALL deemat(dee,prop(1,etype(iel)),prop(2,etype(iel)))
   num=g_num(:,iel); coord=TRANSPOSE(g_coord(:,num)); g=g_g(:,iel)
   eld=xnew(g)
   gauss_pts_2: DO i=1,nip
     CALL shape_der(der,points,i); CALL shape_fun(fun,points,i)
     gc=MATMUL(fun,coord); jac=MATMUL(der,coord)
     CALL invert(jac); deriv=MATMUL(jac,der)
     CALL beemat(bee,deriv); sigma=MATMUL(dee,MATMUL(bee,eld))
     IF(iel<4) THEN 
       WRITE(11,'(I8,4X,3E12.4)')iel,gc; WRITE(11,'(6E12.4)')sigma
     END IF
   END DO gauss_pts_2
 END DO elements_3
 WRITE(11,'(/A,F12.4,A)') "Time for setup was                        ",   &
                           timest(2)-timest(1),"s"
 WRITE(11,'(A,F12.4,A)') "Total time for this analysis was          ",    &
                           elap_time()-timest(1),"s"
 STOP
END PROGRAM p1211
