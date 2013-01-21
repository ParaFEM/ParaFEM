PROGRAM P128
!-------------------------------------------------------------------------
!      program p12.8 eigenvalues and eigenvectors of a cuboidal elastic
!      solid in 3d using uniform 8-node hexahedral elements ; Arnoldi "dp" 
!      for lumped mass this is done element by element : parallel version
!-------------------------------------------------------------------------
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE maths; USE gather_scatter; USE partition
 USE elements; USE steering; USE pcg !; USE timing
 IMPLICIT NONE
! neq,ntot are global variables - not declared
 INTEGER::nn,nr,nip,nodof=3,nod,nst=6,i,j,jj,k,l,inode,iel,ndim=3,iters, &
   model,nconv,ncv,nev,maxitr,argc,iargc,partitioner=1,nels,ndof,        &
   npes_pp,meshgen,first,last,ic,ido,ierr,info,iparam(11),ipntr(11),     &
   ishfts,lworkl,node_end,node_start,nodes_pp,nlen  
 REAL(iwp)::rho,e,nu,det,sigma,tol
 REAL(iwp),PARAMETER::zero=0.0_iwp, one=1.0_iwp  
 CHARACTER (LEN=15)::element; CHARACTER(LEN=50)::argv 
 CHARACTER(LEN=6)::ch; CHARACTER::bmat*1, which*2 
 LOGICAL::rvec
!---------------------------- dynamic arrays -----------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),dee(:,:),coord(:,:),vdiag(:),        &
   resid(:),fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),bee(:,:),     &
   store_km_pp(:,:,:),emm(:,:),ecm(:,:),diag(:),pmul(:),d(:,:),v(:,:),   &
   workl(:),workd(:),g_coord_pp(:,:,:),diag1(:),udiag(:),utemp(:),       &
   eigv(:,:),eigv1(:,:),timest(:),temp(:)
 INTEGER, ALLOCATABLE::rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:),  &
   g_g_local(:,:),g_g(:,:),g_num(:,:),g_num_local(:,:)
 LOGICAL, ALLOCATABLE::select(:)   
 INCLUDE  'debug.h'     ! From ARPACK library
!------------------------ input and initialisation -----------------------
 ALLOCATE(timest(20)); CALL cpu_time(timest(1))
!timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen) 
 CALL read_xx6(argv,numpe,bmat,e,element,maxitr,meshgen,ncv,nels,    &
   nev,nip,nn,nod,nr,rho,tol,nu,which)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof
 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),           &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp2(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest)! ; timest(2)=elap_time()
!----------- initialise values of parameters required by ARPACK ----------
 lworkl=ncv*(ncv+8); ndigit=-3; logfil=6; msgets=0; mseupd=0; msaitr=0
 msaupd=1; msaup2=0; mseigt=0; msapps=0; info=0; ishfts=1; model=1; ido=0    
 iparam(1)=ishfts; iparam(3)=maxitr; iparam(7)=model
!-------------- allocate dynamic arrays used in main program -------------
 ALLOCATE(points(nip,ndim),pmul(ntot),coord(nod,ndim),fun(nod),          &
   jac(ndim,ndim), weights(nip),der(ndim,nod),deriv(ndim,nod),           &
   dee(nst,nst),num(nod),store_km_pp(ntot,ntot,nels_pp),                 &
   g_g_pp(ntot,nels_pp),ecm(ntot,ntot),utemp(ntot),d(ncv,2),             &
   workl(lworkl),select(ncv),bee(nst,ntot),emm(ntot,ntot),g(ntot))    
!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_1: DO iel=1,nels_pp
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_1
 elements_2: DO iel=1,nels_pp  
    i=MAXVAL(g_g_pp(:,iel)); IF(i>neq) neq=i
 END DO elements_2  
 neq=MAX_INTEGER_P(neq); CALL calc_neq_pp; CALL calc_npes_pp(npes,npes_pp)
 CALL make_ggl2(npes_pp,npes,g_g_pp)  
 ALLOCATE(v(neq,ncv),diag(neq),udiag(neq),vdiag(neq),workd(3*neq),       &
   resid(neq),diag1(neq)); diag1=zero; v=zero; diag=zero; udiag=zero
   vdiag=zero; workd=zero; resid=zero; diag1=zero
!------ element stiffness integration and build the preconditioner -------
 CALL sample(element,points,weights); CALL deemat(e,nu,dee)
 store_km_pp=zero
 elements_3: DO iel=1,nels_pp
   emm=zero
   integrating_pts_1:  DO i=1,nip
     CALL shape_fun(fun,points,i); CALL shape_der(der,points,i)
     jac=MATMUL(der,g_coord_pp(:,:,iel)); det=determinant(jac) 
     CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(deriv,bee)
     store_km_pp(:,:,iel)=store_km_pp(:,:,iel)+                          & 
                          MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*        &
                          det*weights(i)
     CALL ecmat(ecm,fun,ntot,nodof); emm=emm+ecm*det*weights(i)*rho
   END DO integrating_pts_1
   DO k=1,ntot; IF(g_g_pp(k,iel)/=0) THEN
     diag1(g_g_pp(k,iel))=diag1(g_g_pp(k,iel))+SUM(emm(k,:)); END IF        
   END DO
 END DO elements_3
 CALL MPI_ALLREDUCE(diag1,diag,neq,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)   
!-------------------------- find the eigenvalues -------------------------
 iters=0; diag=one/sqrt(diag) ! diag holds l**(-1/2) 
 DO
   iters=iters+1
   CALL dsaupd(ido,bmat,neq,which,nev,tol,resid,                 &
               ncv,v,neq,iparam,ipntr,workd,workl,lworkl,info) 
   IF(ido/=-1 .AND. ido/=1) EXIT
   diag1=zero; vdiag=workd(ipntr(1):ipntr(1)+neq-1)*diag
   elements_4: DO iel=1,nels_pp
     DO i=1,ntot
       IF(g_g_pp(i,iel)==0) pmul(i)=zero
       IF(g_g_pp(i,iel)/=0) pmul(i)=vdiag(g_g_pp(i,iel))
     END DO
     utemp=MATMUL(store_km_pp(:,:,iel),pmul)   
!    CALL dgemv('n',ntot,ntot,one,km,ntot,pmul,1,zero,utemp,1) 
     DO i=1,ntot
       IF(g_g_pp(i,iel)/=0) THEN
         diag1(g_g_pp(i,iel))=diag1(g_g_pp(i,iel))+utemp(i)
       END IF
     END DO
   END DO elements_4
   CALL MPI_ALLREDUCE(diag1,udiag,neq,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,  &
     ier); udiag=udiag*diag; workd(ipntr(2):ipntr(2)+neq-1)=udiag
 END DO
 IF(numpe==1) THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A,I8,A)') "It took ",iters,"  iterations"
 END IF
 IF(info < 0) THEN ! Either we have convergence or there is an error
   IF(numpe==1) WRITE(11,*) "Fatal error in dsaupd "
 ELSE
   rvec=.TRUE.
   CALL dseupd(rvec,'All',select,d,v,neq,sigma,bmat,neq,which,nev,tol,   &
     resid,ncv,v,neq,iparam,ipntr,workd,workl,lworkl,ierr)
   IF(ierr/=0) THEN; IF(numpe==1) THEN
     WRITE(11,*) "Fatal error in dseupd "
     WRITE(11,*) "Should compute residuals "
   END IF; END IF
   IF(numpe==1) THEN
     WRITE(11,'(A)') "The eigenvalues are  :"
     WRITE(11,'(5E12.4)') d(1:nev,1)
     WRITE(11,'(A)') "The eigenvectors are  :"
     DO i=1,nev
        udiag(:)=v(1:neq,i); udiag=udiag*diag  
        WRITE(11,'("Eigenvector number  ",I4," is: ")')  i
        WRITE(11,'(6E12.4)') udiag(1:6)
     END DO
   END IF
 END IF
!---------------------- write out the eigenvectors -----------------------   
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 ALLOCATE(eigv(ndim,nn),g_num(nod,nels),g_num_local(nod,nels),           &
   g_g(ntot,nels),g_g_local(ntot,nels))
 eigv=zero; g_num=0; g_num_local=0; g_g=0; g_g_local=0
 first=iel_start; last=iel_start+nels_pp-1
 g_num_local(:,first:last)=g_num_pp(:,:)
 g_g_local(:,first:last)=g_g_pp(:,:); ic=nod*nels 
 CALL MPI_ALLREDUCE(g_num_local,g_num,ic,MPI_INTEGER,MPI_SUM,            &
   MPI_COMM_WORLD,ier)
 ic=ntot*nels
 CALL MPI_ALLREDUCE(g_g_local,g_g,ic,MPI_INTEGER,MPI_SUM,                &
   MPI_COMM_WORLD,ier)
 DO i=1,nev
   IF(numpe==1) THEN; WRITE(ch,'(I6.6)') i
     OPEN(12,file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',       &
       action='write')
     WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
     WRITE(12,'(A/A/A)') "part", "     1","coordinates"
   END IF
   eigv=zero; udiag=zero; udiag(:)=v(1:neq,i); udiag=udiag*diag  
   DO iel=1,nels
     DO j=1,nod; inode=g_num(j,iel)
       DO k=1,ndim; l=((j-1)*ndim)+k
         IF(g_g(l,iel)/=0) eigv(k,inode)=udiag(g_g(l,iel))
   END DO; END DO; END DO
   IF(numpe==1) THEN
     DO jj=1,nn; k=i+(ndim*(jj-1)); temp(jj)=eigv(k,j); END DO
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
   END IF
   IF(numpe==1) CLOSE(12) ! *** need to open and close ***
 END DO
!IF(numpe==1) WRITE(11,*) "This analysis took  :", elap_time( ) - timest(1)
 CALL shutdown()
 END PROGRAM P128
