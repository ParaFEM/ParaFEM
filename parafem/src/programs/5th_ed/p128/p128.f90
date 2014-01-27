PROGRAM p128    
!-------------------------------------------------------------------------
!      Program 10.4 eigenvalues and eigenvectors of a cuboidal elastic
!      solid in 3d using  uniform 8-node hexahedral elements  
!      for lumped mass this is done element by element : parallel
!-------------------------------------------------------------------------
!USE mpi_wrapper  !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE new_library; USE eigen; IMPLICIT NONE
! neq,ntot are now global variables - not declared 
 INTEGER,PARAMETER::nodof=3,nod=8,nst=6,ndim=3
 INTEGER::nels,nn,nr,nip,i,j,k,l,iel,nmodes,jflag,lalfa,leig,lx,lz,      &
   iters,nlen,ndof,partitioner,meshgen,npes_pp,neig=0,iflag=-1,lp=11,    &
   node_end,node_start,nodes_pp
 REAL(iwp)::rho,e,v,det,el,er,acc ; REAL(iwp),PARAMETER::zero=0.0_iwp 
 CHARACTER(LEN=15)::element='hexahedron';  CHARACTER(LEN=50)::argv
 CHARACTER(LEN=6)::ch
!--------------------------- dynamic arrays-------------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),dee(:,:),vdiag_pp(:),eig(:),del(:),  &  
   fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),bee(:,:),val_pp(:),    &
   emm(:,:),ecm(:,:),utemp_pp(:,:),ua_pp(:),storkm_pp(:,:,:),timest(:),  &
   udiag_pp(:),diag_pp(:),alfa(:),beta(:),w1_pp(:),y_pp(:,:),z_pp(:,:),  &
   pmul_pp(:,:),v_store_pp(:,:),g_coord_pp(:,:,:),diag_tmp(:,:),x(:),    &
   va_pp(:),eigv_pp(:),temp(:)
 INTEGER,ALLOCATABLE::rest(:,:),g_num_pp(:,:),g_g_pp (:,:),nu(:),jeig(:,:)
!----------------------input and initialisation---------------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen) 
 CALL read_p128(argv,numpe,acc,e,el,er,lalfa,leig,lx,lz,meshgen,nels,    &
   nip,nmodes,nn,nr,partitioner,rho,v)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof
 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),           &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest); timest(2)=elap_time()
 ALLOCATE(points(nip,ndim),pmul_pp(ntot,nels_pp),fun(nod),dee(nst,nst),  &
   jac(ndim,ndim),weights(nip),der(ndim,nod),deriv(ndim,nod),x(lx),      &
   g_g_pp(ntot,nels_pp),ecm(ntot,ntot),eig(leig),del(lx),nu(lx),         &
   jeig(2,leig),alfa(lalfa),beta(lalfa),z_pp(lz,leig),bee(nst,ntot),     &
   utemp_pp(ntot,nels_pp),emm(ntot,ntot),diag_tmp(ntot,nels_pp),         &
   storkm_pp(ntot,ntot,nels_pp))
!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_0: DO iel=1,nels_pp
   CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_0
 neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
 IF(numpe==1)THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A,I7,A)') "This job ran on ",npes," processes"
   WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes", nr, &
                           " restrained and ",neq," equations"
   WRITE(11,'(A,F10.4)') "Time to read input is:",timest(2)-timest(1)
   WRITE(11,'(A,F10.4)') "Time after setup is:",elap_time()-timest(1)
 END IF
 ALLOCATE(ua_pp(neq_pp),va_pp(neq_pp),vdiag_pp(neq_pp),udiag_pp(neq_pp), &
   v_store_pp(neq_pp,lalfa),diag_pp(neq_pp),w1_pp(neq_pp),               &
   y_pp(neq_pp,leig)); ua_pp=zero; va_pp=zero; eig=zero
 diag_tmp=zero; jeig=0; x=zero; del=zero; nu=0; alfa=zero; beta=zero
 diag_pp=zero; udiag_pp=zero; w1_pp=zero; y_pp=zero; z_pp=zero
 CALL sample(element,points,weights); CALL deemat(dee,e,v)
!--------------- element stiffness integration and assembly---------------
 elements_2: DO iel=1,nels_pp; emm=zero
   integrating_pts_1: DO i=1,nip
     CALL shape_fun(fun,points,i); CALL shape_der(der,points,i)
     jac=MATMUL(der,g_coord_pp(:,:,iel)); det=determinant(jac)
     CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
     storkm_pp(:,:,iel)=storkm_pp(:,:,iel)+                              &
                   MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
     CALL ecmat(ecm,fun,ntot,nodof); emm=emm+ecm*det*weights(i)*rho
   END DO integrating_pts_1
   DO k=1,ntot; diag_tmp(k,iel)=diag_tmp(k,iel)+sum(emm(k,:)); END DO
 END DO elements_2
 CALL scatter(diag_pp,diag_tmp); DEALLOCATE(diag_tmp)
!---------------------------- find eigenvalues ---------------------------
 diag_pp=1._iwp/sqrt(diag_pp) ! diag_pp holds l**(-1/2)
 DO iters=1,lalfa
   CALL lancz1(neq_pp,el,er,acc,leig,lx,lalfa,lp,iflag,ua_pp,va_pp,eig,  &
     jeig,neig,x,del,nu,alfa,beta,v_store_pp)
   IF(iflag==0)EXIT
   IF(iflag>1)THEN; 
     IF(numpe==1)THEN        
       WRITE(11,'(A,I5)')                                                &
         " Lancz1 is signalling failure, with iflag = ",iflag; EXIT
     END IF 
   END IF           
!----- iflag = 1 therefore form u + a * v  (done element by element) -----
   vdiag_pp=va_pp; vdiag_pp=vdiag_pp*diag_pp  ! vdiag is l**(-1/2).va
   udiag_pp=zero; pmul_pp=zero
   CALL gather(vdiag_pp,pmul_pp)
   elements_3: DO iel=1,nels_pp
     utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
   END DO elements_3
   CALL scatter(udiag_pp,utemp_pp)   ! udiag is A.l**(-1/2).va
   udiag_pp=udiag_pp*diag_pp; ua_pp=ua_pp+udiag_pp
 END DO
!-------------- iflag = 0 therefore write out the spectrum --------------- 
 IF(numpe==1)THEN 
   WRITE(11,'(2(A,E12.4))')"The range is",el,"  to ",er
   WRITE(11,'(A,I8,A)')"There are ",neig," eigenvalues in the range"
   WRITE(11,'(A,I8,A)')"It took ",iters,"  iterations"
   WRITE(11,'(A)')"The eigenvalues are:"; WRITE(11,'(6E12.4)')eig(1:neig)  
 END IF     
!  calculate the eigenvectors
 IF(neig>10)neig=10
 CALL lancz2(neq_pp,lalfa,lp,eig,jeig,neig,alfa,beta,lz,jflag,y_pp,      &
   w1_pp,z_pp,v_store_pp)
!-------------- if jflag is zero calculate the eigenvectors --------------
 IF(jflag==0)THEN; timest(3)=elap_time()
   CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
   ALLOCATE(eigv_pp(nodes_pp*ndim),temp(nodes_pp))
   eigv_pp=zero; utemp_pp=zero; temp=0   
   IF(numpe==1) WRITE(11,'(A)')"The eigenvectors are  :" 
   DO i=1,nmodes
     udiag_pp(:)=y_pp(:,i); udiag_pp=udiag_pp*diag_pp
     IF(numpe==1) THEN
       WRITE(11,'("Eigenvector number  ",I4," is: ")')i
       WRITE(11,'(6E12.4)')udiag_pp(1:6)
       WRITE(ch,'(I6.6)') i
       OPEN(12,file=argv(1:nlen)//".ensi.EIGV-"//ch,status='replace',    &
         action='write'); WRITE(12,'(A)')                                &
         "Alya Ensight Gold --- Vector per-node variable file"
       WRITE(12,'(A/A/A)') "part", "     1","coordinates"
     END IF
     CALL gather(udiag_pp(1:),utemp_pp); eigv_pp=zero
     CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,      &
                        node_start,node_end,utemp_pp,eigv_pp,1)
     DO j=1,ndim ; temp=zero
       DO l=1,nodes_pp; k=j+(ndim*(l-1)); temp(l)=eigv_pp(k); END DO
         CALL dismsh_ensi_p(12,l,nodes_pp,npes,numpe,1,temp)
     END DO ; IF(numpe==1) CLOSE(12)
   END DO
   IF(numpe==1) THEN
     WRITE(11,'(A,F10.4)')"Writing out eigenvectors took:",              &
       elap_time()-timest(3)
   END IF
 ELSE
   IF(numpe==1) THEN  ! lancz2 fails
     WRITE(11,'(A,I5)')"Lancz2 is signalling failure with jflag= ", jflag
   END IF
 END IF
 IF(numpe==1)WRITE(11,'(A,F10.4)')                                       &
   "This analysis took:",elap_time()-timest(1)
 CALL SHUTDOWN()  
END PROGRAM p128
