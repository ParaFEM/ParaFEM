PROGRAM p128    
!-------------------------------------------------------------------------
!      Program 10.4 eigenvalues and eigenvectors of a cuboidal elastic
!      solid in 3d using  uniform 8-node hexahedral elements  
!      for lumped mass this is done element by element : parallel version
!-------------------------------------------------------------------------
 USE new_library
 USE geometry_lib
 USE lancz_lib
 USE precision
 USE timing
 USE utility
 USE mp_module
 USE global_variables1
 USE gather_scatter6
 IMPLICIT NONE
!------------ ndof,nels,neq,ntot are global - not declared----------------
 INTEGER::nxe,nye,nze,nn,nr,nip,nodof=3,nod=8,nst=6,i,j,k,iel,ndim=3,     &
   nmodes,jflag,iflag=-1,lp=11,lalfa,leig,lx,lz,iters,neig=0,nlen
 REAL(iwp)::aa,bb,cc,rho,e,v,det,el,er,acc  
 CHARACTER(LEN=15)::element='hexahedron',argv
!--------------------------- dynamic arrays-------------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),dee(:,:),coord(:,:),vdiag_pp(:),      &  
   fun(:),jac(:,:),der(:,:),deriv(:,:),weights(:),bee(:,:),km(:,:),       &
   emm(:,:),ecm(:,:),utemp_pp(:,:),ua_pp(:),va_pp(:),eig(:),del(:),       &
   udiag_pp(:),diag_pp(:),alfa(:),beta(:),w1_pp(:),y_pp(:,:),z_pp(:,:),   &
   pmul_pp(:,:),v_store_pp(:,:),p_g_co_pp(:,:,:),diag_tmp(:,:),x(:)
 INTEGER,ALLOCATABLE::rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp (:,:),   &
   nu(:),jeig(:,:)
!----------------------input and initialisation---------------------------
 timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes)
 IF(numpe==npes)THEN
   CALL getname(argv,nlen)
   OPEN(10,FILE=argv(1:nlen)//'.dat',STATUS='OLD',ACTION='READ')
   OPEN(11,FILE=argv(1:nlen)//'.res',STATUS='REPLACE',ACTION='WRITE')              
   READ(10,*)nels,nxe,nze,nip,aa,bb,cc,rho,e,v,nmodes,el,er,lalfa,leig,lx,&
     lz,acc
 END IF
 CALL bcast_inputdata_p128(numpe,npes,nels,nxe,nze,nip,aa,bb,cc,rho,e,v,  &
   nmodes,el,er,lalfa,leig,lx,lz,acc)
 CALL calc_nels_pp
 ndof=nod*nodof
 ntot=ndof
 nn=0
 neq=0
 nr=(nxe+1)*(nze+1)
 nye=nels/nxe/nze
 ALLOCATE(rest(nr,nodof+1),points(nip,ndim),pmul_pp(ntot,nels_pp),        &
   coord(nod,ndim),fun(nod),jac(ndim,ndim),weights(nip),                  &
   g_num_pp(nod,nels_pp),der(ndim,nod),deriv(ndim,nod),dee(nst,nst),      &
   num(nod),km(ntot,ntot),g(ntot),g_g_pp(ntot,nels_pp),ecm(ntot,ntot),    &
   eig(leig),x(lx),del(lx),nu(lx),jeig(2,leig),alfa(lalfa),beta(lalfa),   &
   z_pp(lz,leig),utemp_pp(ntot,nels_pp),p_g_co_pp(nod,ndim,nels_pp),      &
   bee(nst,ntot),emm(ntot,ntot),diag_tmp(ntot,nels_pp))
 rest=0
 DO i=1,nr
   rest(i,1)=i
 END DO
 ielpe=iel_start
!------------------ loop the elements to set up global arrays  -----------
 elements_1: DO iel=1,nels_pp
   CALL geometry_8bxz(ielpe,nxe,nze,aa,bb,cc,coord,num)
   CALL find_g(num,g,rest)
   g_num_pp(:,iel)=num
   p_g_co_pp(:,:,iel)=coord
   g_g_pp(:,iel)=g
   ielpe=ielpe+1
 END DO elements_1
 nn=(nxe+1)*(nze+1)*(nye+1)
 neq=(nn-nr)*nodof
 CALL calc_neq_pp
 CALL make_ggl(g_g_pp)
 IF(numpe==npes)THEN
   WRITE(11,'(A,I5,A)')"This job ran on ",npes,"   processors"
   WRITE(11,'(A)')"Global coordinates and node numbers"
   DO i=1,nels_pp,nels_pp-1
     WRITE(11,'(A,I8)')"Element ",i
     num=g_num_pp(:,i)
     DO k=1,nod
       WRITE(11,'(A,I8,3E12.4)')"  Node",num(k),p_g_co_pp(k,:,i)
     END DO
   END DO
   WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nr,                    &
     " restrained  and",neq," equations"
   WRITE(11,*)"Time after setup is  :",elap_time()-timest(1)
 END IF
 ALLOCATE(ua_pp(neq_pp),va_pp(neq_pp),vdiag_pp(neq_pp),                   &
   v_store_pp(neq_pp,lalfa),diag_pp(neq_pp),udiag_pp(neq_pp),             &
   w1_pp(neq_pp),y_pp(neq_pp,leig))
 ua_pp=.0_iwp
 va_pp=.0_iwp
 eig=.0_iwp
 diag_tmp=.0_iwp
 jeig=0
 x=.0_iwp
 del=.0_iwp
 nu=0
 alfa=.0_iwp
 beta=.0_iwp
 diag_pp=.0_iwp
 udiag_pp=.0_iwp
 w1_pp=.0_iwp
 y_pp=.0_iwp
 z_pp=.0_iwp
 CALL sample(element,points,weights)
 CALL deemat(dee,e,v)
!--------------- element stiffness integration and assembly---------------
 elements_2: DO iel=1,nels_pp
   coord=p_g_co_pp(:,:,iel)
   g=g_g_pp(:,iel)
   km=0.0_iwp
   emm=0.0_iwp
   integrating_pts_1: DO i=1,nip
     CALL shape_fun(fun,points,i)
     CALL shape_der(der,points,i)
     jac=MATMUL(der,coord)
     det=determinant(jac)
     CALL invert(jac)
     deriv=MATMUL(jac,der)
     CALL beemat(bee,deriv)
     km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
     CALL ecmat(ecm,fun,ntot,nodof)
     emm=emm+ecm*det*weights(i)*rho
   END DO integrating_pts_1
   DO k=1,ntot
     diag_tmp(k,iel)=diag_tmp(k,iel)+sum(emm(k,:))
   END DO
 END DO elements_2
 CALL scatter(diag_pp,diag_tmp)
 DEALLOCATE(diag_tmp)
!------------------------------find eigenvalues---------------------------
 diag_pp=1._iwp/sqrt(diag_pp) ! diag_pp holds l**(-1/2)
 DO iters=1,lalfa
   CALL lancz1(neq_pp,el,er,acc,leig,lx,lalfa,lp,iflag,ua_pp,va_pp,eig,   &
     jeig,neig,x,del,nu,alfa,beta,v_store_pp)
   IF(iflag==0)EXIT
   IF(iflag>1)THEN  
     IF(numpe==npes)THEN        
       WRITE(11,'(A,I5)')                                                 &
         " Lancz1 is signalling failure, with iflag = ",iflag
       EXIT
     END IF 
   END IF           
!---- iflag = 1 therefore form u + a * v  ( done element by element )-----
   vdiag_pp=va_pp
   vdiag_pp=vdiag_pp*diag_pp  ! vdiag is l**(-1/2).va
   udiag_pp=.0_iwp
   pmul_pp=.0_iwp
   CALL gather(vdiag_pp,pmul_pp)
   elements_3: DO iel=1,nels_pp
!    utemp_pp(:,iel) = MATMUL(km,pmul_pp(:,iel))
     CALL dgemv('n',ntot,ntot,1.0,km,ntot,pmul_pp(:,iel),1,0.0,           &
       utemp_pp(:,iel),1)
   END DO elements_3
   CALL scatter(udiag_pp,utemp_pp)   ! udiag is A.l**(-1/2).va
   udiag_pp=udiag_pp*diag_pp
   ua_pp=ua_pp+udiag_pp
 END DO
!-------------- iflag = 0 therefore write out the spectrum --------------- 
 IF(numpe==npes)THEN 
   WRITE(11,'(2(A,E12.4))')"The range is",el,"  to ",er
   WRITE(11,'(A,I8,A)')"There are ",neig," eigenvalues in the range"
   WRITE(11,'(A,I8,A)')"It took ",iters,"  iterations"
   WRITE(11,'(A)')"The eigenvalues are   :"
   WRITE(11,'(6E12.4)')eig(1:neig)  
 END IF     
!  calculate the eigenvectors
 IF(neig>10)neig=10
 CALL lancz2(neq_pp,lalfa,lp,eig,jeig,neig,alfa,beta,lz,jflag,y_pp,w1_pp, &
   z_pp,v_store_pp)
!------------------if jflag is zero  calculate the eigenvectors ----------
 IF(jflag==0)THEN   
   IF(numpe==npes)THEN 
     WRITE(11,'(A)')"The eigenvectors are  :"  
     DO i=1,nmodes
       udiag_pp(:)=y_pp(:,i)
       udiag_pp=udiag_pp*diag_pp
       WRITE(11,'("Eigenvector number  ",I4," is: ")')i
       WRITE(11,'(6E12.4)')udiag_pp(1:6)
     END DO
   ELSE
! lancz2 fails
     WRITE(11,'(A,I5)')" Lancz2 is signalling failure with jflag = ",jflag  
   END IF 
 END IF
 IF(numpe==npes)WRITE(11,*)"This analysis took  :",elap_time()-timest(1)
 CALL shutdown()  
END PROGRAM p128
