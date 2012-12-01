PROGRAM p125      
!-------------------------------------------------------------------------
!   Program 8.4 conduction equation on a 3-d box volume using 8-node
!   hexahedral elements and a simple explicit algorithm : parallel version
!   box_bc : write on processor it at freedom nres
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
!---- ndof, nels, neq , ntot  are now global variables - not declared-----
 INTEGER::nxe,nye,nze,nn,nr,nip,nodof=1,nod=8,ndim=3,i,j,k,iel,neq_temp,  &
   nn_temp,nstep,npri,nres,it,is,nlen
 REAL(iwp)::aa,bb,cc,kx,ky,kz,det,dtim,val0,real_time
 CHARACTER(LEN=15)::element='hexahedron',argv
!------------------------- dynamic arrays---------------------------------
 REAL(iwp),ALLOCATABLE::loads_pp(:),points(:,:),kay(:,:),coord(:,:),      &
   jac(:,:),der(:,:),deriv(:,:),weights(:),kc(:,:),pm(:,:),funny(:,:),    &
   p_g_co_pp(:,:,:),globma_pp(:),fun(:),store_pm_pp(:,:,:),newlo_pp(:),   &
   mass(:),globma_tmp(:,:),pmul_pp(:,:),utemp_pp(:,:)
 INTEGER,ALLOCATABLE::rest(:,:),g(:),num(:),g_num_pp(:,:),g_g_pp(:,:)         
!-----------------------input and initialisation--------------------------
 timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes)
 IF(numpe==npes)THEN
   CALL getname(argv,nlen)
   OPEN(10,FILE=argv(1:nlen)//'.dat',STATUS='OLD',ACTION='READ')
   READ(10,*)nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz,dtim,nstep,npri,val0
 END IF
 CALL bcast_inputdata_p125(numpe,npes,nels,nxe,nze,nip,aa,bb,cc,kx,ky,kz, &
   dtim,nstep,npri,val0)
 CALL calc_nels_pp
 ndof=nod*nodof
 nn_temp=0
 neq_temp=0
 ntot=ndof
 nye=nels/nxe/nze
 nr=(nxe+1)*(nye+1)+(nxe+1)*nze+nye*nze
 ALLOCATE(rest(nr,nodof+1),points(nip,ndim),weights(nip),kay(ndim,ndim),  &
   coord(nod,ndim),jac(ndim,ndim),p_g_co_pp(nod,ndim,nels_pp),            &
   der(ndim,nod),deriv(ndim,nod),g_num_pp(nod,nels_pp),kc(ntot,ntot),     &
   g(ntot),funny(1,nod),num(nod),g_g_pp(ntot,nels_pp),                    &
   store_pm_pp(ntot,ntot,nels_pp),mass(ntot),fun(nod),pm(ntot,ntot),      &
   globma_tmp(ntot,nels_pp),pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp))
 kay=0.0_iwp
 kay(1,1)=kx
 kay(2,2)=ky
 kay(3,3)=kz
 ielpe=iel_start
 CALL box_bc8(nxe,nye,nze,rest)
 CALL rearrange_2(rest)
!---------------loop the elements to  set up global arrays --------------- 
 elements_0: DO iel=1,nels_pp
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
 END DO elements_0
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
   WRITE(11,'(A,I5,A)')"This job ran on ",npes,"  processors" 
   WRITE(11,'(A)')"Global coordinates and node numbers" 
   DO i=1,nels_pp,nels_pp-1
     WRITE(11,'(A,I8)')"Element ",i
     num=g_num_pp(:,i)
     DO k=1,nod
       WRITE(11,'(A,I8,3E12.4)')"  Node",num(k),p_g_co_pp(k,:,i)
     END DO
   END DO
   WRITE(11,'(A,3(I8,A))')"There are ",nn," nodes",nr," restrained and",  &
     neq," equations"
   WRITE(11,*)"Time after setup  is  :",elap_time()-timest(1)
 END IF
 CALL sample(element,points,weights)
 globma_tmp=.0_iwp
 ALLOCATE(loads_pp(neq_pp),newlo_pp(neq_pp),globma_pp(neq_pp))
 loads_pp=.0_iwp
 newlo_pp=.0_iwp
 globma_pp=.0_iwp
!------------ loop the elements for integration and invert mass ----------  
 elements_1: DO iel=1,nels_pp
   coord=p_g_co_pp(:,:,iel)
   kc=0.0_iwp
   pm=0.0_iwp
   gauss_pts: DO i=1,nip
     CALL shape_der(der,points,i)
     CALL shape_fun(fun,points,i)
     funny(1,:)=fun(:)
     jac=MATMUL(der,coord)
     det=determinant(jac)
     CALL invert(jac)
     deriv=MATMUL(jac,der)
     kc=kc+MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
     pm=pm+MATMUL(TRANSPOSE(funny),funny)*det*weights(i) 
   END DO gauss_pts
   DO i=1,ntot
     mass(i)=sum(pm(i,:))
   END DO
   pm=.0_iwp
   DO i=1,ntot
     pm(i,i)=mass(i)
   END DO
   store_pm_pp(:,:,iel)=pm-kc*dtim
   DO i=1,ntot
     globma_tmp(i,iel)=globma_tmp(i,iel)+mass(i)
   END DO
 END DO elements_1
 IF(numpe==it)                                                            &
   WRITE(11,*)"Time after element integration is :",elap_time()-timest(1)
 CALL scatter(globma_pp,globma_tmp)
 globma_pp=1._iwp/globma_pp
 loads_pp=val0
 DEALLOCATE(globma_tmp)
!-------------------time stepping recursion-------------------------------
 IF(numpe==it)THEN
   WRITE(11,'(A)')"    Time     Pressure"
 END IF
 timesteps: DO j=1,nstep
   real_time=j*dtim              
!---------------  go round the elements  ---------------------------------
   pmul_pp=.0_iwp
   CALL gather(loads_pp,pmul_pp)
   utemp_pp=.0_iwp
   elements_2: DO iel=1,nels_pp  
     pm=store_pm_pp(:,:,iel) 
     utemp_pp(:,iel)=utemp_pp(:,iel)+MATMUL(pm,pmul_pp(:,iel)) 
   END DO elements_2 
   CALL scatter(newlo_pp,utemp_pp)
   loads_pp=newlo_pp*globma_pp
   newlo_pp=.0_iwp 
   IF(j/npri*npri==j.AND.numpe==it)                                       &
     WRITE(11,'(2E12.4)')real_time,loads_pp(is) 
 END DO timesteps
 IF(numpe==it)WRITE(11,*)"This analysis took  :",elap_time()-timest(1)
 CALL shutdown() 
END PROGRAM p125
