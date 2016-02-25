PROGRAM p125      
!-------------------------------------------------------------------------
!   Program 12.5 conduction equation on a 3-d box volume using 8-node
!   hexahedral elements and a simple explicit algorithm : parallel version
!   write on processor it at freedom nres
!-------------------------------------------------------------------------
!USE mpi_wrapper !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE geometry; USE new_library; IMPLICIT NONE
! neq,ntot are now global variables - not declared
 INTEGER,PARAMETER::nodof=1,ndim=3
 INTEGER::nels,ndof,npes_pp,nn,nr,nip,nod,i,j,k,iel,nstep,nlen,          &
   npri,nres,it,is,meshgen,partitioner,loaded_nodes=0,fixed_freedoms,    &
   nodes_pp,node_start,node_end
 REAL(iwp)::kx,ky,kz,det,dtim,val0,real_time
 REAL(iwp),PARAMETER::zero=0.0_iwp; CHARACTER(LEN=6)::ch
 CHARACTER(LEN=15)::element; CHARACTER(LEN=50)::argv
!------------------------------ dynamic arrays --------------------------- 
 REAL(iwp),ALLOCATABLE::loads_pp(:),points(:,:),kay(:,:),jac(:,:),       &
   der(:,:),deriv(:,:),weights(:),kc(:,:),funny(:,:),globma_pp(:),       &
   fun(:),store_pm_pp(:,:,:),newlo_pp(:),mass(:),globma_tmp(:,:),        &
   pmul_pp(:,:),utemp_pp(:,:),g_coord_pp(:,:,:),timest(:),pm(:,:),       &
   ptl_pp(:)
 INTEGER,ALLOCATABLE::rest(:,:),g_num_pp(:,:),g_g_pp(:,:)         
!-------------------------- input and initialisation ---------------------
 ALLOCATE(timest(25)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)
 CALL read_p125(argv,numpe,dtim,element,fixed_freedoms,kx,ky,kz,         &
   loaded_nodes,meshgen,nels,nip,nn,nod,npri,nr,nres,nstep,partitioner,  &
   val0)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof
 ALLOCATE(g_num_pp(nod,nels_pp),g_coord_pp(nod,ndim,nels_pp),            &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen==2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest);timest(2)=elap_time()
 ALLOCATE (points(nip,ndim),weights(nip),kay(ndim,ndim),jac(ndim,ndim),  &
   der(ndim,nod),deriv(ndim,nod),kc(ntot,ntot),funny(1,nod),             &
   g_g_pp(ntot,nels_pp),store_pm_pp(ntot,ntot,nels_pp),                  &
   mass(ntot),fun(nod),pm(ntot,ntot),globma_tmp(ntot,nels_pp),           &
   pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp))
!----------  find the steering array and equations per process -----------
 CALL rearrange_2(rest); g_g_pp = 0; neq=0
 elements_0: DO iel = 1, nels_pp
   CALL find_g4(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_0
 neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
 DO i=1,neq_pp; IF(nres==ieq_start+i-1)THEN;it=numpe;is=i;END IF;END DO
 IF(numpe==it)THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')              
   WRITE(11,'(A,I5,A)')"This job ran on ",npes,"  processes" 
   WRITE(11,'(A,3(I12,A))')"There are ",nn," nodes",nr," restrained and",&
     neq," equations"
   WRITE(11,'(A,F10.4)') "Time to read input is:",timest(2)-timest(1)
   WRITE(11,'(A,F10.4)') "Time after setup is:",elap_time()-timest(1)
 END IF
 ALLOCATE(loads_pp(neq_pp),newlo_pp(neq_pp),globma_pp(neq_pp))
 loads_pp=zero; newlo_pp=zero; globma_pp=zero; globma_tmp=zero
!------------- loop the elements for integration and invert mass ---------
 timest(3)=elap_time(); CALL sample(element,points,weights)
 kay=zero; kay(1,1)=kx; kay(2,2)=ky; kay(3,3)=kz
 elements_1: DO iel=1,nels_pp
   kc=zero; pm=zero
    gauss_pts: DO i=1,nip
      CALL shape_der(der,points,i); CALL shape_fun(fun,points,i)
      funny(1,:)=fun(:); jac=MATMUL(der,g_coord_pp(:,:,iel))
      det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
      kc=kc+MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i)
      pm=pm+MATMUL(TRANSPOSE(funny),funny)*det*weights(i) 
    END DO gauss_pts
    DO i=1,ntot; mass(i)=sum(pm(i,:)); END DO
    pm=zero; DO i=1,ntot; pm(i,i)=mass(i); END DO
    store_pm_pp(:,:,iel)=pm-kc*dtim
    DO i=1,ntot; globma_tmp(i,iel)=globma_tmp(i,iel)+mass(i); END DO
  END DO elements_1
  IF(numpe==it) THEN; WRITE(11,'(A,F10.4)')                              &
    "Time for element integration is :",elap_time()-timest(3); END IF
  CALL scatter(globma_pp,globma_tmp); globma_pp = 1._iwp/globma_pp
  loads_pp=val0; DEALLOCATE(globma_tmp); timest(4)=elap_time()
!------------------------- time stepping recursion -----------------------
 IF(numpe==it)THEN
   WRITE(11,'(A)')"  Time        Pressure"
   WRITE(11,'(2E12.4)') 0.0_iwp, loads_pp(is)
 END IF
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 ALLOCATE(ptl_pp(nodes_pp*ndim))
 timesteps: DO j=1,nstep
   real_time = j*dtim              
!------------------------- go round the elements -------------------------
   utemp_pp=zero; pmul_pp=zero; CALL gather(loads_pp,pmul_pp)
   elements_2: DO iel=1,nels_pp  
     pm = store_pm_pp(:,:,iel) 
     utemp_pp(:,iel) = utemp_pp(:,iel)+MATMUL(pm,pmul_pp(:,iel)) 
   END DO elements_2; CALL scatter(newlo_pp,utemp_pp)
   loads_pp=newlo_pp*globma_pp; newlo_pp=zero 
   IF(j/npri*npri==j) THEN
     IF(numpe==it) WRITE(11,'(2E12.4)')real_time,loads_pp(is)
!--------------------- output pressures for ParaView ---------------------
     IF(numpe==1)THEN; WRITE(ch,'(I6.6)') j
       OPEN(12,file=argv(1:nlen)//".ensi.NDPRE-"//ch,status='replace',   &
         action='write')
       WRITE(12,'(A)')"Alya Ensight Gold --- Scalar per-node variable file"
       WRITE(12,'(A/A/A)') "part", "    1","coordinates"
     END IF
     ptl_pp=zero; utemp_pp=zero; CALL gather(loads_pp(1:),utemp_pp)
     CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,     &
       node_start,node_end,utemp_pp,ptl_pp,1)
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,ptl_pp); CLOSE(12)
   END IF
 END DO timesteps; timest(5)=elap_time()
 IF(numpe==it) THEN
   WRITE(11,'(A,F10.4)')"Time stepping recursion took  :",               &
     timest(5)-timest(4)
   WRITE(11,'(A,F10.4)')"This analysis took  :",elap_time()-timest(1)
 END IF; CALL SHUTDOWN() 
END PROGRAM p125
