PROGRAM xx8       
!------------------------------------------------------------------------- 
!      Program xx8 three dimensional analysis of an elastic solid
!      using 20-node brick elements, preconditioned conjugate gradient
!      solver; diagonal preconditioner diag_precon; parallel version
!      loaded_nodes only
!------------------------------------------------------------------------- 
!USE mpi_wrapper  !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE steering; USE new_library; IMPLICIT NONE
! neq,ntot are now global variables - must not be declared
 INTEGER,PARAMETER::nodof=3,ndim=3,nst=6
 INTEGER::loaded_nodes,iel,i,j,k,iters,limit,nn,nr,nip,nod,nels,ndof,    &
   npes_pp,node_end,node_start,nodes_pp,meshgen,partitioner,nlen
 REAL(iwp),PARAMETER::zero=0.0_iwp
 REAL(iwp)::e,v,det,tol,up,alpha,beta,q; LOGICAL::converged=.false.
 CHARACTER(LEN=50)::argv; CHARACTER(LEN=15)::element; CHARACTER(LEN=6)::ch 
!---------------------------- dynamic arrays -----------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),dee(:,:),weights(:),val(:,:),        &
   disp_pp(:),g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:),bee(:,:),   &
   storkm_pp(:,:,:),eps(:),sigma(:),diag_precon_pp(:),p_pp(:),r_pp(:),   &
   x_pp(:),xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),        &
   timest(:),diag_precon_tmp(:,:),eld_pp(:,:),temp(:)
 INTEGER,ALLOCATABLE::rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)
!------------------------ input and initialisation -----------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen) 
 CALL read_p121(argv,numpe,e,element,limit,loaded_nodes,meshgen,nels,    &
   nip,nn,nod,nr,partitioner,tol,v)
 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof
 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),           &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest); timest(2)=elap_time()
 ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),der(ndim,nod),    &
   deriv(ndim,nod),bee(nst,ntot),weights(nip),eps(nst),sigma(nst),       &
   storkm_pp(ntot,ntot,nels_pp),pmul_pp(ntot,nels_pp),                   &
   utemp_pp(ntot,nels_pp),g_g_pp(ntot,nels_pp))
!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_0: DO iel=1,nels_pp
   CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_0
 neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
 ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),        &
   u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp)); diag_precon_pp=zero
 p_pp=zero;  r_pp=zero;  x_pp=zero; xnew_pp=zero; u_pp=zero; d_pp=zero
!------ element stiffness integration and build the preconditioner -------
 dee=zero; CALL deemat(dee,e,v); CALL sample(element,points,weights)
 storkm_pp=zero
 elements_1: DO iel=1,nels_pp
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
     det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
     CALL beemat(bee,deriv)
     storkm_pp(:,:,iel)=storkm_pp(:,:,iel) +                             &
                    MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)   
   END DO gauss_pts_1
 END DO elements_1
 ALLOCATE(diag_precon_tmp(ntot,nels_pp)); diag_precon_tmp=zero
 elements_2: DO iel=1,nels_pp ; DO i=1,ndof
   diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel)+storkm_pp(i,i,iel)
 END DO;  END DO elements_2
 CALL scatter(diag_precon_pp,diag_precon_tmp); DEALLOCATE(diag_precon_tmp)
 IF(numpe==1)THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A,I7,A)') "This job ran on ",npes," processes"
   WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes", nr, &
                           " restrained and ",neq," equations"
   WRITE(11,'(A,F10.4)') "Time to read input is:",timest(2)-timest(1)
   WRITE(11,'(A,F10.4)') "Time after setup is:",elap_time()-timest(1)
 END IF
!----------------------------- get starting r ----------------------------
 IF(loaded_nodes>0) THEN
   ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes)); node=0; val=zero
   CALL read_loads(argv,numpe,node,val)
   CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:)); q=SUM_P(r_pp(1:))
   IF(numpe==1) WRITE(11,'(A,E12.4)') "The total load is:",q
   DEALLOCATE(node,val)
 END IF
 DEALLOCATE(g_g_pp); diag_precon_pp=1._iwp/diag_precon_pp
 d_pp=diag_precon_pp*r_pp; p_pp=d_pp; x_pp=zero
!--------------------- preconditioned cg iterations ----------------------
 iters=0; timest(3)=elap_time()
 iterations: DO 
   iters=iters+1; u_pp=zero; pmul_pp=zero; utemp_pp=zero
   CALL gather(p_pp,pmul_pp)
   elements_3: DO iel=1,nels_pp
     utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
!    CALL dgemv('n',ntot,ntot,1.0_iwp,storkm_pp(:,:,iel),ntot,           &
!               pmul_pp(:,iel),1,0.0_iwp,utemp_pp(:,iel),1)
   END DO elements_3 ;CALL scatter(u_pp,utemp_pp)
!-------------------------- pcg equation solution ------------------------
   up=DOT_PRODUCT_P(r_pp,d_pp); alpha=up/DOT_PRODUCT_P(p_pp,u_pp)
   xnew_pp=x_pp+p_pp*alpha; r_pp=r_pp-u_pp*alpha
   d_pp=diag_precon_pp*r_pp; beta=DOT_PRODUCT_P(r_pp,d_pp)/up
   p_pp=d_pp+p_pp*beta; CALL checon_par(xnew_pp,tol,converged,x_pp)    
   IF(converged.OR.iters==limit)EXIT
 END DO iterations
 IF(numpe==1)THEN
   WRITE(11,'(A,I6)')"The number of iterations to convergence was ",iters
   WRITE(11,'(A,F10.4)')"Time to solve equations was  :",                &
                         elap_time()-timest(3)  
   WRITE(11,'(A,E12.4)')"The central nodal displacement is :",xnew_pp(1)
 END IF
 DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,storkm_pp,pmul_pp) 
!--------------- recover stresses at centroidal gauss point --------------
 ALLOCATE(eld_pp(ntot,nels_pp)); eld_pp=zero; points=zero; nip=1; iel=1
 CALL gather(xnew_pp(1:),eld_pp); DEALLOCATE(xnew_pp)
 IF(numpe==1)WRITE(11,'(A)')"The Centroid point stresses for element 1 are"
 gauss_pts_2: DO i=1,nip
   CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
   CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
   eps=MATMUL(bee,eld_pp(:,iel)); sigma=MATMUL(dee,eps)
   IF(numpe==1.AND.i==1) THEN
     WRITE(11,'(A,I5)')"Point ",i ; WRITE(11,'(6E12.4)') sigma
   END IF
 END DO gauss_pts_2; DEALLOCATE(g_coord_pp)
!------------------------ write out displacements ------------------------
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 IF(numpe==1) THEN;  WRITE(ch,'(I6.6)') numpe
   OPEN(12,file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',       &
     action='write')
   WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
   WRITE(12,'(A/A/A)') "part", "     1","coordinates"
 END IF
 ALLOCATE(disp_pp(nodes_pp*ndim),temp(nodes_pp)); disp_pp=zero; temp=zero
 CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,          &
                    node_start,node_end,eld_pp,disp_pp,1)
 DO i=1,ndim ; temp=zero
   DO j=1,nodes_pp; k=i+(ndim*(j-1)); temp(j)=disp_pp(k); END DO
   CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
 END DO ; IF(numpe==1) CLOSE(12)
 IF(numpe==1) WRITE(11,'(A,F10.4)')"This analysis took  :",              &
   elap_time()-timest(1)  
 CALL SHUTDOWN() 
END PROGRAM xx8
