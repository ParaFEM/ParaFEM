PROGRAM p121       
!------------------------------------------------------------------------- 
!      Program 12.1 three dimensional analysis of an elastic solid
!      using 20-node brick elements, preconditioned conjugate gradient
!      solver; diagonal preconditioner diag_precon; parallel version
!      loaded_nodes only
!------------------------------------------------------------------------- 
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE partition; USE elements; USE steering; USE pcg; IMPLICIT NONE
! neq,ntot are now global variables - must not be declared
 INTEGER,PARAMETER::nodof=3,ndim=3,nst=6
 INTEGER::loaded_nodes,fixed_freedoms,iel,i,j,k,l,idx1,idx2,iters,limit, &
   nn,nr,nip,nod,nels,ndof,npes_pp,node_end,node_start,nodes_pp,argc,    &
   iargc,meshgen,partitioner,fixed_freedoms_pp,fixed_freedoms_start
 REAL(iwp),PARAMETER::zero=0.0_iwp,penalty=1.0e20_iwp
 REAL(iwp)::e,v,det,tol,up,alpha,beta,tload
 CHARACTER(LEN=50)::program_name='p121',fname,job_name,label
 CHARACTER(LEN=15)::element ;LOGICAL::converged=.false.
!---------------------------- dynamic arrays -----------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),dee(:,:),weights(:),val(:,:),        &
   disp_pp(:),g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:),bee(:,:),   &
   storkm_pp(:,:,:),eld(:),eps(:),sigma(:),diag_precon_pp(:),p_pp(:),    &
   r_pp(:),x_pp(:),xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:),        &
   d_pp(:),timest(:),diag_precon_tmp(:,:),eld_pp(:,:),tensor_pp(:,:,:),  &
   valf(:),store_pp(:),fun(:),shape_integral_pp(:,:),                    &
   stress_integral_pp(:,:),stressnodes_pp(:),principal_integral_pp(:,:), &
   princinodes_pp(:),principal(:),reacnodes_pp(:),rest(:,:),             &
   g_num_pp(:,:),g_g_pp(:,:),node(:)
 INTEGER,  ALLOCATABLE::no(:),no_pp(:),no_pp_temp(:),no_global(:),sense(:)
!------------------------ input and initialisation -----------------------
 ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
 CALL find_pe_procs(numpe,npes); argc=iargc(); CALL GETARG(1, job_name) 
 CALL read_p121(job_name,numpe,e,element,fixed_freedoms,limit,           &
   loaded_nodes,meshgen,nels,nip,nn,nod,nr,partitioner,tol,v)
 CALL calc_nels_pp(job_name,nels,npes,numpe,partitioner,nels_pp)
 ndof=nod*nodof; ntot=ndof
 ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),           &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp2(job_name,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(job_name,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(job_name,numpe,rest); timest(2)=elap_time()
 ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),principal(ndim),  &
   der(ndim,nod),deriv(ndim,nod),bee(nst,ntot),weights(nip),fun(nod),    &
   storkm_pp(ntot,ntot,nels_pp),eld(ntot),eps(nst),sigma(nst),           &
   pmul_pp(ntot,nels_pp),utemp_pp(ntot,nels_pp),g_g_pp(ntot,nels_pp))
!----------  find the steering array and equations per process ----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_1: DO iel = 1, nels_pp
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_1
!CALL MPI_BARRIER(MPI_COMM_WORLD,ier) ! remove
 elements_2: DO iel = 1, nels_pp  
    i=MAXVAL(g_g_pp(:,iel)); IF(i>neq) neq=i
 END DO elements_2  
 neq=MAX_INTEGER_P(neq); CALL calc_neq_pp; CALL calc_npes_pp(npes,npes_pp)
 CALL make_ggl2(npes_pp,npes,g_g_pp)
!CALL MPI_BARRIER(MPI_COMM_WORLD,ier) ! remove
 ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),        &
   u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp))
 p_pp=zero;  r_pp=zero;  x_pp=zero; xnew_pp=zero; u_pp=zero; d_pp=zero
 diag_precon_pp=zero
!------ element stiffness integration and build the preconditioner ------
 dee=zero; CALL deemat(e,v,dee); CALL sample(element,points,weights)
 storkm_pp=zero
 elements_3: DO iel=1,nels_pp
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
     det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
     CALL beemat(deriv,bee)
     storkm_pp(:,:,iel)=storkm_pp(:,:,iel) +                             &
                    MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)   
   END DO gauss_pts_1
 END DO elements_3
 ALLOCATE(diag_precon_tmp(ntot,nels_pp)); diag_precon_tmp=zero
 elements_4: DO iel=1,nels_pp ; DO i=1,ndof
   diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel)+storkm_pp(i,i,iel)
 END DO;  END DO elements_4
 CALL scatter(diag_precon_pp,diag_precon_tmp); DEALLOCATE(diag_precon_tmp)
 IF(numpe==1)THEN
   OPEN(11,FILE='p121.res',STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A,I7,A)') "This job ran on ",npes," processes"
   WRITE(11,'(A,3(I8,A))') "There are ",nn," nodes", nr, &
                           " restrained and ",neq," equations"
   WRITE(11,*) "The time to read input is:",timest(2)-timest(1)
   WRITE(11,*) "The time to after setup is:",elap_time()-timest(1)
 END IF
!----------------------------- get starting r ----------------------------
 IF(loaded_nodes>0) THEN
   ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes)); val=zero; node=0
   CALL read_loads(job_name,numpe,node,val)
   CALL load(g_g_pp,g_num_pp,node,val,r_pp(1:))
   tload = SUM_P(r_pp(1:)); DEALLOCATE(node,val)
 END IF
 DEALLOCATE(g_g_pp); diag_precon_pp=1._iwp/diag_precon_pp
 d_pp=diag_precon_pp*r_pp; p_pp=d_pp; x_pp=zero
!--------------------- preconditioned cg iterations ----------------------
 iters=0; timest(3)=elap_time()
 iterations: DO 
   iters=iters+1; u_pp=zero; pmul_pp=zero; utemp_pp=zero
   CALL gather(p_pp,pmul_pp)
   elements_5: DO iel=1,nels_pp
     utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
   END DO elements_5 ;CALL scatter(u_pp,utemp_pp)
!-------------------------- pcg equation solution ------------------------
   up=DOT_PRODUCT_P(r_pp,d_pp); alpha=up/DOT_PRODUCT_P(p_pp,u_pp)
   xnew_pp=x_pp+p_pp*alpha; r_pp=r_pp-u_pp*alpha
   d_pp=diag_precon_pp*r_pp; beta=DOT_PRODUCT_P(r_pp,d_pp)/up
   p_pp=d_pp+p_pp*beta  
   CALL checon_par(xnew_pp,tol,converged,x_pp)    
   IF(converged.OR.iters==limit)EXIT
 END DO iterations
 IF(numpe==1)THEN
   WRITE(11,'(A,I6)')"The number of iterations to convergence was ",iters
   WRITE(11,*)"Time to solve equations was  :", elap_time()-timest(3)  
   WRITE(11,'(A,E12.4)')"The central nodal displacement is :",xnew_pp(1)
 END IF
 DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,storkm_pp,pmul_pp) 
!------------------------ write out displacements ------------------------
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 IF(numpe==1) THEN
   fname = job_name(1:INDEX(job_name, " ")-1)//".dis"
   OPEN(24, file=fname, status='replace', action='write')
 END IF
 ALLOCATE(eld_pp(ntot,nels_pp)); eld_pp=zero
 CALL gather(xnew_pp(1:),eld_pp); DEALLOCATE(xnew_pp)
 ALLOCATE(disp_pp(nodes_pp*ndim)); disp_pp = zero
 label="*DISPLACEMENT"
 CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,          &
                    node_start,node_end,eld_pp,disp_pp,1)
 CALL write_nodal_variable(label,24,1,nodes_pp,npes,numpe,ndim,disp_pp)
 DEALLOCATE(disp_pp); IF(numpe==1) CLOSE(24)
!--------------- recover stresses at centroidal gauss point --------------
 nip=1; points=zero; iel=1
 IF(numpe==1)WRITE(11,'(A)')"The Centroid point stresses for element 1 are"
 gauss_pts_2: DO i=1,nip
   CALL shape_der(der,points,i); jac=MATMUL(der,coord)
   CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
   eps=MATMUL(bee,eld); sigma=MATMUL(dee,eps)
   IF(numpe==1.AND.i==1)THEN
     WRITE(11,'(A,I5)')"Point ",i ; WRITE(11,'(6E12.4)') sigma
   END IF
 END DO gauss_pts_2
 IF(numpe==1)WRITE(11,*)"This analysis took  :", elap_time()-timest(1)  
 CALL shutdown() 
END PROGRAM p121
