PROGRAM xx17    

 !/****f* dev/xx17
 !*  NAME
 !*    PROGRAM: xx17
 !*  SYNOPSIS
 !*    Usage:   main program
 !*  FUNCTION
 !*    Steady state 3-d Navier-Stokes equation using 20-node velocity
 !*    hexahedral elements coupled to 8-node pressure hexahedral elements:
 !*    u-p-v-w order.  Choice between built in BiCGStab(l) solver and PETSc
 !*    library (BiCGStab(l) or other.  ARCHER eCSE06 project.
 !*    
 !*    Parallel version.  See program p126.
 !*    
 !*  AUTHORS
 !*    Lee Margetts, Mark Filipiak
 !*  CREATION DATE
 !*    11.06.2015
 !*  MODIFICATION HISTORY
 !*    Version 2, 11.01.2017, Mark Filipiak
 !*  COPYRIGHT
 !*    (c) University of Manchester 2007-2015, University of Edinburgh 2017
 !******
 !*  Place remarks that should not be included in the documentation here.
 !*
 !*/

 USE choose_solvers
 USE parafem_petsc
 USE precision; USE global_variables; USE mp_interface; USE input; 
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE steering; USE fluid; USE new_library; USE bicg
 IMPLICIT NONE

 ! To use PETSc directly and use parafem_petsc, include only the PETSc Fortran
 ! types.  The Fortran interfaces for the most useful PETSc subroutines are
 ! included in the parafem_petsc module: see parafem_petsc.F90 for the
 ! interfaces that are included.  If you need the interfaces for others then add
 ! the includes after the PETSc types.
 ! PETSc types
 !#include <petsc/finclude/petscdef.h>
 ! Use Fortran interfaces for the PETSc routines not included by the
 ! parafem_petsc module, for example, to use PETSc viewers:
 !#include <petsc/finclude/petscviewer.h>
 !#include <petsc/finclude/petscviewer.h90>

!neq,ntot are now global variables - not declared
 INTEGER,PARAMETER::nodof=4,nod=20,nodf=8,ndim=3
 INTEGER::nn,nip,cj_tot,i,j,k,l,iel,ell,limit,fixed_freedoms,iters,      &
   cjits,nr,n_t,fixed_freedoms_pp,nres,is,it,nlen,nels,ndof,     &
   npes_pp,meshgen,partitioner,node_end,node_start,nodes_pp,             &
   fixed_freedoms_start,cjiters=0  
 REAL(iwp):: visc,rho,det,ubar,vbar,wbar,tol,cjtol,alpha,      &
   penalty,x0,pp,kappa
 REAL(iwp),PARAMETER::zero=0.0_iwp,one=1.0_iwp
 LOGICAL::converged
 CHARACTER(LEN=15)::element='hexahedron'; CHARACTER(LEN=50)::argv
 CHARACTER(len=2)::iters_s

 CHARACTER(len=choose_solvers_string_length) :: solvers
 LOGICAL                  :: error
 CHARACTER(:),ALLOCATABLE :: message
 CHARACTER,PARAMETER      :: tab = ACHAR(9)
 REAL                     :: peak_memory_use

!--------------------------- dynamic arrays ------------------------------
 REAL(iwp),ALLOCATABLE::points(:,:),derivf(:,:),fun(:),store_pp(:),      &
   jac(:,:),kay(:,:),der(:,:),deriv(:,:),weights(:),derf(:,:),funf(:),   &
   coordf(:,:),g_coord_pp(:,:,:),c11(:,:),c21(:,:),c12(:,:),val(:),      &
   wvel(:),c23(:,:),c32(:,:),x_pp(:),b_pp(:),temp(:),                    &
   funny(:,:),row1(:,:),row2(:,:),uvel(:),vvel(:),funnyf(:,:),rowf(:,:), &
   storke_pp(:,:,:),diag_pp(:),utemp_pp(:,:),xold_pp(:),c24(:,:),        &
   c42(:,:),row3(:,:),                                                   &
   diag_tmp(:,:),pmul_pp(:,:),timest(:),upvw_pp(:),ke(:,:,:),val_pp(:)
 INTEGER,ALLOCATABLE::rest(:,:),g_num_pp(:,:),g_g_pp(:,:),no(:),g_t(:),  &
   no_pp(:),no_pp_temp(:)
!---------------------- input and initialisation -------------------------
 ALLOCATE(timest(40)); timest=zero

 timest(1) = elap_time()
 timest(2) = elap_time()

 CALL find_pe_procs(numpe,npes)

 peak_memory_use = p_memory_peak()
 IF (numpe == 1) WRITE(*,'(A,F7.2,A)')                                        &
   "peak memory use at start:           ", peak_memory_use, " GB"

 CALL getname(argv,nlen)
 CALL read_p126(argv,numpe,cjits,cjtol,ell,fixed_freedoms,kappa,limit,   &
   meshgen,nels,nip,nn,nr,nres,partitioner,penalty,rho,tol,x0,visc) 

 solvers = get_solvers()
 IF (.NOT. solvers_valid(solvers)) THEN
   CALL shutdown
 END IF

 CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
 ntot=nod+nodf+nod+nod; n_t=nod*nodof
 ALLOCATE(g_num_pp(nod,nels_pp),g_coord_pp(nod,ndim,nels_pp),            &
   rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
 CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
 IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
 CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
 CALL read_rest(argv,numpe,rest)

 timest(30) = timest(30) + elap_time()-timest(2) ! 30 = read
 timest(2) = elap_time()

 ALLOCATE(points(nip,ndim),derivf(ndim,nodf),pmul_pp(ntot,nels_pp),      &
   jac(ndim,ndim),kay(ndim,ndim),der(ndim,nod),deriv(ndim,nod),          &
   derf(ndim,nodf),funf(nodf),coordf(nodf,ndim),funny(nod,1),            &
   g_g_pp(ntot,nels_pp),c11(nod,nod),c12(nod,nodf),c21(nodf,nod),        &
   c24(nodf,nod),c42(nod,nodf),c32(nod,nodf),fun(nod),row2(1,nod),       &
   c23(nodf,nod),uvel(nod),vvel(nod),row1(1,nod),funnyf(nodf,1),         &
   rowf(1,nodf),no_pp_temp(fixed_freedoms),wvel(nod),row3(1,nod),        &
   ! ke is (ntot,ntot,1) so that formupvw can be used (formupvw cannot be
   ! changed until a new edition of the book).
   ke(ntot,ntot,1),g_t(n_t),                                &
   no(fixed_freedoms),val(fixed_freedoms),weights(nip),                  &
   diag_tmp(ntot,nels_pp),utemp_pp(ntot,nels_pp))

 IF (solvers == parafem_solvers) THEN
   ALLOCATE(storke_pp(ntot,ntot,nels_pp))
 END IF

!----------  find the steering array and equations per process -----------
 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_1: DO iel=1,nels_pp
   CALL find_g3(g_num_pp(:,iel),g_t,rest)
   CALL g_t_g_ns(nod,g_t,g_g_pp(:,iel))
 END DO elements_1
 neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp)
 CALL make_ggl(npes_pp,npes,g_g_pp)
 IF (solvers == parafem_solvers) THEN
   DEALLOCATE(g_g_pp)
 END IF

 DO i=1,neq_pp; IF(nres==ieq_start+i-1)THEN;it=numpe;is=i;END IF;END DO
 IF(numpe==it) THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
 END IF

 ALLOCATE(x_pp(neq_pp),b_pp(neq_pp),                                     &
   diag_pp(neq_pp),xold_pp(neq_pp),                                      &
   store_pp(neq_pp))
 x_pp=zero; b_pp=zero; diag_pp=zero
 xold_pp=zero; store_pp=zero

 !-----------------------------------------------------------------------------
 ! 4. Start up PETSc after find_pe_procs (so that MPI has been started) and
 !    after find_g3 so that neq_pp and g_g_pp are set up.
 !-----------------------------------------------------------------------------
 IF (solvers == petsc_solvers) THEN
   CALL p_initialize(argv,error)
   IF (error) THEN
     CALL shutdown
   END IF
   ! Set up PETSc.
   CALL p_setup(ntot,g_g_pp,error)
   IF (error) THEN
     CALL p_finalize
     CALL shutdown
   END IF
 END IF
 peak_memory_use = p_memory_peak()
 IF (numpe == 1) WRITE(*,'(A,F7.2,A)')                                         &
   "peak memory use after setup:        ", peak_memory_use," GB "

 timest(31) = timest(31) + elap_time()-timest(2) ! 31 = setup
 timest(2) = elap_time()

!-------------------------- organise fixed equations ---------------------
 CALL read_loads_ns(argv,numpe,no,val)
 CALL reindex(ieq_start,no,no_pp_temp,fixed_freedoms_pp,                       &
   fixed_freedoms_start,neq_pp)
 ALLOCATE(no_pp(1:fixed_freedoms_pp),val_pp(1:fixed_freedoms_pp))
 no_pp = no_pp_temp(1:fixed_freedoms_pp)
 val_pp = val(fixed_freedoms_start:fixed_freedoms_start+fixed_freedoms_pp-1)
 DEALLOCATE(no_pp_temp,val)

 timest(30) = timest(30) + elap_time()-timest(2) ! 30 = read
 timest(2) = elap_time()

!------------------------- main iteration loop ---------------------------
 CALL sample(element,points,weights); uvel=zero; vvel=zero; wvel=zero
 kay=zero; iters=0; cj_tot=0; kay(1,1)=visc/rho; kay(2,2)=visc/rho
 kay(3,3)=visc/rho; timest(3)=elap_time()

 timest(31) = timest(31) + elap_time()-timest(2) ! 31 = setup
 timest(2) = elap_time()

 iterations: DO

   timest(2) = elap_time()

   IF(numpe==1) THEN
     WRITE(iters_s,'(i2.2)') iters
     OPEN(12,file=argv(1:nlen)//".ensi.VEL-"//iters_s,status='replace',        &
       action='write')
     WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
     WRITE(12,'(A/A/A)') "part", "     1","coordinates"
   END IF
   CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
   ALLOCATE(upvw_pp(nodes_pp*nodof),temp(nodes_pp)); upvw_pp=zero          
   temp=zero; utemp_pp=zero; CALL gather(x_pp(1:),utemp_pp)
   CALL scatter_nodes_ns(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,          &
     node_start,node_end,utemp_pp,upvw_pp,1)
   DO i=1,nodof ; temp=zero
     IF(i/=2) THEN
       DO j=1,nodes_pp; k=i+(nodof*(j-1)); temp(j)=upvw_pp(k); END DO
       CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp); END IF
   END DO; IF(numpe==1) CLOSE(12)
   DEALLOCATE(upvw_pp,temp)         
     
   timest(35) = timest(35) + elap_time()-timest(2) ! 35 = write
   timest(2) = elap_time()

   iters=iters+1

   IF (solvers == parafem_solvers) THEN
     storke_pp=zero
   ELSE IF (solvers == petsc_solvers) THEN
     CALL p_zero_matrix
   END IF

   diag_pp=zero; utemp_pp=zero
   b_pp=zero; pmul_pp=zero; CALL gather(x_pp,utemp_pp)
   CALL gather(xold_pp,pmul_pp)

   timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in main loop
   timest(2) = elap_time()

!-------------------- element stiffness integration ----------------------
   elements_2: DO iel=1,nels_pp
     uvel=(utemp_pp(1:nod,iel)+pmul_pp(1:nod,iel))*.5_iwp
     DO i=nod+nodf+1,nod+nodf+nod
       vvel(i-nod-nodf)=(utemp_pp(i,iel)+pmul_pp(i,iel))*.5_iwp
     END DO
     DO i=nod+nodf+nod+1,ntot
       wvel(i-nod-nodf-nod)=(utemp_pp(i,iel)+pmul_pp(i,iel))*.5_iwp
     END DO                                                        
     c11=zero; c12=zero; c21=zero; c23=zero; c32=zero; c24=zero; c42=zero
     gauss_points_1: DO i=1,nip
!------------------------ velocity contribution --------------------------
       CALL shape_fun(funny(:,1),points,i)
       ubar=DOT_PRODUCT(funny(:,1),uvel);vbar=DOT_PRODUCT(funny(:,1),vvel)
       wbar=DOT_PRODUCT(funny(:,1),wvel)
       IF(iters==1)THEN; ubar=one; vbar=zero; wbar=zero; END IF
       CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
       det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
       row1(1,:)=deriv(1,:); row2(1,:)=deriv(2,:); row3(1,:)=deriv(3,:)
       c11=c11+                                                          &
           MATMUL(MATMUL(TRANSPOSE(deriv),kay),deriv)*det*weights(i) +   &
           MATMUL(funny,row1)*det*weights(i)*ubar +                      &
           MATMUL(funny,row2)*det*weights(i)*vbar +                      &
           MATMUL(funny,row3)*det*weights(i)*wbar
!------------------------ pressure contribution --------------------------
       CALL shape_fun(funnyf(:,1),points,i); CALL shape_der(derf,points,i)
       coordf(1:4,:)=g_coord_pp(1:7:2,:,iel)
       coordf(5:8,:)=g_coord_pp(13:19:2,:,iel); jac=MATMUL(derf,coordf)
       det=determinant(jac); CALL invert(jac); derivf=MATMUL(jac,derf)
       rowf(1,:)=derivf(1,:)
       c12=c12+MATMUL(funny,rowf)*det*weights(i)/rho
       rowf(1,:)=derivf(2,:)
       c32=c32+MATMUL(funny,rowf)*det*weights(i)/rho
       rowf(1,:)=derivf(3,:)
       c42=c42+MATMUL(funny,rowf)*det*weights(i)/rho
       c21=c21+MATMUL(funnyf,row1)*det*weights(i)
       c23=c23+MATMUL(funnyf,row2)*det*weights(i) 
       c24=c24+MATMUL(funnyf,row3)*det*weights(i)
     END DO gauss_points_1
      ! ke is (ntot,ntot,1) so that formupvw can be used (formupvw cannot be
      ! changed until a new edition of the book).
      CALL formupvw(ke,1,c11,c12,c21,c23,c32,c24,c42)
      IF (solvers == parafem_solvers) THEN
        storke_pp(:,:,iel) = ke(:,:,1)
      ELSE IF (solvers == petsc_solvers) THEN
        CALL p_add_element(g_g_pp(:,iel),ke(:,:,1))
      END IF
   END DO elements_2
   peak_memory_use = p_memory_peak()
   IF (numpe == 1) WRITE(*,'(A,F7.2,A)')                                        &
     "peak memory use after add elements: ", peak_memory_use," GB "
 
   IF (solvers == petsc_solvers) THEN
     CALL p_assemble
   END IF
   peak_memory_use = p_memory_peak()
   IF (numpe == 1) WRITE(*,'(A,F7.2,A)')                                        &
     "peak memory use after assemble:     ", peak_memory_use," GB "

   timest(33) = timest(33) + elap_time()-timest(2) ! 33 = matrix assemble
   timest(2) = elap_time()

!----------------------- build the preconditioner ------------------------
   IF (solvers == parafem_solvers) THEN
     diag_tmp=zero
     elements_2a: DO iel=1,nels_pp; DO k=1,ntot
       diag_tmp(k,iel)=diag_tmp(k,iel)+storke_pp(k,k,iel); END DO
     END DO elements_2a; CALL scatter(diag_pp,diag_tmp)
   END IF

   timest(34) = timest(34) + elap_time()-timest(2) ! 34 = solve
   timest(2) = elap_time()

!------------------- prescribed values of velocity and pressure ----------
   IF (solvers == parafem_solvers) THEN
     DO i=1,fixed_freedoms_pp; k=no_pp(i)-ieq_start+1
       diag_pp(k)=diag_pp(k)+penalty
       b_pp(k)=diag_pp(k)*val_pp(i)
       store_pp(k)=diag_pp(k)
     END DO   
   ELSE IF (solvers == petsc_solvers) THEN
     ! Notwithstanding the name, the penalty method is not used for the ParaFEM
     ! solver;  in fact the fixed freedom rows are zeroed apart from the
     ! diagonal, the penalty added to the diagonal entry and the RHS adjusted.
     ! So use PETSc's equivalent, which SETS the diagonal entry to the penalty
     ! and therefore requires a different adjustment to the RHS.
     CALL p_zero_rows(no_pp,penalty,val_pp,b_pp)
   END IF

   timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in main loop
   timest(2) = elap_time()

!---------- solve the equations element-by-element using BiCGSTAB --------
   IF(iters==1) x_pp = x0

   IF (solvers == parafem_solvers) THEN
     CALL bicgstabl_p(b_pp,cjiters,cjits,cjtol,ell,kappa,ieq_start,no_pp,      &
                      pmul_pp,store_pp,storke_pp,x_pp)

     b_pp=x_pp-xold_pp; pp=norm_p(b_pp); cj_tot=cj_tot+cjiters

     IF(numpe==it) THEN
       WRITE(11,'(A,I6,A)') "It took BiCGSTAB(L) ", cjiters,               &
          " iterations to converge"
     END IF
   ELSE IF (solvers == petsc_solvers) THEN
     CALL p_solve(b_pp,x_pp,initial_guess_nonzero=.TRUE.)
     CALL p_print_info(it,11)

     b_pp=x_pp-xold_pp; pp=norm_p(b_pp)

   END IF
   IF(numpe==it) THEN
     WRITE(11,'(A,E12.4)') "Norm of the error is:", pp
     FLUSH(11)
   END IF

   timest(34) = timest(34) + elap_time()-timest(2) ! 34 = solve
   timest(2) = elap_time()

   CALL checon_par(x_pp,tol,converged,xold_pp)

   timest(32) = timest(32) + elap_time()-timest(2) ! 32 = other work in main loop
   timest(2) = elap_time()

   IF(converged.OR.iters==limit)EXIT

 END DO iterations

 DEALLOCATE(b_pp,diag_pp,xold_pp,store_pp)
 DEALLOCATE(pmul_pp)
 IF (solvers == parafem_solvers) THEN
   DEALLOCATE(storke_pp)
 END IF
 IF (solvers == petsc_solvers) THEN
   ! g_g_pp is needed by PETSc, so don't deallocate until finished.
   DEALLOCATE(g_g_pp)
 END IF
!------------------------- output results --------------------------------
 IF(numpe==1) THEN
   WRITE(iters_s,'(i2.2)') iters
   OPEN(12,file=argv(1:nlen)//".ensi.VEL-"//iters_s,status='replace',          &
     action='write')
   WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
   WRITE(12,'(A/A/A)') "part", "     1","coordinates"
 END IF
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 ALLOCATE(upvw_pp(nodes_pp*nodof),temp(nodes_pp)); upvw_pp=zero          
 temp=zero; utemp_pp=zero; CALL gather(x_pp(1:),utemp_pp)
 CALL scatter_nodes_ns(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,            &
   node_start,node_end,utemp_pp,upvw_pp,1)
 DO i=1,nodof ; temp=zero
   IF(i/=2) THEN
     DO j=1,nodes_pp; k=i+(nodof*(j-1)); temp(j)=upvw_pp(k); END DO
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp); END IF
 END DO; IF(numpe==1) CLOSE(12)
 DEALLOCATE(upvw_pp,temp)         

 IF(numpe==1) THEN
   OPEN(12,file=argv(1:nlen)//".ensi.VEL-final",status='replace',              &
     action='write')
   WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
   WRITE(12,'(A/A/A)') "part", "     1","coordinates"
 END IF
 CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
 ALLOCATE(upvw_pp(nodes_pp*nodof),temp(nodes_pp)); upvw_pp=zero          
 temp=zero; utemp_pp=zero; CALL gather(x_pp(1:),utemp_pp)
 CALL scatter_nodes_ns(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,            &
   node_start,node_end,utemp_pp,upvw_pp,1)
 DO i=1,nodof ; temp=zero
   IF(i/=2) THEN
     DO j=1,nodes_pp; k=i+(nodof*(j-1)); temp(j)=upvw_pp(k); END DO
     CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp); END IF
 END DO; IF(numpe==1) CLOSE(12)
 DEALLOCATE(upvw_pp,temp)         

 IF(numpe==1) THEN
   WRITE(iters_s,'(i2.0)') iters + 1
   OPEN(12,file=argv(1:nlen)//".ensi.case",status='replace',action='write')
   WRITE(12,'(A)') "#"
   WRITE(12,'(A)') "# Post-processing file generated by subroutine WRITE_ENSI in "
   WRITE(12,'(A)') "# Smith, Griffiths and Margetts, 'Programming the Finite Element Method',"
   WRITE(12,'(A)') "# Wiley, 2014."
   WRITE(12,'(A)') "#"
   WRITE(12,'(A)') "# Ensight Gold Format"
   WRITE(12,'(A)') "#"
   WRITE(12,'(A)') "# Problem name: "//argv(1:nlen)
   WRITE(12,'(A)') "#"
   WRITE(12,'(A)') "FORMAT"
   WRITE(12,'(A)') "type:  ensight gold"
   WRITE(12,'(A)') "GEOMETRY"
   WRITE(12,'(A)') "model: 1  "//argv(1:nlen)//".ensi.geo"
   WRITE(12,'(A)') "VARIABLE"
   WRITE(12,'(A)') "scalar per element:  material      "//argv(1:nlen)//".ensi.MATID"
   WRITE(12,'(A)') "scalar per node:     restraint     "//argv(1:nlen)//".ensi.NDBND"
   WRITE(12,'(A)') "vector per node:     velocity      "//argv(1:nlen)//".ensi.VEL-**"
   WRITE(12,'(A)') "vector per node:     load          "//argv(1:nlen)//".ensi.NDLDS"
   WRITE(12,'(A)') "TIME"
   WRITE(12,'(A)') "time set:     1"
   WRITE(12,'(A)') "number of steps:    "//iters_s
   WRITE(12,'(A)') "filename start number:    0"
   WRITE(12,'(A)') "filename increment:    1"
   WRITE(12,'(A)') "time values:"
   WRITE(12,'(F5.1)') (REAL(i),i=0,iters)
   CLOSE(12)
 END IF
     
 timest(35) = timest(35) + elap_time()-timest(2) ! 35 = write
 timest(2) = elap_time()

 peak_memory_use = p_memory_peak()

 IF (numpe==it) THEN
   WRITE(11,'(A)') "The pressure at the corner of the box is: "
   WRITE(11,'(A)') "Freedom  Pressure "
   WRITE(11,'(I6,E12.4)') nres, x_pp(is)
   IF (solvers == parafem_solvers) THEN
     WRITE(11,'(A,I6)')"The total number of BiCGSTAB iterations was:",cj_tot
   END IF
   WRITE(11,'(A,I5,A)')"The solution took",iters," iterations to converge"

   IF (solvers == parafem_solvers) THEN
     WRITE(11,'(A)') "ParaFEM revision "//REVISION
   ELSE IF (solvers == petsc_solvers) THEN
     WRITE(11,'(A)') "ParaFEM revision "//REVISION//"; "//p_version()
   END IF
   WRITE(11,'(A,I5,A)') "This job ran on ",npes," processors"
   WRITE(11,'(A,3(I8,A))') "There are ",nn," nodes",nels," elements and ",    &
                           neq," equations"
   WRITE(11,'(A,F10.4)') "Time to read input:       ", timest(30)
   WRITE(11,'(A,F10.4)') "Time for setup:           ", timest(31)
   WRITE(11,'(A,F10.4)') "Time for matrix assemble: ", timest(33)
   WRITE(11,'(A,F10.4)') "Time for linear solve:    ", timest(34)
   WRITE(11,'(A,F10.4)') "Other time in main loop:  ", timest(32)
   WRITE(11,'(A,F10.4)') "Time to write results:    ", timest(35)
   WRITE(11,'(A)')       "                          ----------"
   WRITE(11,'(A,F10.4)') "Total:                    ", SUM(timest(30:35))
   WRITE(11,'(A)')       "                          ----------"
   WRITE(11,'(A,F10.4)') "This analysis took:       ", elap_time()-timest(1)
   WRITE(11,*)
   WRITE(11,'(A,F10.2,A)') "Peak memory use: ",peak_memory_use," GB"
   WRITE(11,*)
   ! REVISION is substituted with a string like "2108" (including the double
   ! quotes) by the preprocessor, see the makefile.
   WRITE(11,'(2A,2(I0,A),4(F0.2,A),F0.2)',advance='no')                       &
     REVISION,tab,                                                            &
     npes,tab,                                                                &
     neq,tab,                                                                 &
     timest(31),tab,                                                          &
     timest(33),tab,                                                          &
     timest(34),tab,                                                          &
     timest(31)+timest(33)+timest(34),tab,                                    &
     peak_memory_use
   IF (solvers == petsc_solvers) THEN
     WRITE(11,'(2A)') tab,p_version()
   END IF
   CLOSE(11)
 END IF

 IF (solvers == petsc_solvers) THEN
   CALL p_shutdown
 END IF
 CALL SHUTDOWN() 

END PROGRAM xx17
