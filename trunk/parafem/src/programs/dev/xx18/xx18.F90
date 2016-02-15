PROGRAM xx18       
!------------------------------------------------------------------------- 
!      Program xx18 three dimensional analysis of an elastic solid
!      using 20-node brick elements, choice between built in 
!      preconditioned conjugate gradient solver and PETSc library
!      parallel version; loaded_nodes only; ARCHER eCSE06 project
!      
!      See program p121
!------------------------------------------------------------------------- 

! PETSc modules
! Use system, vectors, matrices, solvers, preconditioners
 USE petscsys; USE petscvec; USE petscmat; USE petscksp; USE petscpc

! PETSc won't work without MPI
!!USE mpi_wrapper  !remove comment for serial compilation
 USE precision; USE global_variables; USE mp_interface; USE input
 USE output; USE loading; USE timing; USE maths; USE gather_scatter
 USE steering; USE new_library; USE large_strain 
 IMPLICIT NONE

! PETSc types
! Use system, vectors, matrices, solvers, preconditioners
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petsckspdef.h>
#include <petsc/finclude/petscpcdef.h>

! neq,ntot are now global variables - must not be declared

 INTEGER,PARAMETER   :: nodof=3,ndim=3,nst=6
 INTEGER             :: loaded_nodes,iel,i,j,k,iters,limit,nn,nr,nip,nod
 INTEGER             :: nels,ndof
 INTEGER             :: npes_pp,node_end,node_start,nodes_pp,meshgen
 INTEGER             :: partitioner,nlen
 INTEGER             :: inewton ! Note - does not appear to be needed
 REAL(iwp),PARAMETER :: zero=0.0_iwp
 REAL(iwp)           :: rn0 ! Note - does not appear to be used
 REAL(iwp)           :: e,v,det,tol,up,alpha,beta,q
 LOGICAL             :: converged=.false.
 CHARACTER(LEN=50)   :: argv
 CHARACTER(LEN=15)   :: element
 CHARACTER(LEN=6)    :: ch 

! PETSc variables
 INTEGER,PARAMETER   :: p_max_string_length=1024
 PetscErrorCode      :: p_ierr
 Vec                 :: p_x,p_b,p_r
 Mat                 :: p_A
 KSP                 :: p_ksp
 PetscScalar         :: p_pr_n2,p_r_n2,p_b_n2
 PetscScalar,POINTER :: p_varray(:)
 PetscInt            :: p_its
 KSPConvergedReason  :: p_reason
 CHARACTER(LEN=p_max_string_length) :: p_description

!-------------------------------------------------------------------------
! 1. Dynamic arrays
!-------------------------------------------------------------------------

 REAL(iwp),ALLOCATABLE::points(:,:),dee(:,:),weights(:),val(:,:),        &
   disp_pp(:),g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:),bee(:,:),   &
   km(:,:),eps(:),sigma(:),diag_precon_pp(:),p_pp(:),r_pp(:),            &
   x_pp(:),xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),        &
   timest(:),diag_precon_tmp(:,:),eld_pp(:,:),temp(:)
 INTEGER,ALLOCATABLE::rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)

!-------------------------------------------------------------------------
! 2a. Input and initialisation
!-------------------------------------------------------------------------

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
   km(ntot,ntot),pmul_pp(ntot,nels_pp),                                  &
   utemp_pp(ntot,nels_pp),g_g_pp(ntot,nels_pp))

!-------------------------------------------------------------------------
! 2b. Start up PETSc after MPI has been started
!-------------------------------------------------------------------------

 CALL PetscInitialize(trim(argv)//".petsc",p_ierr)

!-------------------------------------------------------------------------
! 3. Find the steering array and equations per process
!-------------------------------------------------------------------------

 CALL rearrange(rest); g_g_pp=0; neq=0
 elements_0: DO iel=1,nels_pp
   CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
 END DO elements_0
 neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
 CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
 ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),        &
   u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp)); diag_precon_pp=zero
 p_pp=zero;  r_pp=zero;  x_pp=zero; xnew_pp=zero; u_pp=zero; d_pp=zero

!-------------------------------------------------------------------------
! 4. Element stiffness integration and global stiffness matrix creation
!-------------------------------------------------------------------------

! PETSc 32-bit indices and 64-bit reals
 CALL MatCreate(PETSC_COMM_WORLD,p_A,p_ierr)
 Call MatSetSizes(p_A,neq_pp,neq_pp,PETSC_DETERMINE,PETSC_DETERMINE,     &
                  p_ierr)
 CALL MatSetType(p_A,MATAIJ,p_ierr)
!- Block size fixed to 1 for just now - this cannot be set to nodof
!- until the restraints are handled block-wise in ParaFEM.
!- CALL MatSetBlockSize(p_A,nodof,p_ierr)

! Hack for now: For 20-node hexahedra and 3 dofs per node, then will
! be an (average) maximum of 81 * 3 = 243 entries per row
 CALL MatSeqAIJSetPreallocation(p_A,243,PETSC_NULL_INTEGER,p_ierr);
 CALL MatMPIAIJSetPreallocation(p_A,243,PETSC_NULL_INTEGER,              &
                                243,PETSC_NULL_INTEGER,p_ierr);

 dee=zero; CALL deemat(dee,e,v); CALL sample(element,points,weights)
 elements_1: DO iel=1,nels_pp
   km=zero
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
     det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
     CALL beemat(bee,deriv)
     km=km + MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
   END DO gauss_pts_1
   ! 1/ Use ubound(g_g_pp,1) instead of ntot? 2/ The following depends
   ! on g_g_pp holding 0 for erased rows/columns, and MatSetValues
   ! ignoring negative indices (PETSc always uses zero-based
   ! indexing). 3/ PETSc uses C array order, so a transpose is needed.
   CALL MatSetValues(p_A,ntot,g_g_pp(:,iel)-1,ntot,g_g_pp(:,iel)-1,      &
                     transpose(km),ADD_VALUES,p_ierr)
 END DO elements_1

 CALL MatAssemblyBegin(p_A,MAT_FINAL_ASSEMBLY,p_ierr)
 CALL MatAssemblyEnd(p_A,MAT_FINAL_ASSEMBLY,p_ierr)

!-------------------------------------------------------------------------
! 5. 
!-------------------------------------------------------------------------

 IF(numpe==1)THEN
   OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
   WRITE(11,'(A,I7,A)') "This job ran on ",npes," processes"
   WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes", nr, &
                           " restrained and ",neq," equations"
   WRITE(11,'(A,F10.4)') "Time to read input is:",timest(2)-timest(1)
   WRITE(11,'(A,F10.4)') "Time after setup is:",elap_time()-timest(1)
 END IF

!-------------------------------------------------------------------------
! 6. Get starting r
!-------------------------------------------------------------------------

 IF(loaded_nodes>0) THEN
   ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes)); node=0; val=zero
   CALL read_loads(argv,numpe,node,val)
   CALL load(g_g_pp,g_num_pp,node,val,r_pp); q=SUM_P(r_pp)
   IF(numpe==1) WRITE(11,'(A,E12.4)') "The total load is:",q
   DEALLOCATE(node,val)
 END IF

! Create a Vec.  For this particular order (create, set size, set type)
! the allocation is done by set type.
 CALL VecCreate(PETSC_COMM_WORLD,p_b,p_ierr)
 CALL VecSetSizes(p_b,neq_pp,PETSC_DECIDE,p_ierr)
 CALL VecSetType(p_b,VECSTANDARD,p_ierr)
!- Block size fixed to 1 for just now - this cannot be set to nodof
!- until the restraints are handled block-wise in ParaFEM.
!- CALL VecSetBlockSize(p_b,nodof,p_ierr)
 CALL VecGetArrayF90(p_b,p_varray,p_ierr)
 p_varray = r_pp
 CALL VecRestoreArrayF90(p_b,p_varray,p_ierr)

 DEALLOCATE(g_g_pp)

!-------------------------------------------------------------------------
! 7. Preconditioned conjugate gradient solver
!-------------------------------------------------------------------------

! d_pp=diag_precon_pp*r_pp; p_pp=d_pp; x_pp=zero

 timest(3) = elap_time()
 iters     = 0
 x_pp      = zero
 inewton   = 1    ! must be one in this program
 rn0       = zero ! note rn0 does not appear to be used in XX15

 CALL KSPCreate(PETSC_COMM_WORLD,p_ksp,p_ierr)
 CALL KSPSetOperators(p_ksp,p_A,p_A,p_ierr)
! KSP type, tolerances (rtol, abstol, dtol, maxits), KSP options, PC
! type, PC options are set in the xx18*.ppetsc file.  Those options are
! used to set up the preconditioned Krylov solver
 CALL KSPSetFromOptions(p_ksp,p_ierr)
! But the relative tolerance and maximum number of iterations are
! overriden by the value in the ParaFEM control file.  Note that relative
! tolerance for PETSc is for the preconditioned residual.
 CALL KSPSetTolerances(p_ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,  &
                       limit,p_ierr)

! Solution vector
 CALL VecDuplicate(p_b,p_x,p_ierr) 
! For non-linear solves, the previous solution will be used as an initial
! guess.  Set p_x to initial guess here.
 CALL KSPSetInitialGuessNonzero(p_ksp,PETSC_FALSE,p_ierr)
 CALL KSPSolve(p_ksp,p_b,p_x,p_ierr)

! Preconditioned residual L2 norm
 CALL KSPGetResidualNorm(p_ksp,p_pr_n2,p_ierr)
! True residual L2 norm
 CALL VecDuplicate(p_b,p_r,p_ierr)
 CALL KSPBuildResidual(p_ksp,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,p_r,    &
                       p_ierr)
 CALL VecNorm(p_r,NORM_2,p_r_n2,p_ierr)
 CALL VecDestroy(p_r,p_ierr)
! L2 norm of load
 CALL VecNorm(p_b,NORM_2,p_b_n2,p_ierr)
 CALL KSPGetIterationNumber(p_ksp,p_its,p_ierr)

 CALL KSPGetConvergedReason(p_ksp,p_reason,p_ierr)
 CALL p_describe_reason(p_reason,p_description)

! Copy PETSc solution vector to ParaFEM
 CALL VecGetArrayF90(p_x,p_varray,p_ierr)
 xnew_pp = p_varray
 CALL VecRestoreArrayF90(p_x,p_varray,p_ierr)

 CALL VecDestroy(p_x,p_ierr)
 CALL KSPDestroy(p_ksp,p_ierr)
 CALL VecDestroy(p_b,p_ierr)
 CALL MatDestroy(p_A,p_ierr)

! CALL PCG_VER1(inewton,limit,tol,storkm_pp,r_pp,diag_precon_pp, &
!               rn0,x_pp,iters)

!--------------------- preconditioned cg iterations ----------------------
! iters=0; timest(3)=elap_time()
! iterations: DO 
!   iters=iters+1; u_pp=zero; pmul_pp=zero; utemp_pp=zero
!   CALL gather(p_pp,pmul_pp)
!   elements_3: DO iel=1,nels_pp
!     utemp_pp(:,iel) = MATMUL(storkm_pp(:,:,iel),pmul_pp(:,iel))
!    CALL dgemv('n',ntot,ntot,1.0_iwp,storkm_pp(:,:,iel),ntot,           &
!               pmul_pp(:,iel),1,0.0_iwp,utemp_pp(:,iel),1)
!   END DO elements_3 ;CALL scatter(u_pp,utemp_pp)
!-------------------------- pcg equation solution ------------------------
!   up=DOT_PRODUCT_P(r_pp,d_pp); alpha=up/DOT_PRODUCT_P(p_pp,u_pp)
!   xnew_pp=x_pp+p_pp*alpha; r_pp=r_pp-u_pp*alpha
!   d_pp=diag_precon_pp*r_pp; beta=DOT_PRODUCT_P(r_pp,d_pp)/up
!   p_pp=d_pp+p_pp*beta; CALL checon_par(xnew_pp,tol,converged,x_pp)    
!   IF(converged.OR.iters==limit)EXIT
! END DO iterations

 IF(numpe==1)THEN
   WRITE(11,'(A,I6,A)')"The reason for convergence was ",p_reason,       &
                       " "//trim(p_description)
   WRITE(11,'(A,I6)')"The number of iterations to convergence was ",p_its
   WRITE(11,'(A,E17.7)')"The preconditioned residual L2 norm was ",p_pr_n2
   WRITE(11,'(A,E17.7)')"The true residual L2 norm ||b-Ax|| was  ",p_r_n2
   WRITE(11,'(A,E17.7)')"The relative error ||b-Ax||/||b|| was   ",      &
                        p_r_n2/p_b_n2
   WRITE(11,'(A,F10.4)')"Time to solve equations was  :",                &
                         elap_time()-timest(3)  
   WRITE(11,'(A,E12.4)')"The central nodal displacement is :",xnew_pp(1)
 END IF
!DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,storkm_pp,pmul_pp) 

!-------------------------------------------------------------------------
! 8. Recover stresses at centroidal gauss point
!-------------------------------------------------------------------------

 ALLOCATE(eld_pp(ntot,nels_pp)); eld_pp=zero; points=zero; nip=1; iel=1
 CALL gather(xnew_pp,eld_pp); DEALLOCATE(xnew_pp)
 IF(numpe==1)WRITE(11,'(A)')"The Centroid point stresses for element 1 are"
 gauss_pts_2: DO i=1,nip
   CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
   CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
   eps=MATMUL(bee,eld_pp(:,iel)); sigma=MATMUL(dee,eps)
   IF(numpe==1.AND.i==1) THEN
     WRITE(11,'(A,I5)')"Point ",i ; WRITE(11,'(6E12.4)') sigma
   END IF
 END DO gauss_pts_2; DEALLOCATE(g_coord_pp)

!-------------------------------------------------------------------------
! 9. Write out displacements
!-------------------------------------------------------------------------

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


!-------------------------------------------------------------------------
! 10. Shut down PETSc and ParaEFM
!-------------------------------------------------------------------------

 CALL PetscFinalize(p_ierr)
 CALL SHUTDOWN() 
 
CONTAINS
  SUBROUTINE p_describe_reason(p_reason,p_description)
    KSPConvergedReason, INTENT(IN)  :: p_reason
    CHARACTER(LEN=*), INTENT(OUT)   :: p_description
    
    ! The string constants in this routine have to be p_max_string_length or
    ! shorter.
    SELECT CASE (p_reason)
    ! Converged.
! not in PETSc Fortran    CASE (KSP_CONVERGED_RTOL_NORMAL)
! not in PETSc Fortran       p_description = "KSP_CONVERGED_RTOL_NORMAL"
! not in PETSc Fortran    CASE (KSP_CONVERGED_ATOL_NORMAL)
! not in PETSc Fortran       p_description = "KSP_CONVERGED_ATOL_NORMAL"
    CASE (KSP_CONVERGED_RTOL)
       p_description = "KSP_CONVERGED_RTOL (residual 2-norm decreased by a factor of rtol, from 2-norm of right hand side)"
    CASE (KSP_CONVERGED_ATOL)
       p_description = "KSP_CONVERGED_ATOL (residual 2-norm less than abstol)"
    CASE (KSP_CONVERGED_ITS)
       p_description = "KSP_CONVERGED_ITS (used by the preonly preconditioner that always uses ONE iteration, or when the KSPConvergedSkip() convergence test routine is set)"
    CASE (KSP_CONVERGED_CG_NEG_CURVE)
       p_description = "KSP_CONVERGED_CG_NEG_CURVE"
    CASE (KSP_CONVERGED_CG_CONSTRAINED)
       p_description = "KSP_CONVERGED_CG_CONSTRAINED"
    CASE (KSP_CONVERGED_STEP_LENGTH)
       p_description = "KSP_CONVERGED_STEP_LENGTH"
    CASE (KSP_CONVERGED_HAPPY_BREAKDOWN)
       p_description = "KSP_CONVERGED_HAPPY_BREAKDOWN"
    ! Diverged.
    CASE (KSP_DIVERGED_NULL)
       p_description = "KSP_DIVERGED_NULL"
    CASE (KSP_DIVERGED_ITS)
       p_description = "KSP_DIVERGED_ITS (required more than its to reach convergence)"
    CASE (KSP_DIVERGED_DTOL)
       p_description = "KSP_DIVERGED_DTOL (residual norm increased by a factor of divtol)"
    CASE (KSP_DIVERGED_BREAKDOWN)
       p_description = "KSP_DIVERGED_BREAKDOWN (generic breakdown in method)"
    CASE (KSP_DIVERGED_BREAKDOWN_BICG)
       p_description = "KSP_DIVERGED_BREAKDOWN_BICG (Initial residual is orthogonal to preconditioned initial residual. Try a different preconditioner, or a different initial level.)"
    CASE (KSP_DIVERGED_NONSYMMETRIC)
       p_description = "KSP_DIVERGED_NONSYMMETRIC"
    CASE (KSP_DIVERGED_INDEFINITE_PC)
       p_description = "KSP_DIVERGED_INDEFINITE_PC"
    CASE (KSP_DIVERGED_NANORINF)
       p_description = "KSP_DIVERGED_NANORINF (residual norm became NaN or Inf likely due to 0/0)"
    CASE (KSP_DIVERGED_INDEFINITE_MAT)
       p_description = "KSP_DIVERGED_INDEFINITE_MAT"
! not yet in PETSc Fortran    CASE (KSP_DIVERGED_PCSETUP_FAILED)
! not yet in PETSc Fortran       p_description = "KSP_DIVERGED_PCSETUP_FAILED"
    ! Still iterating.
    CASE (KSP_CONVERGED_ITERATING)
       p_description = "KSP_CONVERGED_ITERATING"
    ! Unknown  
    CASE DEFAULT
       p_description = "reason not known"
    END SELECT
  END SUBROUTINE p_describe_reason
END PROGRAM xx18
