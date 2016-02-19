PROGRAM xx18       
  
  !/****f* dev/xx18
  !*  NAME
  !*    PROGRAM: xx18
  !*  SYNOPSIS
  !*    Usage:   main program
  !*  FUNCTION
  !*    Three dimensional analysis of an elastic solid using 20-node brick
  !*    elements, choice between built in preconditioned conjugate gradient
  !*    solver and PETSc library.  ARCHER eCSE06 project.
  !*    
  !*    Parallel version.  Loaded_nodes only.  See program p121.
  !*    
  !*    Subroutine                 Purpose
  !*    
  !*    p_describe_reason          Returns a description of the reason
  !*                               for convergence in PETSc
  !*    p_row_nnz                  Returns an approximate upper estimate
  !*                               of the number of non zeroes in a row
  !*                               of the global matrix
  !*  AUTHORS
  !*    Lee Margetts, Mark Filipiak
  !*  CREATION DATE
  !*    02.01.2007
  !*  MODIFICATION HISTORY
  !*    Version 2, 19.02.2016, Mark Filipiak
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010, University of Edinburgh 2016
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  ! PETSc modules
  ! Use system, vectors, matrices, solvers, preconditioners
  USE petscsys
  USE petscvec
  USE petscmat
  USE petscksp
  USE petscpc
  
  ! PETSc will always use MPI (unless PETSc itself has been compiled without
  ! MPI)
  !USE mpi_wrapper  !remove comment for serial compilation
  USE PRECISION; USE global_variables; USE mp_interface; USE input
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
  INTEGER             :: inewton
  REAL(iwp),PARAMETER :: zero=0.0_iwp
  REAL(iwp)           :: rn0
  REAL(iwp)           :: e,v,det,tol,up,alpha,beta,q
  LOGICAL             :: converged=.FALSE.
  CHARACTER(LEN=50)   :: argv
  CHARACTER(LEN=15)   :: element
  CHARACTER(LEN=6)    :: ch 

  ! PETSc variables
  INTEGER,PARAMETER   :: p_max_string_length=1024
  ! p_over_allocation should be a run time setting with a default set in the
  ! program.  Choosing too small a value (e.g. 0.99) slows MatSetValues to a
  ! crawl.
  REAL,PARAMETER      :: p_over_allocation=1.3
  PetscErrorCode      :: p_ierr
  Vec                 :: p_x,p_b,p_r
  Mat                 :: p_A
  KSP                 :: p_ksp
  PetscScalar         :: p_pr_n2,p_r_n2,p_b_n2
  PetscInt,ALLOCATABLE    :: p_rows(:),p_cols(:)
  PetscScalar,ALLOCATABLE :: p_values(:)
  PetscScalar,POINTER :: p_varray(:)
  PetscInt            :: p_its,p_nnz
  DOUBLE PRECISION    :: p_info(MAT_INFO_SIZE)
  KSPConvergedReason  :: p_reason
  CHARACTER(LEN=p_max_string_length) :: p_description
  
  !-----------------------------------------------------------------------------
  ! 1. Dynamic arrays
  !-----------------------------------------------------------------------------
  
  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:),val(:,:),           &
    disp_pp(:),g_coord_pp(:,:,:),jac(:,:),der(:,:),deriv(:,:),bee(:,:),        &
    ! storkm_pp is used only if the built in preconditioned conjugate gradient
    ! is used
!!$    storkm_pp(:,:,:),km(:,:),eps(:),sigma(:),diag_precon_pp(:),p_pp(:),        &
    km(:,:),eps(:),sigma(:),diag_precon_pp(:),p_pp(:),                         &
    r_pp(:),x_pp(:),xnew_pp(:),u_pp(:),pmul_pp(:,:),utemp_pp(:,:),d_pp(:),     &
    timest(:),diag_precon_tmp(:,:),eld_pp(:,:),temp(:)
  INTEGER,ALLOCATABLE :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)
  
  !-----------------------------------------------------------------------------
  ! 2. Input and initialisation
  !-----------------------------------------------------------------------------
  
  ALLOCATE(timest(20)); timest=zero; timest(1)=elap_time()
  CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen) 
  CALL read_p121(argv,numpe,e,element,limit,loaded_nodes,meshgen,nels,         &
                 nip,nn,nod,nr,partitioner,tol,v)
  CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)
  ndof=nod*nodof; ntot=ndof
  ALLOCATE(g_num_pp(nod, nels_pp),g_coord_pp(nod,ndim,nels_pp),                &
    rest(nr,nodof+1)); g_num_pp=0; g_coord_pp=zero; rest=0
  CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
  IF(meshgen == 2) CALL abaqus2sg(element,g_num_pp)
  CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
  CALL read_rest(argv,numpe,rest); timest(2)=elap_time()
  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim),der(ndim,nod),         &
    deriv(ndim,nod),bee(nst,ntot),weights(nip),eps(nst),sigma(nst),            &
!!$    storkm_pp(ntot,ntot,nels_pp),km(ntot,ntot),pmul_pp(ntot,nels_pp),          &
    km(ntot,ntot),pmul_pp(ntot,nels_pp),                                       &
    utemp_pp(ntot,nels_pp),g_g_pp(ntot,nels_pp))
  
  !-----------------------------------------------------------------------------
  ! 3. Start up PETSc after MPI has been started
  !-----------------------------------------------------------------------------
  
  CALL PetscInitialize(TRIM(argv)//".petsc",p_ierr)
  
  !-----------------------------------------------------------------------------
  ! 4. Find the steering array and equations per process
  !-----------------------------------------------------------------------------

  CALL rearrange(rest); g_g_pp=0; neq=0
  elements_0: DO iel=1,nels_pp
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest)
  END DO elements_0
  neq=MAXVAL(g_g_pp); neq=max_p(neq); CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp); CALL make_ggl(npes_pp,npes,g_g_pp)
  ALLOCATE(p_pp(neq_pp),r_pp(neq_pp),x_pp(neq_pp),xnew_pp(neq_pp),             &
    u_pp(neq_pp),d_pp(neq_pp),diag_precon_pp(neq_pp)); diag_precon_pp=zero
  p_pp=zero;  r_pp=zero;  x_pp=zero; xnew_pp=zero; u_pp=zero; d_pp=zero
  
  !-----------------------------------------------------------------------------
  ! 5. Element stiffness integration and global stiffness matrix creation
  !-----------------------------------------------------------------------------
  
  ! PETSc 64-bit indices and 64-bit reals.  In most (all?) places, passing a
  ! 32-bit integer where an intent(in) 64-integer is required is safe, because
  ! the PETSc Fortran-C interface de-references the pointer that is actually
  ! passed, even though it does not specify any intent in the Fortran interface.
  ! This is not safe when passing arrays, they need to be copied to PetscInt
  ! arrays.  And for safety the same should be done for the PetscScalar arrays.
  CALL MatCreate(PETSC_COMM_WORLD,p_A,p_ierr)
  CALL MatSetSizes(p_A,neq_pp,neq_pp,PETSC_DETERMINE,PETSC_DETERMINE,p_ierr)
  CALL MatSetType(p_A,MATAIJ,p_ierr)
  !- Block size fixed to 1 for just now - this cannot be set to nodof until the
  !- restraints are handled block-wise in ParaFEM.
  !- CALL MatSetBlockSize(p_A,nodof,p_ierr)
  
  ! Find an approximate number of zeroes per row for the matrix size
  ! pre-allocation.
  CALL p_row_nnz(nodof,ndim,nod,p_over_allocation,p_nnz)
  CALL MatSeqAIJSetPreallocation(p_A,p_nnz,PETSC_NULL_INTEGER,p_ierr)
  CALL MatMPIAIJSetPreallocation(p_A,p_nnz,PETSC_NULL_INTEGER,                 &
                                     p_nnz,PETSC_NULL_INTEGER,p_ierr)
  ! If the allocation is too small, PETSc will produce reams of information and
  ! not construct the matrix properly.  We output some information at the end if
  ! p_over_allocation should be increased.
  CALL MatSetOption(p_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,p_ierr)
  
  ! Arrays of PETSc types to be proof against changes in index and scalar sizes.
  ALLOCATE(p_rows(ntot),p_cols(ntot),p_values(ntot*ntot))
  dee=zero; CALL deemat(dee,e,v); CALL sample(element,points,weights)
!!$  storkm_pp=zero
  elements_1: DO iel=1,nels_pp
    km=zero
    gauss_pts_1: DO i=1,nip
      CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
      det=determinant(jac); CALL invert(jac); deriv=MATMUL(jac,der)
      CALL beemat(bee,deriv)
      km=km + MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
    END DO gauss_pts_1
!!$    storkm_pp(:,:,iel)=km(:,:)
    ! 1/ Use ubound(g_g_pp,1) instead of ntot? 2/ The following depends on
    ! g_g_pp holding 0 for erased rows/columns, and MatSetValues ignoring
    ! negative indices (PETSc always uses zero-based indexing). 3/ PETSc uses C
    ! array order, so a transpose is needed.
    p_rows   = g_g_pp(:,iel) - 1
    p_cols   = p_rows
    p_values = RESHAPE(TRANSPOSE(km),(/ntot*ntot/))
    CALL MatSetValues(p_A,ntot,p_rows,ntot,p_cols,p_values,ADD_VALUES,p_ierr)
  END DO elements_1
  DEALLOCATE(p_rows,p_cols,p_values)
  
  CALL MatAssemblyBegin(p_A,MAT_FINAL_ASSEMBLY,p_ierr)
  CALL MatAssemblyEnd(p_A,MAT_FINAL_ASSEMBLY,p_ierr)
  
  CALL MatGetInfo(p_A,MAT_GLOBAL_SUM,p_info,p_ierr)
  IF(numpe==1)THEN
    IF(p_info(MAT_INFO_MALLOCS)/=0.0)THEN
      WRITE(*,'(A,I0,A)') "The matrix assembly required ",                     &
                          NINT(p_info(MAT_INFO_MALLOCS)),                      &
                          " mallocs.  Increase p_over_allocation to speed up " &
                          //"the assembly."
    END IF
  END IF
  
  !-----------------------------------------------------------------------------
  ! 6. Build the preconditioner (not for PETSc)
  !-----------------------------------------------------------------------------

!!$  ALLOCATE(diag_precon_tmp(ntot,nels_pp)); diag_precon_tmp=zero
!!$  elements_2: DO iel=1,nels_pp ; DO i=1,ndof
!!$    diag_precon_tmp(i,iel) = diag_precon_tmp(i,iel)+storkm_pp(i,i,iel)
!!$  END DO elements_2;  END DO elements_2
!!$  CALL scatter(diag_precon_pp,diag_precon_tmp); DEALLOCATE(diag_precon_tmp)
  
  IF(numpe==1)THEN
    OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
    WRITE(11,'(A,I7,A)') "This job ran on ",npes," processes"
    WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes", nr,                     &
                             " restrained and ",neq," equations"
    WRITE(11,'(A,F10.4)') "Time to read input is:",timest(2)-timest(1)
    WRITE(11,'(A,F10.4)') "Time after setup is:",elap_time()-timest(1)
  END IF
  
  !-----------------------------------------------------------------------------
  ! 7. Get starting r
  !-----------------------------------------------------------------------------
  
  IF(loaded_nodes>0) THEN
    ALLOCATE(node(loaded_nodes),val(ndim,loaded_nodes)); node=0; val=zero
    CALL read_loads(argv,numpe,node,val)
    CALL load(g_g_pp,g_num_pp,node,val,r_pp); q=SUM_P(r_pp)
    IF(numpe==1) WRITE(11,'(A,E12.4)') "The total load is:",q
    DEALLOCATE(node,val)
  END IF
  
  ! RHS vector.  For this particular order (create, set size, set type) the
  ! allocation is done by set type.
  CALL VecCreate(PETSC_COMM_WORLD,p_b,p_ierr)
  CALL VecSetSizes(p_b,neq_pp,PETSC_DECIDE,p_ierr)
  CALL VecSetType(p_b,VECSTANDARD,p_ierr)
  !- Block size fixed to 1 for just now - this cannot be set to nodof until the
  !- restraints are handled block-wise in ParaFEM.
  !- CALL VecSetBlockSize(p_b,nodof,p_ierr)
  CALL VecGetArrayF90(p_b,p_varray,p_ierr)
  ! This is OK as long as PetscScalars not smaller than ParaFEM reals.  There
  ! should be a test for sizes of PetscScalars (and PetscInts) and ParaFEM reals
  ! (and indices).
  p_varray = r_pp
  CALL VecRestoreArrayF90(p_b,p_varray,p_ierr)
  
!!$  DEALLOCATE(g_g_pp); diag_precon_pp=1._iwp/diag_precon_pp
  DEALLOCATE(g_g_pp)
  
  !-----------------------------------------------------------------------------
  ! 8. Preconditioned Krylov solver
  !-----------------------------------------------------------------------------
  
  timest(3) = elap_time()
  x_pp      = zero
!!$  iters     = 0
!!$  inewton   = 1    ! must be one in this program
!!$  rn0       = zero
  
  CALL KSPCreate(PETSC_COMM_WORLD,p_ksp,p_ierr)
  CALL KSPSetOperators(p_ksp,p_A,p_A,p_ierr)
  ! KSP type, tolerances (rtol, abstol, dtol, maxits), KSP options, PC type, PC
  ! options are set in the xx18*.ppetsc file.  Those options are used to set up
  ! the preconditioned Krylov solver
  CALL KSPSetFromOptions(p_ksp,p_ierr)
  ! But the relative tolerance and maximum number of iterations are overriden by
  ! the value in the ParaFEM control file.  Note that relative tolerance for
  ! PETSc is for the preconditioned residual.  What is the ParaFEM iteration
  ! limit: non-linear iterations or solver iterations?
  CALL KSPSetTolerances(p_ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,       &
                        limit,p_ierr)
  
  ! Solution vector
  CALL VecDuplicate(p_b,p_x,p_ierr) 
  ! For non-linear solves, the previous solution will be used as an initial
  ! guess: copy the ParaFEM solution vector to PETSc.
  CALL VecGetArrayF90(p_x,p_varray,p_ierr)
  ! This is OK as long as PetscScalars are not smaller than ParaFEM reals.
  ! There should be a test for sizes of PetscScalars (and PetscInts) and ParaFEM
  ! reals (and indices).
  p_varray = x_pp
  CALL VecRestoreArrayF90(p_x,p_varray,p_ierr)
  CALL KSPSetInitialGuessNonzero(p_ksp,PETSC_TRUE,p_ierr)
  CALL KSPSolve(p_ksp,p_b,p_x,p_ierr)
  
  ! Preconditioned residual L2 norm
  CALL KSPGetResidualNorm(p_ksp,p_pr_n2,p_ierr)
  ! True residual L2 norm
  CALL VecDuplicate(p_b,p_r,p_ierr)
  CALL KSPBuildResidual(p_ksp,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,p_r,p_ierr)
  CALL VecNorm(p_r,NORM_2,p_r_n2,p_ierr)
  CALL VecDestroy(p_r,p_ierr)
  ! L2 norm of load
  CALL VecNorm(p_b,NORM_2,p_b_n2,p_ierr)

  CALL KSPGetIterationNumber(p_ksp,p_its,p_ierr)
  CALL KSPGetConvergedReason(p_ksp,p_reason,p_ierr)
  CALL p_describe_reason(p_reason,p_description)
  
  ! Copy PETSc solution vector to ParaFEM
  CALL VecGetArrayF90(p_x,p_varray,p_ierr)
  ! This is OK as long as ParaFEM reals not smaller than PetscScalars.  There
  ! should be a test for sizes of PetscScalars (and PetscInts) and ParaFEM reals
  ! (and indices).
  xnew_pp = p_varray
  CALL VecRestoreArrayF90(p_x,p_varray,p_ierr)
  
  CALL VecDestroy(p_x,p_ierr)
  CALL KSPDestroy(p_ksp,p_ierr)
  CALL VecDestroy(p_b,p_ierr)
  CALL MatDestroy(p_A,p_ierr)
  
!!$  CALL PCG_VER1(inewton,limit,tol,storkm_pp,r_pp,diag_precon_pp,rn0,x_pp,iters)
!!$  xnew_pp = x_pp

  IF(numpe==1)THEN
    WRITE(11,'(A,I0,A)') "The reason for convergence was ",p_reason,           &
                         " "//TRIM(p_description)
    WRITE(11,'(A,I0)') "The number of iterations to convergence was ",p_its
    WRITE(11,'(A,E17.7)') "The preconditioned residual L2 norm was ",p_pr_n2
    WRITE(11,'(A,E17.7)') "The true residual L2 norm ||b-Ax|| was  ",p_r_n2
    WRITE(11,'(A,E17.7)') "The relative error ||b-Ax||/||b|| was   ",          &
                          p_r_n2/p_b_n2
!!$    WRITE(11,'(A,I6)') "The number of iterations to convergence was ",iters
    WRITE(11,'(A,F10.4)') "Time to solve equations was  :",                    &
                          elap_time()-timest(3)  
    WRITE(11,'(A,E12.4)') "The central nodal displacement is :",xnew_pp(1)
  END IF

!!$  DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,storkm_pp,pmul_pp) 
  DEALLOCATE(p_pp,r_pp,x_pp,u_pp,d_pp,diag_precon_pp,pmul_pp) 
  
  !-----------------------------------------------------------------------------
  ! 9. Recover stresses at centroidal gauss point
  !-----------------------------------------------------------------------------
  
  ALLOCATE(eld_pp(ntot,nels_pp)); eld_pp=zero; points=zero; nip=1; iel=1
  CALL gather(xnew_pp,eld_pp); DEALLOCATE(xnew_pp)
  IF(numpe==1)WRITE(11,'(A)') "The Centroid point stresses for element 1 are"
  gauss_pts_2: DO i=1,nip
    CALL shape_der(der,points,i); jac=MATMUL(der,g_coord_pp(:,:,iel))
    CALL invert(jac); deriv=MATMUL(jac,der); CALL beemat(bee,deriv)
    eps=MATMUL(bee,eld_pp(:,iel)); sigma=MATMUL(dee,eps)
    IF(numpe==1.AND.i==1) THEN
      WRITE(11,'(A,I5)') "Point ",i ; WRITE(11,'(6E12.4)') sigma
    END IF
  END DO gauss_pts_2; DEALLOCATE(g_coord_pp)
  
  !-----------------------------------------------------------------------------
  ! 10. Write out displacements
  !-----------------------------------------------------------------------------
  
  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)
  IF(numpe==1) THEN;  WRITE(ch,'(I6.6)') numpe
    OPEN(12,file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',            &
         action='write')
    WRITE(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
    WRITE(12,'(A/A/A)') "part", "     1","coordinates"
  END IF
  ALLOCATE(disp_pp(nodes_pp*ndim),temp(nodes_pp)); disp_pp=zero; temp=zero
  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,               &
                     node_start,node_end,eld_pp,disp_pp,1)
  DO i=1,ndim ; temp=zero
    DO j=1,nodes_pp
      k=i+(ndim*(j-1)); temp(j)=disp_pp(k)
    END DO
    CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
  END DO; IF(numpe==1) CLOSE(12)
  IF(numpe==1) WRITE(11,'(A,F10.4)') "This analysis took  :",                  &
                                     elap_time()-timest(1)  
  
  !-----------------------------------------------------------------------------
  ! 11. Shut down PETSc and ParaEFM
  !-----------------------------------------------------------------------------
  
  CALL PetscFinalize(p_ierr)
  CALL SHUTDOWN() 
  
CONTAINS

  !/****if* petsc/p_describe_reason
  !*  NAME
  !*    SUBROUTINE: p_describe_reason
  !*  SYNOPSIS
  !*    Usage:      p_describe_reason(p_reason,p_description)
  !*  FUNCTION
  !*      Returns a description of the reason for convergence in PETSc
  !*  ARGUMENTS
  !*    INTENT(IN)
  !*
  !*    p_reason           : KSPConvergedReason
  !*                         The PETSc code for the convergence reason
  !*    INTENT(OUT)
  !*
  !*    p_description      : Character
  !*                         Description of the convergence reason
  !*  AUTHOR
  !*    Mark Filipiak
  !*  CREATION DATE
  !*    19.02.2016
  !*  MODIFICATION HISTORY
  !*    Version 1, 19.02.2016, Mark Filipiak
  !*  COPYRIGHT
  !*    (c) University of Edinburgh 2016
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  This will change with PETSc version.  The alternative is to use PETSc's
  !*  KSPConvergedReasons C string array but that only gives the name, not the
  !*  extended description.
  !*/
  SUBROUTINE p_describe_reason(p_reason,p_description)
    IMPLICIT NONE
    KSPConvergedReason, INTENT(IN)  :: p_reason
    CHARACTER(LEN=*),   INTENT(OUT) :: p_description

    ! The string constants in this routine have to be p_max_string_length or
    ! shorter.
    SELECT CASE (p_reason)
    ! Converged.
    ! not in PETSc Fortran CASE (KSP_CONVERGED_RTOL_NORMAL)
    ! not in PETSc Fortran   p_description = "KSP_CONVERGED_RTOL_NORMAL"
    ! not in PETSc Fortran CASE (KSP_CONVERGED_ATOL_NORMAL)
    ! not in PETSc Fortran   p_description = "KSP_CONVERGED_ATOL_NORMAL"
    CASE (KSP_CONVERGED_RTOL)
      p_description = "KSP_CONVERGED_RTOL "                                    &
                      //"(residual norm <= rtol * initial residual norm.  "    &
                      //"The default norm for left preconditioning is the "    &
                      //"2-norm of the preconditioned residual.  "             &
                      //"The default norm for right preconditioning is the "   &
                      //"2-norm of the residual.)"
    CASE (KSP_CONVERGED_ATOL)
      p_description = "KSP_CONVERGED_ATOL (residual norm <= abstol.  "         &
                      //"The default norm for left preconditioning is the "    &
                      //"2-norm of the preconditioned residual.  "             &
                      //"The default norm for right preconditioning is the "   &
                      //"2-norm of the residual.)"
    CASE (KSP_CONVERGED_ITS)
      p_description = "KSP_CONVERGED_ITS "                                     &
                      //"(Used by the KSPPREONLY solver after the single "     &
                      //"iteration of the preconditioner is applied.  Also "   &
                      //"used when the KSPConvergedSkip() convergence test "   &
                      //"routine is set in KSP.)"
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
      p_description = "KSP_DIVERGED_ITS "                                      &
                      //"(Ran out of iterations before any convergence "       &
                      //"criteria was reached"
    CASE (KSP_DIVERGED_DTOL)
      p_description = "KSP_DIVERGED_DTOL "                                     &
                      //"(residual norm >= dtol * initial residual norm.  "    &
                      //"The default norm for left preconditioning is the "    &
                      //"2-norm of the preconditioned residual.  "             &
                      //"The default norm for right preconditioning is the "   &
                      //"2-norm of the residual.)"
    CASE (KSP_DIVERGED_BREAKDOWN)
      p_description = "KSP_DIVERGED_BREAKDOWN "                                &
                      //"(A breakdown in the Krylov method was detected so "   &
                      //"the method could not continue to enlarge the Krylov " &
                      //"space. Could be due to a singular matrix or "         &
                      //"preconditioner.)"
    CASE (KSP_DIVERGED_BREAKDOWN_BICG)
      p_description = "KSP_DIVERGED_BREAKDOWN_BICG "                           &
                      //"(A breakdown in the KSPBICG method was detected so "  &
                      //"the method could not continue to enlarge the Krylov " &
                      //"space.)"
    CASE (KSP_DIVERGED_NONSYMMETRIC)
      p_description = "KSP_DIVERGED_NONSYMMETRIC "                             &
                      //"(It appears the operator or preconditioner is not "   &
                      //"symmetric and this Krylov method (KSPCG, KSPMINRES, " &
                      //"KSPCR) requires symmetry.)"
    CASE (KSP_DIVERGED_INDEFINITE_PC)
      p_description = "KSP_DIVERGED_INDEFINITE_PC "                            &
                      //"(It appears the preconditioner is indefinite (has "   &
                      //"both positive and negative eigenvalues) and this "    &
                      //"Krylov method (KSPCG) requires it to be positive "    &
                      //"definite.  This can happen with the PCICC "           &
                      //"preconditioner; use "                                 &
                      //"-pc_factor_shift_positive_definite to force the "     &
                      //"PCICC preconditioner to generate a positive definite "&
                      //"preconditioner.)"
    CASE (KSP_DIVERGED_NANORINF)
      p_description = "KSP_DIVERGED_NANORINF "                                 &
                      //"(residual norm became NaN or Inf likely due to 0/0.)"
    CASE (KSP_DIVERGED_INDEFINITE_MAT)
      p_description = "KSP_DIVERGED_INDEFINITE_MAT"
    ! not in PETSc Fortran CASE (KSP_DIVERGED_PCSETUP_FAILED)
    ! not in PETSc Fortran   p_description = "KSP_DIVERGED_PCSETUP_FAILED"
    ! Still iterating.
    CASE (KSP_CONVERGED_ITERATING)
      p_description = "KSP_CONVERGED_ITERATING "                               &
                      //"(This flag is returned if you call "                  &
                      //"KSPGetConvergedReason() while KSPSolve() is still "   &
                      //"running.)"
    ! Unknown  
    CASE DEFAULT
      p_description = "reason not known"
    END SELECT
  END SUBROUTINE p_describe_reason
  
  !/****if* petsc/p_row_nnz
  !*  NAME
  !*    SUBROUTINE: p_row_nnz
  !*  SYNOPSIS
  !*    Usage:      p_row_nnz(nodof,ndim,nod,over_allocation,nnz)
  !*  FUNCTION
  !*    Returns an approximate upper estimate of the number of non-zeroes in a
  !*    row of the global matrix.  The number of non-zeroes in a row is used to
  !*    get an approximate pre-allocation for the PETSc matrix assembly.  Too
  !*    small a value will slow down the assembly, too large a value will waste
  !*    memory. 
  !*  ARGUMENTS
  !*    INTENT(IN)
  !*
  !*    nodof              : Integer
  !*                         Number of dofs per node
  !*    ndim               : Integer
  !*                         Number of dimensions of the problem
  !*    nod                : Integer
  !*                         Number of nodes per element
  !*    over_allocation    : Real
  !*                         over_allocation is a safety factor and allows for
  !*                         distorted grids with more than the simple number of
  !*                         elements around a point or edge.  Chose this to be
  !*                         larger than 1.25, which would correspond to 5
  !*                         hexahedra about an edge.  (Note that the average
  !*                         number of tetrahedra around a point is about 22, so
  !*                         the safety-factor will allow for that as well.
  !*
  !*    INTENT(OUT)
  !*
  !*    nnz                : PetscInt
  !*                         Approximate upper estimate of the number of
  !*                         non-zeroes in a row of the global matrix
  !*  AUTHOR
  !*    Mark Filipiak
  !*  CREATION DATE
  !*    19.02.2016
  !*  MODIFICATION HISTORY
  !*    Version 1, 19.02.2016, Mark Filipiak
  !*  COPYRIGHT
  !*    (c) University of Edinburgh 2016
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  This will not work if there are several different elements in the mesh
  !*  (or would have to called for each element and the largest nnz used).
  !*/
  SUBROUTINE p_row_nnz(nodof,ndim,nod,over_allocation,nnz)
    IMPLICIT NONE
    INTEGER,  INTENT(IN)   :: nodof
    INTEGER,  INTENT(IN)   :: ndim
    INTEGER,  INTENT(IN)   :: nod
    REAL,     INTENT(IN)   :: over_allocation
    PetscInt, INTENT(OUT)  :: nnz

    INTEGER                :: v,e,f,i,el_per_v,el_per_e,el_per_f,el_per_i

    SELECT CASE (ndim)

    CASE(1) ! one dimensional case
      SELECT CASE (nod)
      CASE(2)
        el_per_v = 2; el_per_i = 1 ! line
      ! two lines     line
        ! two lines
        !  1 interior vertices
        !  2 exterior vertices
        !  2 interior lines
        ! line
        !  0 interior vertices
        !  2 exterior vertices
        !  1 interior lines
        v = 3; i = 0
      CASE(3)
        el_per_v = 2; el_per_i = 1
        v = 5; i = 3
      CASE(4)
        el_per_v = 2; el_per_i = 1
        v = 7; i = 4
      CASE(5)
        el_per_v = 2; el_per_i = 1
        v = 9; i = 5
      CASE DEFAULT
        WRITE(*,'(A)') "wrong number of nodes in p_row_nnz"        
      END SELECT
    CASE(2)      ! two dimensional elements
      SELECT CASE (nod)
      CASE(3)
        el_per_v = 6; el_per_e = 2; el_per_i = 1 ! triangle
      ! hexagon       rhombus       triangle
        ! hexagon
        !  1 interior vertices
        !  6 exterior vertices
        !  6 interior edges
        !  6 exterior edges
        !  6 interior triangles
        ! rhombus
        !  0 interior vertices
        !  4 exterior vertices
        !  1 interior edges
        !  4 exterior edges
        !  2 interior triangles
        ! triangle
        !  0 interior vertices
        !  3 exterior vertices
        !  0 interior edges
        !  3 exterior edges
        !  1 interior triangles
        v = 7; e = 0; i = 0
      CASE(6) 
        el_per_v = 6; el_per_e = 2; el_per_i = 1 ! triangle
        v = 19; e = 9; i = 0
      CASE(10) 
        el_per_v = 6; el_per_e = 2; el_per_i = 1 ! triangle
        v = 37; e = 16; i = 10
      CASE(15)  
        el_per_v = 6; el_per_e = 2; el_per_i = 1 ! triangle
        v = 61; e = 25; i = 15
      CASE(4)  
        el_per_v = 4; el_per_e = 2; el_per_i = 1 ! quadrilateral
      ! 4 squares     2 squares     square
        ! 4 squares
        !  1 interior vertices
        !  8 exterior vertices
        !  4 interior edges
        !  8 exterior edges
        !  4 interior squares
        ! 2 squares
        !  0 interior vertices
        !  6 exterior vertices
        !  1 interior edges
        !  6 exterior edges
        !  2 interior squares
        ! square
        !  0 interior vertices
        !  4 exterior vertices
        !  0 interior edges
        !  4 exterior edges
        !  1 interior squares
        v = 9; e = 0; i = 0
      CASE(8)
        el_per_v = 4; el_per_e = 2; el_per_i = 1 ! quadrilateral
        v = 21; e = 13; i = 0
      CASE(9)
        el_per_v = 4; el_per_e = 2; el_per_i = 1 ! quadrilateral
        v = 25; e = 15; i = 9
      CASE DEFAULT
        WRITE(*,'(A)') "wrong number of nodes in p_row_nnz"        
      END SELECT
    CASE(3)  ! three dimensional elements
      SELECT CASE (nod)
      CASE(4)
        el_per_v = 20; el_per_e = 5; el_per_f = 2; el_per_i = 1 ! tetrahedron
      ! icosahedron    pentagonal    triangular    tetrahedron
      !                bipyramid     bipyramid
        ! icosahedron
        !  1 interior vertices
        ! 12 exterior vertices
        ! 12 interior edges
        ! 30 exterior edges
        ! 30 interior faces
        ! 20 exterior faces
        ! 20 interior tetrahedra
        ! pentagonal bipyramid
        !  0 interior vertices
        !  7 exterior vertices
        !  1 interior edges
        ! 15 exterior edges
        !  5 interior faces
        ! 10 exterior faces
        !  5 interior tetrahedra
        ! triangular bipyramid
        !  0 interior vertices
        !  5 exterior vertices
        !  0 interior edges
        !  9 exterior edges
        !  1 interior faces
        !  6 exterior faces
        !  2 interior tetrahedra
        ! tetrahedron
        !  0 interior vertices
        !  4 exterior vertices
        !  0 interior edges
        !  6 exterior edges
        !  0 interior faces
        !  4 exterior faces
        !  1 interior tetrahedra
        v = 13; e = 0; f = 0; i = 0
      CASE(8)
        el_per_v = 8; el_per_e = 4; el_per_f = 2; el_per_i = 1 ! hexahedron
      ! 8 cubes       4 cubes       2 cubes       cube
        ! 8 cubes
        !  1 interior vertices
        ! 26 exterior vertices
        !  6 interior edges
        ! 48 exterior edges
        ! 12 interior faces
        ! 24 exterior faces
        !  8 interior cubes
        ! 4 cubes
        !  0 interior vertices
        ! 18 exterior vertices
        !  1 interior edges
        ! 32 exterior edges
        !  4 interior faces
        ! 16 exterior faces
        !  4 interior cubes
        ! 2 cubes
        !  0 interior vertices
        ! 12 exterior vertices
        !  0 interior edges
        ! 20 exterior edges
        !  1 interior faces
        ! 10 exterior faces
        !  2 interior cubes
        ! cube
        !  0 interior vertices
        !  8 exterior vertices
        !  0 interior edges
        ! 12 exterior edges
        !  0 interior faces
        !  6 exterior faces
        !  1 interior cubes
        v = 27; e = 0; f = 0; i = 0
      CASE(14) ! type 6 element
        el_per_v = 8; el_per_e = 4; el_per_f = 2; el_per_i = 1 ! hexahedron
        v = 63; e = 0; f = 23; i = 0
      CASE(20)
        el_per_v = 8; el_per_e = 4; el_per_f = 2; el_per_i = 1 ! hexahedron
        v = 81; e = 50; f = 0; i = 0
      CASE DEFAULT
        WRITE(*,'(A)') "wrong number of nodes in p_row_nnz"        
      END SELECT
    CASE DEFAULT
      WRITE(*,'(A)') "wrong number of dimensions in p_row_nnz"
    END SELECT
    
    nnz = NINT(over_allocation * v) * nodof
  END SUBROUTINE p_row_nnz

END PROGRAM xx18
  
