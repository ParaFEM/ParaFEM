MODULE parafem_petsc
  
  !/****h* /petsc
  !*  NAME
  !*    MODULE: petsc
  !*  SYNOPSIS
  !*    Usage:      USE petsc
  !*  FUNCTION
  !*    Contains data and subroutines that define the PETSc interface.  These
  !*    subroutines are parallel and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    p_describe_reason      Describes reason for solver convergence
  !*    p_row_nnz              Estimate the number of non-zeroes in a row
  !*  AUTHOR
  !*    Mark Filipiak
  !*  COPYRIGHT
  !*    2016 University of Edinburgh
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  If you contribute to this module, add your author name.
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

  IMPLICIT NONE

  PUBLIC

  ! PETSc types
#include <petsc/finclude/petscdef.h>
  
  ! PETSc variables
  INTEGER,PARAMETER   :: p_max_string_length=1024
  ! p_over_allocation should be a run time setting with a default set in the
  ! program.  Choosing too small a value (e.g. 0.99) slows MatSetValues to a
  ! crawl.
  REAL,PARAMETER      :: p_over_allocation=1.3

  ! the following will become a global data structure variable (not p_ierr)

!!$  PetscErrorCode      :: p_ierr
!!$  ! The PETSc objects cannot be initialised here because PETSC_NULL_OBJECT is a
!!$  ! common-block-object and not a constant.
!!$  Vec                 :: p_x,p_b,p_r
!!$  Mat                 :: p_A
!!$  KSP                 :: p_ksp
!!$  PetscScalar         :: p_pr_n2,p_r_n2,p_b_n2
!!$  PetscInt,ALLOCATABLE    :: p_rows(:),p_cols(:)
!!$  PetscScalar,ALLOCATABLE :: p_values(:)
!!$  PetscScalar,POINTER :: p_varray(:)
!!$  PetscInt            :: p_its,p_nnz
!!$  DOUBLE PRECISION    :: p_info(MAT_INFO_SIZE)
!!$  KSPConvergedReason  :: p_reason
!!$  CHARACTER(LEN=p_max_string_length) :: p_description
  
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

END MODULE parafem_petsc
