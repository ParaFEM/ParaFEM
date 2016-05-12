MODULE parafem_petsc
  
  !/****h* /parafem_petsc
  !*  NAME
  !*    MODULE: parafem_petsc
  !*  SYNOPSIS
  !*    Usage:      USE parafem_petsc
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
  
  ! PETSc will always use MPI (unless PETSc itself has been compiled without
  ! MPI)
  USE mp_interface

  IMPLICIT NONE

  PUBLIC

  ! PETSc types
#include <petsc/finclude/petscdef.h>

  ! PETSc modules could be used, but the ParaFEM mp_interface module INCLUDEs
  ! mpif.h rather than USEing mpi, so there are multiple definitions from the
  ! USE mpi that is in the PETSc modules.  So include the PETSc headers and
  ! don't include mpif.h again.
#define PETSC_AVOID_MPIF_H
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscmat.h90>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscksp.h90>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscpc.h90>
  
  ! Private parameters
  INTEGER,PARAMETER,PRIVATE :: string_length=1024
  
  ! Variables collected as one derived type
  TYPE p_type
    ! The PETSc objects cannot be initialised here because PETSC_NULL_OBJECT is
    ! a common-block-object and not a constant.
    ! over_allocation should be a run time setting with a default set in the
    ! program, e.g. as an option.  Choosing too small a value (e.g. 0.99) slows
    ! MatSetValues to a crawl.
    REAL                :: over_allocation=1.3
    Vec                 :: x,b,r
    Mat                 :: A
    KSP                 :: ksp
    PetscScalar         :: pr_n2,r_n2,b_n2
    PetscInt,ALLOCATABLE    :: rows(:),cols(:)
    PetscScalar,ALLOCATABLE :: values(:)
    PetscScalar,POINTER :: varray(:)
    PetscInt            :: its,nnz
    PetscReal           :: rtol
    DOUBLE PRECISION    :: info(MAT_INFO_SIZE)
    KSPConvergedReason  :: reason
    CHARACTER(len=string_length) :: description
    PetscErrorCode      :: ierr
  END TYPE p_type

  TYPE(p_type) :: p_object
 
CONTAINS
  
  SUBROUTINE p_initialize(numpe,fname_base)

    !/****if* petsc/p_initialize
    !*  NAME
    !*    SUBROUTINE: p_initialize
    !*  SYNOPSIS
    !*    Usage:      p_initialize(numpe,argv)
    !*  FUNCTION
    !*      Initialises PETSc and its matrices, vectors and solvers
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    numpe              : Integer
    !*                         Number of this process (starting at 1)
    !*    fname_base         : Character
    !*                         Base name of the data file
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
    !*/

    INTEGER,          INTENT(in) :: numpe
    CHARACTER(len=*), INTENT(in) :: fname_base

    CHARACTER(len=string_length) :: fname
    LOGICAL :: exist
    INTEGER :: ierr

    IF(numpe == 1) THEN
      fname = TRIM(fname_base) // ".petsc"
      INQUIRE(file=TRIM(fname),exist=exist)
      IF (.NOT. exist) THEN
        fname = ""
      END IF
    END IF
    ! The communicator is MPI_COMM_WORLD, set by find_pe_procs(), and numpe ==
    ! 1 corresponds to rank == 0.
    CALL MPI_BCAST(fname,LEN(fname),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    
    CALL PetscInitialize(fname,p_object%ierr)
    p_object%x   = PETSC_NULL_OBJECT
    p_object%b   = PETSC_NULL_OBJECT
    p_object%r   = PETSC_NULL_OBJECT
    p_object%A   = PETSC_NULL_OBJECT
    p_object%ksp = PETSC_NULL_OBJECT
  END SUBROUTINE p_initialize

  SUBROUTINE p_create_matrix(neq_pp,nodof,ndim,nod,ntot)

    !/****if* petsc/p_create_matrix
    !*  NAME
    !*    SUBROUTINE: p_create_matrix
    !*  SYNOPSIS
    !*    Usage:      p_create_matrix(neq_pp,nodof,ndim,nod,ntot)
    !*  FUNCTION
    !*      Creates the global matrix (but doesn't assemble it)
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    neq_pp             : Integer
    !*                         Number of equations on this process = number of
    !*                         rows in the global matrix on this process
    !*    nodof              : Integer
    !*                         Number of dofs per node
    !*    ndim               : Integer
    !*                         Number of dimensions of the problem
    !*    nod                : Integer
    !*                         Number of nodes per element
    !*    ntot_max           : Integer
    !*                         Maximum number of dofs in an element (this will
    !*                         be ntot if there is only one element type).  This
    !*                         is used to set the sizes of some workspaces.
    !*                         
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.04.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 28.04.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  PETSc 64-bit indices and 64-bit reals are used.  In most (all?) places,
    !*  passing a 32-bit integer where an intent(in) 64-integer is required is
    !*  safe, because the PETSc Fortran-C interface de-references the pointer
    !*  that is actually passed, even though it does not specify any intent in
    !*  the Fortran interface.  This is not safe when passing arrays, they need
    !*  to be copied to PetscInt arrays.  And for safety the same should be done
    !*  for the PetscScalar arrays.
    !*/

    INTEGER :: nnz

    CALL MatCreate(PETSC_COMM_WORLD,p_object%A,p_object%ierr)
    CALL MatSetSizes(p_object%A,neq_pp,neq_pp,PETSC_DETERMINE,PETSC_DETERMINE, &
                     p_object%ierr)
    CALL MatSetType(p_object%A,MATAIJ,p_object%ierr)
    !- Block size fixed to 1 for just now - this cannot be set to nodof until
    !- the restraints are handled block-wise in ParaFEM.
    !- CALL MatSetBlockSize(p_object%A,nodof,p_object%ierr)
    
    ! Find an approximate number of zeroes per row for the matrix size
    ! pre-allocation.
    nnz = p_row_nnz(nodof,ndim,nod,p_object%over_allocation)
    CALL MatSeqAIJSetPreallocation(p_object%A,nnz,PETSC_NULL_INTEGER,          &
                                   p_object%ierr)
    CALL MatMPIAIJSetPreallocation(p_object%A,nnz,PETSC_NULL_INTEGER,          &
                                   nnz,PETSC_NULL_INTEGER,p_object%ierr)
    ! If the allocation is too small, PETSc will produce reams of information
    ! and not assemble the matrix properly.  We output some information at the
    ! end of the assembly routine if p_object%over_allocation should be
    ! increased.
    CALL MatSetOption(p_object%A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,   &
                      p_object%ierr)

    ! Arrays of PETSc types to be proof against changes in index and scalar
    ! sizes.
    ALLOCATE(p_object%rows(ntot),p_object%cols(ntot),p_object%values(ntot*ntot))

    CALL MatZeroEntries(p_object%A,p_object%ierr)
  END SUBROUTINE p_create_matrix

  SUBROUTINE p_add_element(g,km)

    !/****if* petsc/p_add_element
    !*  NAME
    !*    SUBROUTINE: p_add_element
    !*  SYNOPSIS
    !*    Usage:      p_add_element(g,km)
    !*  FUNCTION
    !*    Adds (assembles) an element matrix into the global matrix.  Once all
    !*    the element matrices are added on each process, the final assembly
    !*    (what PETSc calls assembly) is done to get the global, distributed
    !*    matrix.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    g(:)               : Integer
    !*                         Array of shape (ntot)
    !*                         Gives the global equation number for each dof in
    !*                         the element. 
    !*    km(:,:)            : Real
    !*                         Array of shape (ntot,ntot)
    !*                         The element matrix.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.04.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 28.04.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  This depends on g holding 0 for erased rows/columns (i.e. equations),
    !*  and MatSetValues ignoring negative indices (PETSc always uses zero-based
    !*  indexing).
    !*
    !*  This also depends on the element matrices being shape (ntot,ntot),
    !*  because rows, cols, values are (ntot), (ntot), (ntot,ntot).  What about
    !differing elelemtns?
    !*/
    
    p_object%rows   = g - 1
    p_object%cols   = p_rows
    ! PETSc uses C array order, so a transpose is needed.
    p_object%values = RESHAPE(TRANSPOSE(km),SHAPE(p_object%values))
    CALL MatSetValues(p_object%A,SIZE(p_object%rows),p_object%rows,SIZE(p_object%cols),p_object%cols,p_object%values,ADD_VALUES,p_object%ierr)

    INTEGER :: nnz

    CALL MatCreate(PETSC_COMM_WORLD,p_object%A,p_object%ierr)
    CALL MatSetSizes(p_object%A,neq_pp,neq_pp,PETSC_DETERMINE,PETSC_DETERMINE, &
                     p_object%ierr)
    CALL MatSetType(p_object%A,MATAIJ,p_object%ierr)
    !- Block size fixed to 1 for just now - this cannot be set to nodof until
    !- the restraints are handled block-wise in ParaFEM.
    !- CALL MatSetBlockSize(p_object%A,nodof,p_object%ierr)
    
    ! Find an approximate number of zeroes per row for the matrix size
    ! pre-allocation.
    nnz = p_row_nnz(nodof,ndim,nod,p_object%over_allocation)
    CALL MatSeqAIJSetPreallocation(p_object%A,nnz,PETSC_NULL_INTEGER,          &
                                   p_object%ierr)
    CALL MatMPIAIJSetPreallocation(p_object%A,nnz,PETSC_NULL_INTEGER,          &
                                   nnz,PETSC_NULL_INTEGER,p_object%ierr)
    ! If the allocation is too small, PETSc will produce reams of information
    ! and not assemble the matrix properly.  We output some information at the
    ! end of the assembly routine if p_object%over_allocation should be
    ! increased.
    CALL MatSetOption(p_object%A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,   &
                      p_object%ierr)

    ! Arrays of PETSc types to be proof against changes in index and scalar
    ! sizes.
    ALLOCATE(p_object%rows(ntot),p_object%cols(ntot),p_object%values(ntot*ntot))

    CALL MatZeroEntries(p_object%A,p_object%ierr)
  END SUBROUTINE p_create_matrix

  FUNCTION p_describe_reason(reason,description)

    !/****if* petsc/p_describe_reason
    !*  NAME
    !*    SUBROUTINE: p_describe_reason
    !*  SYNOPSIS
    !*    Usage:      p_describe_reason(p_reason)
    !*  FUNCTION
    !*      Returns a description of the reason for convergence in PETSc
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    p_reason           : KSPConvergedReason
    !*                         The PETSc code for the convergence reason
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

    CHARACTER(:),ALLOCATABLE        :: p_describe_reason
    KSPConvergedReason, INTENT(in)  :: reason
    
    ! The string constants in this routine have to be p_string_length or
    ! shorter.
    SELECT CASE (p_reason)
    ! Converged.
    ! not in PETSc Fortran CASE (KSP_CONVERGED_RTOL_NORMAL)
    ! not in PETSc Fortran   p_describe_reason = "KSP_CONVERGED_RTOL_NORMAL"
    ! not in PETSc Fortran CASE (KSP_CONVERGED_ATOL_NORMAL)
    ! not in PETSc Fortran   p_describe_reason = "KSP_CONVERGED_ATOL_NORMAL"
    CASE (KSP_CONVERGED_RTOL)
      p_describe_reason = "KSP_CONVERGED_RTOL "                                &
                      //"(residual norm <= rtol * initial residual norm.  "    &
                      //"The default norm for left preconditioning is the "    &
                      //"2-norm of the preconditioned residual.  "             &
                      //"The default norm for right preconditioning is the "   &
                      //"2-norm of the residual.)"
    CASE (KSP_CONVERGED_ATOL)
      p_describe_reason = "KSP_CONVERGED_ATOL (residual norm <= abstol.  "     &
                      //"The default norm for left preconditioning is the "    &
                      //"2-norm of the preconditioned residual.  "             &
                      //"The default norm for right preconditioning is the "   &
                      //"2-norm of the residual.)"
    CASE (KSP_CONVERGED_ITS)
      p_describe_reason = "KSP_CONVERGED_ITS "                                 &
                      //"(Used by the KSPPREONLY solver after the single "     &
                      //"iteration of the preconditioner is applied.  Also "   &
                      //"used when the KSPConvergedSkip() convergence test "   &
                      //"routine is set in KSP.)"
    CASE (KSP_CONVERGED_CG_NEG_CURVE)
      p_describe_reason = "KSP_CONVERGED_CG_NEG_CURVE"
    CASE (KSP_CONVERGED_CG_CONSTRAINED)
      p_describe_reason = "KSP_CONVERGED_CG_CONSTRAINED"
    CASE (KSP_CONVERGED_STEP_LENGTH)
      p_describe_reason = "KSP_CONVERGED_STEP_LENGTH"
    CASE (KSP_CONVERGED_HAPPY_BREAKDOWN)
      p_describe_reason = "KSP_CONVERGED_HAPPY_BREAKDOWN"
    ! Diverged.
    CASE (KSP_DIVERGED_NULL)
      p_describe_reason = "KSP_DIVERGED_NULL"
    CASE (KSP_DIVERGED_ITS)
      p_describe_reason = "KSP_DIVERGED_ITS "                                  &
                      //"(Ran out of iterations before any convergence "       &
                      //"criteria was reached"
    CASE (KSP_DIVERGED_DTOL)
      p_describe_reason = "KSP_DIVERGED_DTOL "                                 &
                      //"(residual norm >= dtol * initial residual norm.  "    &
                      //"The default norm for left preconditioning is the "    &
                      //"2-norm of the preconditioned residual.  "             &
                      //"The default norm for right preconditioning is the "   &
                      //"2-norm of the residual.)"
    CASE (KSP_DIVERGED_BREAKDOWN)
      p_describe_reason = "KSP_DIVERGED_BREAKDOWN "                            &
                      //"(A breakdown in the Krylov method was detected so "   &
                      //"the method could not continue to enlarge the Krylov " &
                      //"space. Could be due to a singular matrix or "         &
                      //"preconditioner.)"
    CASE (KSP_DIVERGED_BREAKDOWN_BICG)
      p_describe_reason = "KSP_DIVERGED_BREAKDOWN_BICG "                       &
                      //"(A breakdown in the KSPBICG method was detected so "  &
                      //"the method could not continue to enlarge the Krylov " &
                      //"space.)"
    CASE (KSP_DIVERGED_NONSYMMETRIC)
      p_describe_reason = "KSP_DIVERGED_NONSYMMETRIC "                         &
                      //"(It appears the operator or preconditioner is not "   &
                      //"symmetric and this Krylov method (KSPCG, KSPMINRES, " &
                      //"KSPCR) requires symmetry.)"
    CASE (KSP_DIVERGED_INDEFINITE_PC)
      p_describe_reason = "KSP_DIVERGED_INDEFINITE_PC "                        &
                      //"(It appears the preconditioner is indefinite (has "   &
                      //"both positive and negative eigenvalues) and this "    &
                      //"Krylov method (KSPCG) requires it to be positive "    &
                      //"definite.  This can happen with the PCICC "           &
                      //"preconditioner; use "                                 &
                      //"-pc_factor_shift_positive_definite to force the "     &
                      //"PCICC preconditioner to generate a positive definite "&
                      //"preconditioner.)"
    CASE (KSP_DIVERGED_NANORINF)
      p_describe_reason = "KSP_DIVERGED_NANORINF "                             &
                      //"(residual norm became NaN or Inf likely due to 0/0.)"
    CASE (KSP_DIVERGED_INDEFINITE_MAT)
      p_describe_reason = "KSP_DIVERGED_INDEFINITE_MAT"
    ! not in PETSc Fortran CASE (KSP_DIVERGED_PCSETUP_FAILED)
    ! not in PETSc Fortran   p_describe_reason = "KSP_DIVERGED_PCSETUP_FAILED"
    ! Still iterating.
    CASE (KSP_CONVERGED_ITERATING)
      p_describe_reason = "KSP_CONVERGED_ITERATING "                           &
                      //"(This flag is returned if you call "                  &
                      //"KSPGetConvergedReason() while KSPSolve() is still "   &
                      //"running.)"
    ! Unknown  
    CASE DEFAULT
      p_describe_reason = "reason not known"
    END SELECT
  END FUNCTION p_describe_reason
  
  FUNCTION p_row_nnz(nodof,ndim,nod,over_allocation)

    !/****if* petsc/p_row_nnz
    !*  NAME
    !*    SUBROUTINE: p_row_nnz
    !*  SYNOPSIS
    !*    Usage:      p_row_nnz(nodof,ndim,nod,over_allocation,nnz)
    !*  FUNCTION
    !*    Returns an approximate upper estimate of the number of non-zeroes in a
    !*    row of the global matrix.  The number of non-zeroes in a row is used
    !*    to get an approximate pre-allocation for the PETSc matrix assembly.
    !*    Too small a value will slow down the assembly, too large a value will
    !*    waste memory.
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
    !*                         distorted grids with more than the simple number
    !*                         of elements around a point or edge.  Chose this
    !*                         to be larger than 1.25, which would correspond to
    !*                         5 hexahedra about an edge.  (Note that the
    !*                         average number of tetrahedra around a point is
    !*                         about 22, so the safety-factor will allow for
    !*                         that as well.
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

    PetscInt               :: nnz
    INTEGER,  INTENT(in)   :: nodof
    INTEGER,  INTENT(in)   :: ndim
    INTEGER,  INTENT(in)   :: nod
    REAL,     INTENT(in)   :: over_allocation

    INTEGER                :: el_per_v,el_per_e,el_per_f,el_per_i
    REAL                   :: v,e,f,i,v_factor,e_factor

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
        ! Note that the average number of tetrahedra per vertex is about 22.79
        ! based on averaging the solid angles and the average number of
        ! tetrahedra per edge is 5.1 based on averaging the dihedral angles
        v_factor=22.70/20; e_factor = 5.1/5
        ! Note also that there can be many vertices with large numbers of
        ! tetrahedra depending on how the mesh is constructed, for example if it
        ! has been refined, see A. Plaza and M-C. Rivera, 'On the adjacencies of
        ! triangular meshes based on skeleton-regular partitions', Journal of
        ! Computational and Applied Mathematics 140 (2002) 673-693
        ! http://dx.doi.org/10.1016/S0377-0427(01)00484-8.  So there can be
        ! possibly several percent of rows in the global matrix with too small a
        ! pre-allocation and thus slow assembly in PETSc.  See also Table X of
        ! M.W. Beall and M.S. Shephard, 'A general topology-based mesh data
        ! structure', Internatioanl Journal for Numerical Methods in Engineering
        ! 40 (1997) 1573-1596 for average node-node adjacencies - these closely
        ! match the values calculated below.
        !
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
        v = 13 * v_factor; e = 0 * e_factor; f = 0; i = 0
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
        v = 81; e = 51; f = 0; i = 0
      CASE DEFAULT
        WRITE(*,'(A)') "wrong number of nodes in p_row_nnz"        
      END SELECT
    CASE DEFAULT
      WRITE(*,'(A)') "wrong number of dimensions in p_row_nnz"
    END SELECT
    
    p_row_nnz = NINT(over_allocation * v) * nodof
  END FUNCTION p_row_nnz

END MODULE parafem_petsc
