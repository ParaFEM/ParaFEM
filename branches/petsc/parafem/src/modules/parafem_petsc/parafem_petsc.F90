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
    PetscInt            :: row_nnz=-1
    ! over_allocation should be a run time setting with a default set in the
    ! program, e.g. as an option.  Choosing too small a value (e.g. 0.99) slows
    ! MatSetValues to a crawl.
    REAL                :: over_allocation=1.3
    Vec                 :: x,b
    Mat                 :: A
    KSP                 :: ksp
    PetscInt,ALLOCATABLE    :: rows(:),cols(:)
    PetscScalar,ALLOCATABLE :: values(:)
    PetscInt            :: its
    PetscReal           :: rtol,r_n2,b_n2
    DOUBLE PRECISION    :: info(MAT_INFO_SIZE)
    KSPConvergedReason  :: reason
    CHARACTER(len=string_length) :: description=""
    PetscErrorCode      :: ierr
  END TYPE p_type

  TYPE(p_type) :: p_object

  ! Private functions
  PRIVATE :: nnod
  PRIVATE :: row_nnz_2
  PRIVATE :: row_nnz_1
  PRIVATE :: row_nnz
 
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
    p_object%A   = PETSC_NULL_OBJECT
    p_object%ksp = PETSC_NULL_OBJECT
  END SUBROUTINE p_initialize

  SUBROUTINE p_finalize()

    !/****if* petsc/p_finalize
    !*  NAME
    !*    SUBROUTINE: p_finalize
    !*  SYNOPSIS
    !*    Usage:      p_finalize()
    !*  FUNCTION
    !*      Destroys matrices, vectors, solvers and workspace, then finalizes
    !*      (i.e., tidies up) PETSc.
    !*  ARGUMENTS
    !*    None
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    18.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 29.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  There should be a p_destroy routine for each p_create routine
    !*/

    DEALLOCATE(p_object%rows,p_object%cols,p_object%values)
    CALL KSPDestroy(p_object%ksp,p_object%ierr)
    CALL VecDestroy(p_object%x,p_object%ierr)
    CALL VecDestroy(p_object%b,p_object%ierr)
    CALL MatDestroy(p_object%A,p_object%ierr)
    CALL PetscFinalize(p_object%ierr)
  END SUBROUTINE p_finalize

  SUBROUTINE p_create_matrix(neq_pp)

    !/****if* petsc/p_create_matrix
    !*  NAME
    !*    SUBROUTINE: p_create_matrix
    !*  SYNOPSIS
    !*    Usage:      p_create_matrix(neq_pp)
    !*  FUNCTION
    !*      Creates the global matrix (but doesn't assemble it)
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    neq_pp             : Integer
    !*                         Number of equations on this process = number of
    !*                         rows in the global matrix on this process
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.04.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 28.05.2016, Mark Filipiak
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
    !*
    !*  Block size fixed to 1 for just now - this cannot be set to nodof until
    !*  the restraints are handled block-wise in ParaFEM.
    !*  What about multi-field?
    !*/

    INTEGER, INTENT(in) :: neq_pp

    CALL MatCreate(PETSC_COMM_WORLD,p_object%A,p_object%ierr)
    CALL MatSetSizes(p_object%A,neq_pp,neq_pp,PETSC_DETERMINE,PETSC_DETERMINE, &
                     p_object%ierr)
    CALL MatSetType(p_object%A,MATAIJ,p_object%ierr)
    ! CALL MatSetBlockSize(p_object%A,nodof,p_object%ierr)
    
    CALL MatSeqAIJSetPreallocation(p_object%A,                                 &
                                   p_object%row_nnz,PETSC_NULL_INTEGER,        &
                                   p_object%ierr)
    CALL MatMPIAIJSetPreallocation(p_object%A,                                 &
                                   p_object%row_nnz,PETSC_NULL_INTEGER,        &
                                   p_object%row_nnz,PETSC_NULL_INTEGER,        &
                                   p_object%ierr)
    ! If the allocation is too small, PETSc will produce reams of information
    ! and not assemble the matrix properly.  We output some information at the
    ! end of the assembly routine if p_object%over_allocation should be
    ! increased.
    CALL MatSetOption(p_object%A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,   &
                      p_object%ierr)

  END SUBROUTINE p_create_matrix

  SUBROUTINE p_create_workspace(ntot_max)

    !/****if* petsc/p_create_workspace
    !*  NAME
    !*    SUBROUTINE: p_create_workspace
    !*  SYNOPSIS
    !*    Usage:      p_create_workspace(ntot_max)
    !*  FUNCTION
    !*      Creates some workspace
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    ntot_max           : Integer
    !*                         Maximum number of dofs in an element.  This will
    !*                         be ntot if there is only one element type.  If
    !*                         there are various types of element (including
    !*                         variations between different fields if there are
    !*                         several fields), this will be the maximum over
    !*                         these.  If the elements are combined to give one
    !*                         matrix, as in p126, then this will be ntot of
    !*                         this 'element'.  This is used to set the sizes of
    !*                         the workspaces.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 28.05.2016, Mark Filipiak
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

    INTEGER, INTENT(in) :: ntot_max

    ! Arrays of PETSc types to be proof against changes in index and scalar
    ! sizes.
    ALLOCATE(p_object%rows(ntot_max),p_object%cols(ntot_max),                  &
             p_object%values(ntot_max*ntot_max))
  END SUBROUTINE p_create_workspace

  SUBROUTINE p_zero_matrix()

    !/****if* petsc/p_zero_matrix
    !*  NAME
    !*    SUBROUTINE: p_zero_matrix
    !*  SYNOPSIS
    !*    Usage:      p_zero_matrix()
    !*  FUNCTION
    !*      Zeroes the global matrix
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 28.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    CALL MatZeroEntries(p_object%A,p_object%ierr)
  END SUBROUTINE p_zero_matrix

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
    !*    16.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 16.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  This depends on g holding 0 for erased rows/columns (i.e. equations),
    !*  and MatSetValues ignoring negative indices (PETSc always uses zero-based
    !*  indexing).
    !*/

    INTEGER, INTENT(in) :: g(:)
    REAL,    INTENT(in) :: km(:,:)

    PetscInt :: ntot
    PetscInt :: rows(60), cols(60)
    PetscScalar :: values (60*60)
    
    ntot = 60
    rows(1:ntot) = g - 1
    cols(1:ntot) = rows(1:ntot)
    ! PETSc uses C array order, so a transpose is needed.  Actually, you can
    ! get PETSc to use Fortran order for MatSetValues by setting
    ! MAT_ROW_ORIENTED false.  But other routines might be affected?
    CALL MatSetOption(p_object%A,MAT_ROW_ORIENTED,PETSC_FALSE,p_object%ierr)
    values(1:ntot*ntot) = RESHAPE(TRANSPOSE(km),(/ntot*ntot/))
    WRITE(*,*) ntot, SIZE(g), SIZE(rows), SIZE(cols)
    WRITE(*,*) ntot*ntot, SIZE(values), SIZE(km)
    CALL MatSetValues(p_object%A,ntot,rows,ntot,cols,        &
                      values,ADD_VALUES,p_object%ierr)

!!$    ntot = SIZE(g)
!!$    p_object%rows(1:ntot) = g - 1
!!$    p_object%cols(1:ntot) = p_object%rows(1:ntot)
!!$    ! PETSc uses C array order, so a transpose is needed.  Actually, you can
!!$    ! get PETSc to use Fortran order for MatSetValues by setting
!!$    ! MAT_ROW_ORIENTED false.  But other routines might be affected?
!!$    CALL MatSetOption(p_object%A,MAT_ROW_ORIENTED,PETSC_TRUE,p_object%ierr)
!!$    p_object%values(1:ntot*ntot) = RESHAPE(TRANSPOSE(km),(/ntot*ntot/))
!!$    CALL MatSetValues(p_object%A,ntot,p_object%rows,ntot,p_object%cols,        &
!!$                      p_object%values,ADD_VALUES,p_object%ierr)
  END SUBROUTINE p_add_element

  SUBROUTINE p_assemble(numpe)

    !/****if* petsc/p_assemble
    !*  NAME
    !*    SUBROUTINE: p_assemble
    !*  SYNOPSIS
    !*    Usage:      p_assemble(numpe)
    !*  FUNCTION
    !*    Assemble the global matrix.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    numpe              : Integer
    !*                         This processer's number.  Only processor 1
    !*                         prints information.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    17.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 17.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER, INTENT(in) :: numpe

    CALL MatAssemblyBegin(p_object%A,MAT_FINAL_ASSEMBLY,p_object%ierr)
    CALL MatAssemblyEnd(p_object%A,MAT_FINAL_ASSEMBLY,p_object%ierr)
    
    CALL MatGetInfo(p_object%A,MAT_GLOBAL_SUM,p_object%info,p_object%ierr)
    IF (numpe==1) THEN
      IF (p_object%info(MAT_INFO_MALLOCS)/=0.0) THEN
        WRITE(*,'(A,I0,A)') "The matrix assembly required ",                   &
                            NINT(p_object%info(MAT_INFO_MALLOCS)),             &
                            " mallocs.  Increase p_object%over_allocation to " &
                            //"speed up the assembly."
      END IF
    END IF
  END SUBROUTINE p_assemble

  SUBROUTINE p_create_vectors(neq_pp)

    !/****if* petsc/p_create_vectors
    !*  NAME
    !*    SUBROUTINE: p_create_vectors
    !*  SYNOPSIS
    !*    Usage:      p_create_vectors(neq_pp)
    !*  FUNCTION
    !*    Create the global RHS and solution vectors
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    neq_pp             : Integer
    !*                         Number of equations on this process = number of
    !*                         rows in the global matrix on this process
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 28.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  Block size fixed to 1 for just now - this cannot be set to nodof until
    !*  the restraints are handled block-wise in ParaFEM.
    !*  What about multi-field?
    !*/

    INTEGER, INTENT(in) :: neq_pp

    ! RHS vector.  For this particular order (create, set size, set type) the
    ! allocation is done by set type.
    CALL VecCreate(PETSC_COMM_WORLD,p_object%b,p_object%ierr)
    CALL VecSetSizes(p_object%b,neq_pp,PETSC_DECIDE,p_object%ierr)
    CALL VecSetType(p_object%b,VECSTANDARD,p_object%ierr)
    ! CALL VecSetBlockSize(p_object%b,nodof,p_object%ierr)

    ! Solution vector
    CALL VecDuplicate(p_object%b,p_object%x,p_object%ierr)
  END SUBROUTINE p_create_vectors
  
  SUBROUTINE p_create_ksp()

    !/****if* petsc/p_create_ksp
    !*  NAME
    !*    SUBROUTINE: p_create_ksp
    !*  SYNOPSIS
    !*    Usage:      p_create_ksp()
    !*  FUNCTION
    !*    Create the Krylov solver and preconditioner
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 28.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    CALL KSPCreate(PETSC_COMM_WORLD,p_object%ksp,p_object%ierr)
    CALL KSPSetOperators(p_object%ksp,p_object%A,p_object%A,p_object%ierr)
  END SUBROUTINE p_create_ksp
  
  SUBROUTINE p_set_tolerances(rtol,max_it)

    !/****if* petsc/p_set_tolerances
    !*  NAME
    !*    SUBROUTINE: p_set_tolerances
    !*  SYNOPSIS
    !*    Usage:      p_set_tolerances(rtol,max_it)
    !*  FUNCTION
    !*    Set the tolerances for the Krylov solver
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    rtol               : Real
    !*                         Relative tolerance for the preconditioned
    !*                         residual (in ParaFEM the relative tolerance
    !*                         is for the true residual).
    !*    max_it             : Integer
    !*                         Maximum number of Krylov solver iterations.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 28.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  KSP type, per-KSP tolerances (rtol, abstol, dtol, maxits), KSP
    !*  options, PC type, PC options are set in the xx*.petsc file.  Those
    !*  options are used to set up the preconditioned Krylov solver.  If there
    !*  are several KSP types to be chosen from, then each one will be
    !*  bracketed by -prefix_push and -prefix_pop.  For example
    !* 
    !*  -prefix_push abc1_
    !*    -ksp_type minres
    !*  -prefix_pop
    !* 
    !*  in the xx*.petsc file and 
    !* 
    !*  CALL KSPSetOptionsPrefix(p_object%ksp,"abc1_",p_object%ierr)
    !* 
    !*  before KSPSetFromOptions before using the 'abc1_' solver.  Thus you
    !*  can switch from CG to GMRES during a simulation.
    !*/

    REAL,    INTENT(in) :: rtol
    INTEGER, INTENT(in) :: max_it

    CALL KSPSetTolerances(p_object%ksp,rtol,PETSC_DEFAULT_REAL,                &
                          PETSC_DEFAULT_REAL,max_it,p_object%ierr)
    CALL KSPSetFromOptions(p_object%ksp,p_object%ierr)
  END SUBROUTINE p_set_tolerances
  
  SUBROUTINE p_set_rhs(r_pp)

    !/****if* petsc/p_set_rhs
    !*  NAME
    !*    SUBROUTINE: p_set_rhs
    !*  SYNOPSIS
    !*    Usage:      p_set_rhs(r_pp)
    !*  FUNCTION
    !*    Assemble the global RHS
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    r_pp(:)            : Real
    !*                         RHS on this process.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    17.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 17.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  No actual assembly is needed because the RHS vector has already
    !*  been distributed in ParaFEM.
    !*/

    REAL,    INTENT(in) :: r_pp(:)

    PetscScalar,POINTER :: varray(:)

    CALL VecGetArrayF90(p_object%b,varray,p_object%ierr)
    ! This is OK as long as PetscScalars are not smaller than ParaFEM reals.
    ! There should be a test for sizes of PetscScalars (and PetscInts) and
    ! ParaFEM reals (and indices).
    varray = r_pp
    CALL VecRestoreArrayF90(p_object%b,varray,p_object%ierr)
  END SUBROUTINE p_set_rhs

  SUBROUTINE p_set_solution(x_pp)

    !/****if* petsc/p_set_solution
    !*  NAME
    !*    SUBROUTINE: p_set_solution
    !*  SYNOPSIS
    !*    Usage:      p_set_solution(x_pp)
    !*  FUNCTION
    !*    Assemble the global solution
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    x_pp(:)            : Real
    !*                         Solution on this process.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    29.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 29.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  No actual assembly is needed because the solution vector has already
    !*  been distributed in ParaFEM.
    !*/

    REAL,    INTENT(in) :: x_pp(:)

    PetscScalar,POINTER :: varray(:)

    CALL VecGetArrayF90(p_object%x,varray,p_object%ierr)
    ! This is OK as long as PetscScalars are not smaller than ParaFEM reals.
    ! There should be a test for sizes of PetscScalars (and PetscInts) and
    ! ParaFEM reals (and indices).
    varray = x_pp
    CALL VecRestoreArrayF90(p_object%x,varray,p_object%ierr)
  END SUBROUTINE p_set_solution

  SUBROUTINE p_get_solution(x_pp)

    !/****if* petsc/p_get_solution
    !*  NAME
    !*    SUBROUTINE: p_get_solution
    !*  SYNOPSIS
    !*    Usage:      p_get_solution(x_pp)
    !*  FUNCTION
    !*    Convert the solution from PETSc structure to ParaFEM array
    !*  ARGUMENTS
    !*    INTENT(OUT)
    !*
    !*    x_pp(:)            : Real
    !*                         Solution on this process.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    29.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 29.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  No actual assembly is needed because the solution vector has already
    !*  been distributed in ParaFEM.
    !*/

    REAL,   INTENT(out) :: x_pp(:)

    PetscScalar,POINTER :: varray(:)

    CALL VecGetArrayF90(p_object%x,varray,p_object%ierr)
    ! This is OK as long as ParaFEM reals are not smaller than PetscScalars.
    ! There should be a test for sizes of PetscScalars (and PetscInts) and
    ! ParaFEM reals (and indices).
    x_pp = varray
    CALL VecRestoreArrayF90(p_object%x,varray,p_object%ierr)
  END SUBROUTINE p_get_solution

  SUBROUTINE p_solve(rtol,max_it,r_pp,x_pp,initial_guess_nonzero)

    !/****if* petsc/p_solve
    !*  NAME
    !*    SUBROUTINE: p_solve
    !*  SYNOPSIS
    !*    Usage:      p_solve(rtol,max_it,r_pp,x_pp,initial_guess_nonzero)
    !*  FUNCTION
    !*    Solve using PETSc.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    rtol               : Real
    !*                         Relative tolerance for the preconditioned
    !*                         residual (in ParaFEM the relative tolerance
    !*                         is for the true residual).
    !*    max_it             : Integer
    !*                         Maximum number of Krylov solver iterations.
    !*    r_pp(:)            : Real
    !*                         RHS on this process.
    !*    INTENT(INOUT)
    !*
    !*    x_pp(:)            : Real
    !*                         Solution on this process.
    !*    INTENT(IN), OPTIONAL
    !*
    !*    initial_guess_nonzero : Logical
    !*                            x_pp contains the initial guess for the
    !*                            solution.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    17.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 17.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    REAL,    INTENT(in)           :: rtol
    INTEGER, INTENT(in)           :: max_it
    REAL,    INTENT(in)           :: r_pp(:)
    REAL,    INTENT(inout)        :: x_pp(:)
    LOGICAL, INTENT(in), OPTIONAL :: initial_guess_nonzero

    ! Note that relative tolerance for PETSc is for the preconditioned
    ! residual.
    CALL p_set_tolerances(rtol,max_it)
    
    ! load vector
    CALL p_set_rhs(r_pp)

    ! Solution vector
    CALL KSPSetInitialGuessNonzero(p_object%ksp,PETSC_FALSE,p_object%ierr)
    IF (PRESENT(initial_guess_nonzero)) THEN
      IF (initial_guess_nonzero) THEN
        ! For non-linear solves, the previous solution is usually used as an
        ! initial guess: copy the ParaFEM solution vector to PETSc.
        CALL KSPSetInitialGuessNonzero(p_object%ksp,PETSC_TRUE,p_object%ierr)
        CALL p_set_solution(x_pp)
      END IF
    END IF

    CALL KSPSolve(p_object%ksp,p_object%b,p_object%x,p_object%ierr)

    CALL p_get_solution(x_pp)
  END SUBROUTINE p_solve

  SUBROUTINE p_print_info(numpe)

    !/****if* petsc/p_print_info
    !*  NAME
    !*    SUBROUTINE: p_print_info
    !*  SYNOPSIS
    !*    Usage:      p_print_info(numpe)
    !*  FUNCTION
    !*    Prints the convergence information
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    numpe              : Integer
    !*                         This processer's number.  Only processor 1
    !*                         prints information.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    29.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 29.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER, INTENT(in) :: numpe

    Vec                 :: r

    ! Preconditioned residual norm tolerance
    CALL KSPGetTolerances(p_object%ksp,p_object%rtol,PETSC_NULL_REAL,          &
                          PETSC_NULL_REAL,PETSC_NULL_INTEGER,p_object%ierr)

    ! True residual L2 norm
    CALL VecDuplicate(p_object%b,r,p_object%ierr)
    CALL KSPBuildResidual(p_object%ksp,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,r,  &
                          p_object%ierr)
    CALL VecNorm(r,NORM_2,p_object%r_n2,p_object%ierr)
    CALL VecDestroy(r,p_object%ierr)

    ! L2 norm of load
    CALL VecNorm(p_object%b,NORM_2,p_object%b_n2,p_object%ierr)
    
    CALL KSPGetIterationNumber(p_object%ksp,p_object%its,p_object%ierr)
    CALL KSPGetConvergedReason(p_object%ksp,p_object%reason,p_object%ierr)
    p_object%description = p_describe_reason(p_object%reason)
  
    IF(numpe == 1)THEN
      WRITE(11,'(A,I0,A)')                                                     &
        "The reason for convergence was ",p_object%reason," "                  &
        //TRIM(p_object%description)
      WRITE(11,'(A,I0)')                                                       &
        "The number of iterations to convergence was ",p_object%its
      WRITE(11,'(A,E17.7)')                                                    &
        "The preconditioned relative error tolerance was ",p_object%rtol
      WRITE(11,'(A,E17.7)')                                                    &
        "The relative error ||b-Ax||/||b|| was           ",                    &
        p_object%r_n2/p_object%b_n2
    END IF
  END SUBROUTINE p_print_info

  FUNCTION p_describe_reason(reason)

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
    SELECT CASE (reason)
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
  
  FUNCTION nnod(ndim,nod,over_allocation)

    !/****if* petsc/nnod
    !*  NAME
    !*    SUBROUTINE: nnod
    !*  SYNOPSIS
    !*    Usage:      nnod(ndim,nod,over_allocation)
    !*  FUNCTION

    !*    Returns an approximate upper estimate of the number of nodes in all
    !*    the neighbouring elements of a node.  This is used to estimate the
    !*    number of non-zeroes in a row of the global matrix.
    !*
    !*    If there is an error (unknown combination of dimension and number of
    !*    nodes per element), zero is returned.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    ndim               : Integer
    !*                         Number of dimensions of the problem
    !*    nod                : Integer
    !*                         Number of nodes per element
    !*    over_allocation    : Real
    !*                         over_allocation is a safety factor and allows for
    !*                         distorted grids with more than the simple number
    !*                         of elements around a point or edge.  Chose this
    !*                         to be larger than 1.25, which would correspond to
    !*                         5 hexahedra about an edge.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    19.02.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 19.02.2016, Mark Filipiak
    !*    Version 2 (was p_row_nnz), 16.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  To allow mixed elements, count the number of neighbouring nodes for all
    !*  cases, even if there isn't a node on the geometric entity (e.g. for
    !*  8-node hexahedra, calculate the number of nodes for edges).
    !*
    !*  The largest number of neighbouring nodes always occurs for a vertex
    !*  node (v).
    !*
    !*  ndim and nod are used to work out the element type, based on book
    !*  Appendix B
    !*/

    INTEGER                :: nnod
    INTEGER,  INTENT(in)   :: ndim
    INTEGER,  INTENT(in)   :: nod
    REAL,     INTENT(in)   :: over_allocation

    INTEGER                :: el_per_v,el_per_e,el_per_f,el_per_i
    REAL                   :: v=0,e=0,f=0,i=0,v_factor,e_factor

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
        v = 3; i = 2
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
        WRITE(*,'(A)') "wrong number of nodes in nnod"        
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
        v = 7; e = 4; i = 3
      CASE(6) 
        el_per_v = 6; el_per_e = 2; el_per_i = 1 ! triangle
        v = 19; e = 9; i = 6
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
        v = 9; e = 6; i = 4
      CASE(8)
        el_per_v = 4; el_per_e = 2; el_per_i = 1 ! quadrilateral
        v = 21; e = 13; i = 8
      CASE(9)
        el_per_v = 4; el_per_e = 2; el_per_i = 1 ! quadrilateral
        v = 25; e = 15; i = 9
      CASE DEFAULT
        WRITE(*,'(A)') "wrong number of nodes in nnod"        
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
        v = 13 * v_factor; e = 7 * e_factor; f = 5; i = 5
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
        v = 27; e = 18; f = 12; i = 8
      CASE(14) ! type 6 element
        el_per_v = 8; el_per_e = 4; el_per_f = 2; el_per_i = 1 ! hexahedron
        v = 63; e = 38; f = 23; i = 14
      CASE(20)
        el_per_v = 8; el_per_e = 4; el_per_f = 2; el_per_i = 1 ! hexahedron
        v = 81; e = 51; f = 32; i = 20
      CASE DEFAULT
        WRITE(*,'(A)') "wrong number of nodes in nnod"        
      END SELECT
    CASE DEFAULT
      WRITE(*,'(A)') "wrong number of dimensions in nnod"
    END SELECT
    
    nnod = NINT(over_allocation * v)
  END FUNCTION nnod

  FUNCTION row_nnz_2(ndim,nfield,nodof,ntype,nod,over_allocation)

    !/****if* petsc/row_nnz_2
    !*  NAME
    !*    SUBROUTINE: row_nnz_2
    !*  SYNOPSIS
    !*    Usage:      row_nnz_2(ndim,nfield,nodof,ntype,nod,over_allocation)
    !*  FUNCTION
    !*    Returns an approximate upper estimate of the number of non-zeroes in a
    !*    row of the global matrix.  The number of non-zeroes in a row is used
    !*    to get an approximate pre-allocation for the PETSc matrix assembly.
    !*    Too small a value will slow down the assembly, too large a value will
    !*    waste memory.
    !*
    !*    This is the most general version, for multiple fields and multiple
    !*    types of element per field.  It is assumed that there are not several
    !*    different dimension meshes (e.g. a 3D mesh coupled to a 2D mesh), thus
    !*    there is only one ndim argument.  (it is also assumed that ndim is the
    !*    dimension of the mesh, not just the dimension of the coordinates
    !*    (e.g. not a 2D mesh embedded in 3D space).
    !*
    !*    If there is an error in the argument array sizes, -1 is returned.
    !*    If some element type is not recognised, 0 is returned.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    ndim               : Integer
    !*                         Number of dimensions of the problem
    !*    nfield             : Integer
    !*                         Number of fields (e.g., velocity, pressure)
    !*    nodof              : Integer rank 1 array
    !*                         Number of dofs per node for each field
    !*    ntype              : Integer rank 1 array
    !*                         Number of different types of element for each field
    !*    nod                : Integer rank 2 array
    !*                         Number of nodes per element for each element
    !*                         type for each field.  Only
    !*                         nod(1:ntype(field),field) is used.
    !*    over_allocation    : Real
    !*                         over_allocation is a safety factor and allows for
    !*                         distorted grids with more than the simple number
    !*                         of elements around a point or edge.  Chose this
    !*                         to be larger than 1.25, which would correspond to
    !*                         5 hexahedra about an edge.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    16.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 16.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER                :: row_nnz_2
    INTEGER,  INTENT(in)   :: ndim
    INTEGER,  INTENT(in)   :: nfield
    INTEGER,  INTENT(in)   :: nodof(:)
    INTEGER,  INTENT(in)   :: ntype(:)
    INTEGER,  INTENT(in)   :: nod(:,:)
    REAL,     INTENT(in)   :: over_allocation

    INTEGER                :: field,etype
    INTEGER                :: sum,n,max_nnod

    ! Test for arrays of the wrong size or of too small a size.
    IF (SIZE(nodof) /= nfield) THEN
      WRITE(*,'(A)') "size of nodof /= nfield"
      row_nnz_2 = -1
      RETURN
    END IF
    IF (SIZE(ntype) /= nfield) THEN
      WRITE(*,'(A)') "size of ntype /= nfield"
      row_nnz_2 = -1
      RETURN
    END IF
    IF (SIZE(nod,2) /= nfield) THEN
      WRITE(*,'(A)') "size of 2nd dimension of nod /= nfield"
      row_nnz_2 = -1
      RETURN
    END IF
    IF (SIZE(nod,1) < MAXVAL(ntype)) THEN
      WRITE(*,'(A)') "size of 1st dimension of nod too small"
      row_nnz_2 = -1
      RETURN
    END IF
    
    ! Sum over all the fields (e.g., velocity and pressure) and find the
    ! maximum over all the element types for a particular field.
    sum = 0
    DO field = 1, nfield
      max_nnod = 0
      DO etype = 1, ntype(field)
        n = nnod(ndim,nod(etype,field),over_allocation)
        IF (n <= 0) THEN
          WRITE(*,'(A)') "Unknown element type"
          row_nnz_2 = 0
          RETURN
        END IF
        max_nnod = MAX(max_nnod,n)
      END DO
      sum = sum + nodof(field) * max_nnod
    END DO

    row_nnz_2 = sum    
  END FUNCTION row_nnz_2

  SUBROUTINE p_row_nnz_2(ndim,nfield,nodof,ntype,nod)

    !/****if* petsc/p_row_nnz_2
    !*  NAME
    !*    SUBROUTINE: p_row_nnz_2
    !*  SYNOPSIS
    !*    Usage:      p_row_nnz_2(ndim,nfield,nodof,ntype,nod)
    !*  FUNCTION
    !*    Sets p_object%row_nnz to an approximate upper estimate of the number
    !*    of non-zeroes in a row of the global matrix.  The number of non-zeroes
    !*    in a row is used to get an approximate pre-allocation for the PETSc
    !*    matrix assembly.  Too small a value will slow down the assembly, too
    !*    large a value will waste memory.
    !*
    !*    This is the most general version, for multiple fields and multiple
    !*    types of element per field.  It is assumed that there are not several
    !*    different dimension meshes (e.g. a 3D mesh coupled to a 2D mesh), thus
    !*    there is only one ndim argument.  (it is also assumed that ndim is the
    !*    dimension of the mesh, not just the dimension of the coordinates
    !*    (e.g. not a 2D mesh embedded in 3D space).
    !*
    !*    If there is an error in the argument array sizes, -1 is returned.
    !*    If some element type is not recognised, 0 is returned.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    ndim               : Integer
    !*                         Number of dimensions of the problem
    !*    nfield             : Integer
    !*                         Number of fields (e.g., velocity, pressure)
    !*    nodof              : Integer rank 1 array
    !*                         Number of dofs per node for each field
    !*    ntype              : Integer rank 1 array
    !*                         Number of different types of element for each field
    !*    nod                : Integer rank 2 array
    !*                         Number of nodes per element for each element
    !*                         type for each field.  Only
    !*                         nod(1:ntype(field),field) is used.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    16.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 16.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER,  INTENT(in)   :: ndim
    INTEGER,  INTENT(in)   :: nfield
    INTEGER,  INTENT(in)   :: nodof(:)
    INTEGER,  INTENT(in)   :: ntype(:)
    INTEGER,  INTENT(in)   :: nod(:,:)

    p_object%row_nnz = row_nnz_2(ndim,nfield,nodof,ntype,nod,                  &
                                 p_object%over_allocation)
  END SUBROUTINE p_row_nnz_2

  FUNCTION row_nnz_1(ndim,nodof,ntype,nod,over_allocation)

    !/****if* petsc/row_nnz_1
    !*  NAME
    !*    SUBROUTINE: row_nnz_1
    !*  SYNOPSIS
    !*    Usage:      row_nnz_1(ndim,nodof,ntype,nod,over_allocation)
    !*  FUNCTION
    !*    This is for a single field and multiple types of element for that
    !*    field.
    !*
    !*    If there is an error in the argument array sizes, -1 is returned.
    !*    If some element type is not recognised, 0 is returned.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    ndim               : Integer
    !*                         Number of dimensions of the problem
    !*    nodof              : Integer
    !*                         Number of dofs per node
    !*    ntype              : Integer
    !*                         Number of different types of element
    !*    nod                : Integer rank 1 array
    !*                         Number of nodes per element for each element
    !*                         type.
    !*    over_allocation    : Real
    !*                         over_allocation is a safety factor and allows for
    !*                         distorted grids with more than the simple number
    !*                         of elements around a point or edge.  Chose this
    !*                         to be larger than 1.25, which would correspond to
    !*                         5 hexahedra about an edge.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    16.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 16.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER                :: row_nnz_1
    INTEGER,  INTENT(in)   :: ndim
    INTEGER,  INTENT(in)   :: nodof
    INTEGER,  INTENT(in)   :: ntype
    INTEGER,  INTENT(in)   :: nod(:)
    REAL,     INTENT(in)   :: over_allocation

    IF (SIZE(nod) /= ntype) THEN
      WRITE(*,'(A)') "size of nod /= ntype"
      row_nnz_1 = -1
      RETURN
    END IF
    
    row_nnz_1 = row_nnz_2(ndim,1,(/nodof/),(/ntype/),                          &
                          RESHAPE(nod,(/SIZE(nod),1/)),over_allocation)
  END FUNCTION row_nnz_1

  SUBROUTINE p_row_nnz_1(ndim,nodof,ntype,nod)

    !/****if* petsc/p_row_nnz_1
    !*  NAME
    !*    SUBROUTINE: p_row_nnz_1
    !*  SYNOPSIS
    !*    Usage:      p_row_nnz_1(ndim,nodof,ntype,nod)
    !*  FUNCTION
    !*    Sets p_object%row_nnz to an approximate upper estimate of the number
    !*    of non-zeroes in a row of the global matrix.  The number of non-zeroes
    !*    in a row is used to get an approximate pre-allocation for the PETSc
    !*    matrix assembly.  Too small a value will slow down the assembly, too
    !*    large a value will waste memory.
    !*
    !*    This is for a single field and multiple types of element for that
    !*    field.
    !*
    !*    If there is an error in the argument array sizes, -1 is returned.
    !*    If some element type is not recognised, 0 is returned.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    ndim               : Integer
    !*                         Number of dimensions of the problem
    !*    nodof              : Integer
    !*                         Number of dofs per node
    !*    ntype              : Integer
    !*                         Number of different types of element
    !*    nod                : Integer rank 1 array
    !*                         Number of nodes per element for each element
    !*                         type.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    16.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 16.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER,  INTENT(in)   :: ndim
    INTEGER,  INTENT(in)   :: nodof
    INTEGER,  INTENT(in)   :: ntype
    INTEGER,  INTENT(in)   :: nod(:)
    
    p_object%row_nnz = row_nnz_1(ndim,nodof,ntype,nod,p_object%over_allocation)
  END SUBROUTINE p_row_nnz_1

  FUNCTION row_nnz(ndim,nodof,nod,over_allocation)

    !/****if* petsc/row_nnz
    !*  NAME
    !*    SUBROUTINE: row_nnz
    !*  SYNOPSIS
    !*    Usage:      row_nnz(ndim,nodof,nod,over_allocation)
    !*  FUNCTION
    !*    Returns an approximate upper estimate of the number of non-zeroes in a
    !*    row of the global matrix.  The number of non-zeroes in a row is used
    !*    to get an approximate pre-allocation for the PETSc matrix assembly.
    !*    Too small a value will slow down the assembly, too large a value will
    !*    waste memory.
    !*
    !*    This is for a single field and single type of element for that field.
    !*
    !*    If the element type is not recognised, 0 is returned.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    ndim               : Integer
    !*                         Number of dimensions of the problem
    !*    nodof              : Integer
    !*                         Number of dofs per node
    !*    nod                : Integer
    !*                         Number of nodes per element
    !*    over_allocation    : Real
    !*                         over_allocation is a safety factor and allows for
    !*                         distorted grids with more than the simple number
    !*                         of elements around a point or edge.  Chose this
    !*                         to be larger than 1.25, which would correspond to
    !*                         5 hexahedra about an edge.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    16.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 16.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER                :: row_nnz
    INTEGER,  INTENT(in)   :: ndim
    INTEGER,  INTENT(in)   :: nodof
    INTEGER,  INTENT(in)   :: nod
    REAL,     INTENT(in)   :: over_allocation

    row_nnz = row_nnz_1(ndim,nodof,1,(/nod/),over_allocation)
  END FUNCTION row_nnz

  SUBROUTINE p_row_nnz(ndim,nodof,nod)

    !/****if* petsc/p_row_nnz
    !*  NAME
    !*    SUBROUTINE: p_row_nnz
    !*  SYNOPSIS
    !*    Usage:      p_row_nnz(ndim,nodof,nod)
    !*  FUNCTION
    !*    Sets p_object%row_nnz to an approximate upper estimate of the number
    !*    of non-zeroes in a row of the global matrix.  The number of non-zeroes
    !*    in a row is used to get an approximate pre-allocation for the PETSc
    !*    matrix assembly.  Too small a value will slow down the assembly, too
    !*    large a value will waste memory.
    !*
    !*    This is for a single field and single type of element for that field.
    !*
    !*    If the element type is not recognised, 0 is returned.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    ndim               : Integer
    !*                         Number of dimensions of the problem
    !*    nodof              : Integer
    !*                         Number of dofs per node
    !*    nod                : Integer
    !*                         Number of nodes per element
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    16.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 16.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER,  INTENT(in)   :: ndim
    INTEGER,  INTENT(in)   :: nodof
    INTEGER,  INTENT(in)   :: nod

    p_object%row_nnz = row_nnz(ndim,nodof,nod,p_object%over_allocation)
  END SUBROUTINE p_row_nnz

END MODULE parafem_petsc
