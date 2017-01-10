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
  !*  AUTHOR
  !*    Mark Filipiak
  !*  COPYRIGHT
  !*    2016 University of Edinburgh
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  If you contribute to this module, add your author name.
  !*
  !*  PETSc will always use MPI (unless PETSc itself has been compiled without
  !*  MPI).
  !*
  !*  The arguments for PETSc routines are mostly PetscInt and PetscScalar.
  !*  Arguments with INTENT(in) are copied to local variables of the correct
  !*  type to pass to the PETSc routines.
  !*/
  
  USE, INTRINSIC :: iso_fortran_env, ONLY: int64, real32, error_unit
  USE mp_interface
  USE global_variables, ONLY: ntot, nels_pp, neq_pp, numpe
  USE gather_scatter,   ONLY: npes, neq_pp1, num_neq_pp1, threshold, ieq_start
  IMPLICIT NONE

  PUBLIC

  ! PETSc types
#include <petsc/finclude/petscdef.h>

  ! Use Fortran interfaces for the PETSc routines:  this helps infinitely in
  ! catching errors, so you get compile time errors like
  !
  ! The kind (4) of this actual argument does not match that of its associated
  ! dummy argument (8).
  !
  ! instead of run time segmentation violations, usually in a place marginally
  ! related to where the cause of the error is.
#define PETSC_USE_FORTRAN_INTERFACES
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
! PETSc Viewers can be used for testing, see p_assemble in revision 2166 for
! commented-out example.
  
  ! Private parameters
  INTEGER,          PARAMETER, PRIVATE :: string_length = 1024
  CHARACTER(len=*), PARAMETER, PRIVATE :: pre_prefix = "solver"

  ! Variables collected as one derived type
  TYPE p_type
    ! The PETSc objects cannot be initialised here because PETSC_NULL_OBJECT is
    ! a common-block-object and not a constant.
    Vec                 :: x,b
    Mat                 :: A
    PetscLogDouble, DIMENSION(MAT_INFO_SIZE) :: info
    PetscInt            :: solver   = 1
    PetscInt            :: nsolvers = 1
    PetscBool           :: nsolvers_set = .false.
    KSP,         DIMENSION(:), ALLOCATABLE :: ksp
    ! KSPGetOptionsPrefix and PCGetOptionsPrefix don't work in Fortran yet
    ! (PETSc 3.7), so keep a record of the prefixes here.
    CHARACTER(len=string_length), DIMENSION(:), ALLOCATABLE :: prefix
    PetscInt,    DIMENSION(:), ALLOCATABLE :: rows,cols
    PetscScalar, DIMENSION(:), ALLOCATABLE :: values
    PetscInt            :: its
    PetscReal           :: rtol,p_n2,r_n2,b_n2
    KSPConvergedReason  :: reason
    CHARACTER(len=string_length) :: description = ""
  END TYPE p_type

  ! Only p_object%ksp is used as a target.
  TYPE(p_type), TARGET :: p_object

  ! Private variables
  PetscMPIInt,    PRIVATE :: m_ierr
  PetscErrorCode, PRIVATE :: p_ierr

  ! Private functions
  PRIVATE :: collect_elements,row_nnz
 
CONTAINS
  
  FUNCTION p_version()

    !/****if* petsc/p_version
    !*  NAME
    !*    FUNCTION: p_version
    !*  SYNOPSIS
    !*    Usage:      p_version()
    !*  FUNCTION
    !*    Return the PETSc version.
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    14.12.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 14.12.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    USE, INTRINSIC :: iso_c_binding, ONLY: c_char,c_null_char,c_size_t

    CHARACTER(:),ALLOCATABLE :: p_version

    INTERFACE
      FUNCTION PetscGetVersion(version,length) BIND(c,name="PetscGetVersion")
        USE, INTRINSIC :: iso_c_binding, ONLY: c_char,c_size_t
        ! this assumes Fortran PetscErrorCode == C PetscErrorCode
        PetscErrorCode :: PetscGetVersion
        INTEGER(c_size_t), VALUE, INTENT(in) :: length
        CHARACTER(len=length,kind=c_char), INTENT(out) :: version
      END FUNCTION PetscGetVersion
    END INTERFACE

    INTEGER :: last
    INTEGER(c_size_t) :: length
    CHARACTER(len=string_length,kind=c_char) :: version

    length = LEN(version)
    p_ierr = PetscGetVersion(version,length)
    ! Convert null termination to spaces filling the end of the string.
    last = INDEX(version,c_null_char) - 1
    IF (last /= -1) THEN
      p_version = version(1:last)
    ELSE ! no null termination
      p_version = version
    END IF
  END FUNCTION p_version
  
  SUBROUTINE p_initialize(fname_base,error)

    !/****if* petsc/p_initialize
    !*  NAME
    !*    SUBROUTINE: p_initialize
    !*  SYNOPSIS
    !*    Usage:      p_initialize(fname_base,error)
    !*  FUNCTION
    !*      Initialises PETSc
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    fname_base         : Character
    !*                         Base name of the data file
    !*    INTENT(OUT)
    !*
    !*    error              : Logical
    !*                         error = .false. if no error occurred
    !*                         error = .true.  if an error occurred
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    19.02.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 02.06.2016, Mark Filipiak
    !*    Version 2, 16.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    CHARACTER(len=*), INTENT(in) :: fname_base
    LOGICAL, INTENT(out)         :: error

    CHARACTER(len=string_length) :: fname
    LOGICAL :: exist
    INTEGER :: ierr

    error = .false.

    fname = TRIM(fname_base)//".petsc"
    IF (numpe == 1) THEN
      INQUIRE(file=TRIM(fname),exist=exist)
    END IF
    ! The communicator is MPI_COMM_WORLD, set by find_pe_procs(), and numpe ==
    ! 1 corresponds to rank == 0.
    CALL MPI_Bcast(exist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

    IF (.NOT. exist) THEN
      IF (numpe == 1) THEN
        WRITE(error_unit,'(A)') "Error:  file "//TRIM(fname)//" not found"
      END IF
      error = .TRUE.
      RETURN
    END IF

    CALL PetscInitialize(fname,p_ierr)
  END SUBROUTINE p_initialize

  SUBROUTINE p_setup(ntot_max,g_g_pp,error)

    !/****if* petsc/p_setup
    !*  NAME
    !*    SUBROUTINE: p_setup
    !*  SYNOPSIS
    !*    Usage:      p_setup(ntot_max,g_g_pp,error)
    !*  FUNCTION
    !*      Initialises PETSc and its matrices, vectors and solvers
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
    !*    g_g_pp(:,:)        : Integer
    !*                         The steering array. g_g_pp(:,iel) are the global
    !*                         equation numbers for the dofs in element with
    !*                         local number iel (global number iel_start+iel-1)
    !*                         Restrained dofs have equation number 0.
    !*    INTENT(OUT)
    !*
    !*    error              : Logical
    !*                         error = .false. if no error occurred
    !*                         error = .true.  if an error occurred
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    31.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 07.06.2016, Mark Filipiak
    !*    Version 2, 08.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  Is ntot_max useful?  ntot is a global variable, so multiple element
    !*  types will require major changes in ParaFEM
    !*/

    INTEGER, INTENT(in)                 :: ntot_max
    INTEGER, DIMENSION(:,:), INTENT(in) :: g_g_pp
    LOGICAL, INTENT(out)                :: error

    error = .FALSE.

    ! Create the objects
    CALL p_create_matrix(g_g_pp)
    CALL p_create_vectors ! RHS and solution
    CALL p_create_ksps(error) ! Krylov solver(s)
    IF (error) THEN
      CALL p_destroy_vectors
      CALL p_destroy_matrix
      RETURN
    END IF
    CALL p_create_workspace(ntot_max)
  END SUBROUTINE p_setup

  SUBROUTINE p_finalize

    !/****if* petsc/p_finalize
    !*  NAME
    !*    SUBROUTINE: p_finalize
    !*  SYNOPSIS
    !*    Usage:      p_finalize
    !*  FUNCTION
    !*      Finalizes (i.e., tidies up) PETSc.
    !*  ARGUMENTS
    !*    None
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    18.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 31.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    CALL PetscFinalize(p_ierr)
  END SUBROUTINE p_finalize

  SUBROUTINE p_shutdown

    !/****if* petsc/p_shutdown
    !*  NAME
    !*    SUBROUTINE: p_shutdown
    !*  SYNOPSIS
    !*    Usage:      p_shutdown
    !*  FUNCTION
    !*      Destroys matrices, vectors, solvers and workspace, then finalizes
    !*      PETSc.
    !*  ARGUMENTS
    !*    None
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    31.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 31.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  There should be a p_destroy routine for each p_create routine
    !*/

    CALL p_destroy_workspace
    CALL p_destroy_ksps
    CALL p_destroy_vectors
    CALL p_destroy_matrix
    CALL p_finalize
  END SUBROUTINE p_shutdown

  SUBROUTINE p_create_matrix(g_g_pp)

    !/****if* petsc/p_create_matrix
    !*  NAME
    !*    SUBROUTINE: p_create_matrix
    !*  SYNOPSIS
    !*    Usage:      p_create_matrix(g_g_pp)
    !*  FUNCTION
    !*      Creates the global matrix (but doesn't assemble it)
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    g_g_pp(:,:)        : Integer
    !*                         The steering array. g_g_pp(:,iel) are the global
    !*                         equation numbers for the dofs in element with
    !*                         local number iel (global number iel_start+iel-1)
    !*                         Restrained dofs have equation number 0.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.04.2016
    !*  MODIFICATION HISTORY
    !*    Version 2, 8.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  PETSc 64-bit indices and 64-bit reals are used.
    !*
    !*  The slow step is calculating the number of non-zeroes in each row.
    !*
    !*  Block size fixed to 1 for just now - this cannot be set to nodof until
    !*  the restraints are handled block-wise in ParaFEM.
    !*  What about multi-field?
    !*/

    INTEGER, DIMENSION(:,:), INTENT(in) :: g_g_pp

    PetscInt :: p_neq_pp
    PetscInt, DIMENSION(:),   ALLOCATABLE :: dnz,onz
    PetscInt, DIMENSION(:,:), ALLOCATABLE :: g_g_all

    p_neq_pp = neq_pp

    ! Set the number of zeroes per row for the matrix size pre-allocation.
    ! These routines are faster and use less peak memory than using a PETSc
    ! MATPREALLOCATOR matrix.
    CALL collect_elements(g_g_pp,g_g_all) ! g_g_all allocated
    CALL row_nnz(g_g_all,dnz,onz) ! dnz, onz allocated
    DEALLOCATE(g_g_all)

    CALL MatCreate(PETSC_COMM_WORLD,p_object%A,p_ierr)
    CALL MatSetSizes(p_object%A,p_neq_pp,p_neq_pp,                             &
                     PETSC_DETERMINE,PETSC_DETERMINE,p_ierr)
    ! The default matrix type set by MatSetFromOptions is MATAIJ, which is what
    ! we want.  MatSetFromOptions is called before specific options are set
    ! here, so that the user cannot override the specific options set here.
    CALL MatSetFromOptions(p_object%A,p_ierr)
    
    CALL MatSeqAIJSetPreallocation(p_object%A,                                 &
                                   PETSC_NULL_INTEGER,dnz,                     &
                                   p_ierr)
    CALL MatMPIAIJSetPreallocation(p_object%A,                                 &
                                   PETSC_NULL_INTEGER,dnz,                     &
                                   PETSC_NULL_INTEGER,onz,                     &
                                   p_ierr)
    DEALLOCATE(dnz,onz)

    ! PETSc uses C array order, so a transpose would be needed when adding
    ! elements, but you can get PETSc to use Fortran order for MatSetValues by
    ! setting MAT_ROW_ORIENTED false (at least for the Mat types used so far:
    ! SEQAIJ and MPIAIJ).
    CALL MatSetOption(p_object%A,MAT_ROW_ORIENTED,PETSC_FALSE,p_ierr)
  END SUBROUTINE p_create_matrix

  SUBROUTINE p_destroy_matrix

    !/****if* petsc/p_destroy_matrix
    !*  NAME
    !*    SUBROUTINE: p_destroy_matrix
    !*  SYNOPSIS
    !*    Usage:      p_destroy_matrix
    !*  FUNCTION
    !*      Destroys the global matrix
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    31.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 31.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    CALL MatDestroy(p_object%A,p_ierr)
  END SUBROUTINE p_destroy_matrix

  SUBROUTINE p_zero_matrix

    !/****if* petsc/p_zero_matrix
    !*  NAME
    !*    SUBROUTINE: p_zero_matrix
    !*  SYNOPSIS
    !*    Usage:      p_zero_matrix
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

    CALL MatZeroEntries(p_object%A,p_ierr)
  END SUBROUTINE p_zero_matrix

  SUBROUTINE collect_elements(g_g_pp,g_g_all)

    !/****if* petsc/collect_elements
    !*  NAME
    !*    SUBROUTINE: collect_elements
    !*  SYNOPSIS
    !*    Usage:      collect_elements(g_g_pp,g_g_all)
    !*  FUNCTION
    !*    Collects onto this process all elements that contain dofs
    !*    corresponding to the equations on this process.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    g_g_pp(:,:)        : Integer
    !*                         The steering array. g_g_pp(:,iel) are the global
    !*                         equation numbers for the dofs in element with
    !*                         local number iel (global number iel_start+iel-1)
    !*                         Restrained dofs have equation number 0.
    !*    INTENT(OUT)
    !*
    !*    g_g_all(:,:)       : PetscInt
    !*                         The equivalent of g_g_pp, but for all elements
    !*                         that have dofs corresponding to the equations on
    !*                         this process.  The element numbers of g_g_all
    !*                         have no relationship with global or local
    !*                         element numbers:  the element numbers are not
    !*                         needed for the row non-zero counting that
    !*                         g_g_all is used for.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    08.08.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 08.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  PETSc 64-bit indices are used.
    !*/

    INTEGER,  DIMENSION(:,:),              INTENT(in)  :: g_g_pp
    PetscInt, DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: g_g_all

    PetscInt    :: i,n,iel,j,low,p,s,nels_send,nels_recv,nels_all,r
    PetscMPIInt :: n_sends,n_recvs,ierr
    PetscInt,    DIMENSION(:),   ALLOCATABLE :: q
    PetscInt,    DIMENSION(:,:), ALLOCATABLE :: el_p,p_el,send_buffer
    ! In theory you could have more than 2**31 - 1 elements on a process, but
    ! the counts in MPI messages can only be MPI_INTEGER size (usually 32-bit).
    PetscInt,    DIMENSION(:),   ALLOCATABLE :: send_start,recv_start
    PetscMPIInt, DIMENSION(:),   ALLOCATABLE :: send_count,recv_count

    PetscMPIInt, DIMENSION(:), ALLOCATABLE :: requests
    PetscMPIInt                            :: PARAFEM_ELEMENT

    ! el_p defines the element to process mapping
    ! el_p(:,iel) are the process numbers of equations corresponding to dofs
    ! contained in element iel
    ALLOCATE(el_p(ntot,nels_pp))

    ! See the gather_scatter module, especially calc_neq_pp and make_ggl.
    WHERE (g_g_pp == 0) ! restraints
      el_p = 0
    ELSEWHERE (g_g_pp <= threshold .OR. threshold == 0)
      el_p = (g_g_pp - 1)/neq_pp1 + 1
    ELSEWHERE
      el_p = num_neq_pp1 + (g_g_pp - threshold - 1)/(neq_pp1 - 1) + 1
    END WHERE

    ! p_el defines one stage of the mapping (process-element pairs)
    ! from: process number of equations corresponding to dofs contained in
    !       elements on this process, but only for equations that are not
    !       on this process.  For equations on this process corresponding
    !       to dofs contained in elements on this process, g_g_pp contains
    !       the information needed, and g_g_pp becomes part of g_g_all
    !       later.
    ! to:   the corresponding local number of the element.
    ! p_el(:,1) are process numbers
    ! p_el(:,2) are local element numbers
    ! Done this way to use Petsc's sort routines.
    ALLOCATE(p_el(ntot*nels_pp,2))
    ALLOCATE(q(ntot)) ! work array

    i = 1
    DO iel = 1, nels_pp
      q = el_p(:,iel)
      n = ntot ! n becomes the number of unique entries
      CALL PetscSortRemoveDupsInt(n,q,p_ierr)
      ! remove a possible zero entry (from restraints)
      IF (q(1) /= 0) THEN
        low = 1
      ELSE
        low = 2
      END IF
      DO j = low, n
        ! Equations that are on this process, corresponding to dofs contained in
        ! elements on this process don't need to be considered: the elements are
        ! already on this process in the g_g_pp array.
        IF (q(j) /= numpe) THEN
          p_el(i,1) = q(j)
          p_el(i,2) = iel
          i = i + 1
        END IF
      END DO
    END DO
    n = i - 1
    DEALLOCATE(el_p,q)
    ! => p_el(1:n,:) gives the pairing from process to element, but is unsorted

    CALL PetscSortIntWithArray(n,p_el(1:n,1),p_el(1:n,2),p_ierr)
    ! => p_el(1:n,:) gives the pairing from process to element, sorted by
    !    process number.

    ! Now complete the mapping.
    ALLOCATE(send_start(npes),send_count(npes))
    send_count = 0
    p = 0 ! an invalid process number to start things off
    DO i = 1, n
      IF (p_el(i,1) /= p) THEN
        ! next process in mapping
        p = p_el(i,1)
      END IF
      send_count(p) = send_count(p) + 1
    END DO
    s = 1 ! s for start
    DO p = 1, npes
      send_start(p) = s
      s = s + send_count(p)
    END DO
    n_sends = COUNT(send_count /= 0)
    ! => if send_count(p) /= 0 then the elements to send to process p are
    !    p_el(send_start(p):send_start(p)+send_count(p)-1,2)
    ! => process to element mapping complete and send information complete.

    ! For each process p, you could sort the element numbers in
    ! p_el(send_start(p):send_start(p)+send_count(p)-1,2) on their own (no need
    ! to include p_el(send_start(p):send_start(p)+send_count(p)-1,1), which are
    ! all p).  Having the elements ordered for each process may improve the
    ! locality of send_buffer = g_g_pp(:,p_el(1:n,2)).  But you would have to
    ! measure if the time for the sorting is compensated by the time
    ! improvement in the copy.

    ! Transpose the send_count to get the receive counts
    ALLOCATE(recv_start(npes),recv_count(npes))
    CALL MPI_Alltoall(send_count,1,MPI_INTEGER,recv_count,1,MPI_INTEGER,       &
                      MPI_COMM_WORLD,ierr)
    
    s = 1 ! s for start
    DO p = 1, npes
      recv_start(p) = s
      s = s + recv_count(p)
    END DO
    n_recvs = COUNT(recv_count /= 0)
    ! => process to element mapping complete and send and receive information
    !    complete

    ! Collect all the elements required by equations that are on this process.
    ! The elements from other processes will be put into the beginning of
    ! g_g_all, then the elements on this process will be copied from g_g_pp to
    ! the end of g_g_all.  g_g_all can be the receive buffer because it is
    ! contiguous, but there needs to be a contiguous send buffer to use
    ! MPI_Isend.  This would be a problem if you are sending lots of data but it
    ! is expected that the ratio of neigbouring elements (or 'surface') to local
    ! elements (the 'volume') is small - otherwise you have used too many
    ! processes for the size of problem and/or chosen a poor partitioning of the
    ! mesh.  To reduce peak memory, allocate g_g_all after send_buffer is set up
    ! and p_el is freed.
    nels_send = SUM(send_count) ! = n
    nels_recv = SUM(recv_count)
    nels_all = nels_recv + nels_pp
    ALLOCATE(send_buffer(ntot,nels_send))
    ! TODO: check that zero sizes for arrays and zero counts are OK.
    ALLOCATE(requests(n_sends + n_recvs))

    ! Copy the elements to the send buffer, in the correct order.
    ! You can do
    send_buffer = g_g_pp(:,p_el(1:n,2))
    ! instead of
    ! DO p = 1,npes
    !   send_buffer(:,send_start(p):send_start(p)+send_count(p)-1)               &
    !     = g_g_pp(:,p_el(send_start(p):send_start(p)+send_count(p)-1,2))
    ! END DO
    ! Isn't Fortran wonderful!
    DEALLOCATE(p_el)

    ALLOCATE(g_g_all(ntot,nels_all))

    ! An MPI datatype to hold an element.  MPIU_INTEGER is the MPI type for a
    ! PetscInt, correctly set to be 32- or 64-bit depending on how PETSc was
    ! built.
    CALL MPI_Type_contiguous(ntot,MPIU_INTEGER,PARAFEM_ELEMENT,ierr)
    CALL MPI_Type_commit(PARAFEM_ELEMENT,ierr)

    ! If you swap the order of the receive and send loops, then make sure r is
    ! handled correctly.
    r = 1 ! r for requests
    DO p = 1, npes
      IF (recv_count(p) /= 0) THEN
        CALL MPI_Irecv(g_g_all(:,recv_start(p):recv_start(p)+recv_count(p)-1), &
                       recv_count(p),PARAFEM_ELEMENT,p-1,0,                    &
                       MPI_COMM_WORLD,requests(r),ierr)
        r = r + 1
      END IF
    END DO
    IF (r - 1 /= n_recvs) THEN ! disaster
      WRITE(error_unit,'(A,I6)') "Error: r - 1 /= n_recvs on process ", numpe
      CALL MPI_Abort(MPI_COMM_WORLD,ierr)
    END IF

    DO p = 1, npes
      IF (send_count(p) /= 0) THEN
        CALL MPI_Isend(send_buffer                                             &
                         (:,send_start(p):send_start(p)+send_count(p)-1),      &
                       send_count(p),PARAFEM_ELEMENT,p-1,0,                    &
                       MPI_COMM_WORLD,requests(r),ierr)
        r = r + 1
      END IF
    END DO
    IF (r - 1 /= n_recvs + n_sends) THEN ! disaster
      WRITE(error_unit,'(A,I6)')                                               &
        "Error: r - 1 /= n_recvs + n_sends on process ", numpe
      CALL MPI_Abort(MPI_COMM_WORLD,ierr)
    END IF

    CALL MPI_Waitall(n_recvs + n_sends,requests,MPI_STATUSES_IGNORE,ierr)
    IF (ierr /= MPI_SUCCESS) THEN ! disaster
      WRITE(error_unit,'(A,I6)') "Error: receive/send error on process ", numpe
      CALL MPI_Abort(MPI_COMM_WORLD,ierr)
    END IF

    CALL MPI_Type_free(PARAFEM_ELEMENT,ierr)
    DEALLOCATE(requests,send_buffer,recv_start,recv_count,send_start,send_count)
    ! => all the off-process elements are in g_g_all

    g_g_all(:,nels_recv+1:) = g_g_pp
    ! => all the elements needed for the equations that are on this process are
    !    in g_g_all.  The element numbers of g_g_all have no relationship with
    !    global or local element numbers.
  END SUBROUTINE collect_elements

  SUBROUTINE row_nnz(g_g_all,dnz,onz)

    !/****if* petsc/row_nnz
    !*  NAME
    !*    SUBROUTINE: row_nnz
    !*  SYNOPSIS
    !*    Usage:      row_nnz(g_g_all,dnz,onz)
    !*  FUNCTION
    !*    Calculates the number of on- and off-process non-zeroes in the global
    !*    matrix for the equations on this process.
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    g_g_all(:,:)       : PetscInt
    !*                         The equivalent of g_g_pp, but for all elements
    !*                         that have dofs corresponding to the equations on
    !*                         this process.  The element numbers in g_g_all
    !*                         have no relationship with global or local
    !*                         element numbers:  the element numbers are not
    !*                         needed for the row non-zero counting that
    !*                         g_g_all is used for.
    !*    INTENT(OUT)
    !*
    !*    dnz(:)             : PetscInt
    !*                         dnz(r) is the number of non-zeroes corresponding
    !*                         to equations on this process in the global
    !*                         matrix for the equation with local number r.
    !*    onz(:)             : PetscInt
    !*                         onz(r) is the number of non-zeroes corresponding
    !*                         to equations not on this process in the global
    !*                         matrix for the equation with local number r.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    08.08.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 08.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  PETSc 64-bit indices are used.
    !*
    !*  The slow steps are the sort of eq_el and the final counting of
    !*  non-zeroes for each equation.  These two steps take about the same
    !*  time.  In the counting loop, the sorting takes the most time.
    !*/

    PetscInt, DIMENSION(:,:),              INTENT(in)  :: g_g_all
    PetscInt, DIMENSION(:),   ALLOCATABLE, INTENT(out) :: dnz,onz

    PetscInt :: nels_all,ieq_finish,i,n,iel,j,low,ieq,s,nnz_pp,nnz
    PetscMPIInt :: ierr
    PetscInt, DIMENSION(:),   ALLOCATABLE :: d,eq_start,eq_count
    PetscInt, DIMENSION(:,:), ALLOCATABLE :: eq_el

    nels_all = SIZE(g_g_all,2)

    ! eq_el defines one stage of the mapping (equation-element pairs)
    ! from: local equation number of equations on this process
    ! to:   the indices in g_g_all of the elements that contain the equation
    ! eq_el(:,1) are local equation numbers (eventually)
    ! eq_el(:,2) are element numbers in g_g_all
    ! Done this way to use Petsc's sort routines.
    ALLOCATE(eq_el(ntot*nels_all,2))
    ALLOCATE(d(ntot)) ! work array

    ieq_finish = ieq_start + neq_pp -1
    i = 1
    DO iel = 1, nels_all
      d = g_g_all(:,iel) ! global equation numbers
#ifdef P_DUPLICATE_DOFS
      ! Is there any possibility at all of duplicate dofs in an element?
      n = ntot ! n becomes the number of unique entries
      CALL PetscSortRemoveDupsInt(n,d,p_ierr)
      ! There is no need to shift the start of the loop if d(1) == 0, since
      ! ieq_start > 0
      DO j = 1, n
        IF (ieq_start <= d(j) .AND. d(j) <= ieq_finish) then
          eq_el(i,1) = d(j)
          eq_el(i,2) = iel
          i = i + 1
        END IF
      END DO
#else
      DO j = 1, ntot
        ! There is no need to test for d(j) /= 0, since ieq_start > 0
        IF (ieq_start <= d(j) .AND. d(j) <= ieq_finish) THEN
          eq_el(i,1) = d(j)
          eq_el(i,2) = iel
          i = i + 1
        END IF
      END DO
#endif
    END DO
    n = i - 1
    DEALLOCATE(d)
    ! => eq_el(1:n,:) gives the pairing from global equation number (of
    !    equations on this process) to element, but is unsorted

    ! Sort eq_el.  This is slow.
    CALL PetscSortIntWithArray(n,eq_el(1:n,1),eq_el(1:n,2),p_ierr)
    ! => eq_el(1:n,:) gives the pairing from global equation number (of
    !    equations on this process) to element, sorted by global equation
    !    number.

    eq_el(1:n,1) = eq_el(1:n,1) - ieq_start + 1
    ! => eq_el(1:n,:) gives the pairing from local equation number (of equations
    !    on this process) to element, sorted by local equation number.

    ! Now complete the mapping.
    ALLOCATE(eq_start(neq_pp),eq_count(neq_pp))
    eq_count = 0
    ieq = 0 ! an invalid equation number to start things off.
    DO i = 1, n
      IF (eq_el(i,1) /= ieq) THEN
        ! next equation in mapping
        ieq = eq_el(i,1)
      END IF
      eq_count(ieq) = eq_count(ieq) + 1
    END DO

    ! Will it be at all possible to have an equation with no elements?
    IF (MINVAL(eq_count) < 1) THEN
      WRITE(error_unit,'(A,I0,A)')                                             &
        "Warning: equation ", ieq_start + MINLOC(eq_count),                    &
        " is not in any element"
    END IF
      
    s = 1 ! s for start
    DO ieq = 1, neq_pp
      eq_start(ieq) = s
      s = s + eq_count(ieq)
    END DO
    ! => if eq_count(ieq) /= 0 then the elements with the dof that corresponds
    !    to ieq are eq_el(eq_start(ieq):eq_start(ieq)+eq_count(ieq)-1,2)
    ! => equation to element mapping complete.

    ! For each equation on this process, count the on- and off-process
    ! non-zeroes in the corresponding row of the global matrix.  This is slow.

    ALLOCATE(dnz(neq_pp),onz(neq_pp))
    ALLOCATE(d(ntot*MAXVAL(eq_count))) ! work array

    dnz = 0
    onz = 0
    nnz_pp = 0 ! number of non-zeroes in the global matrix on this process
    DO ieq = 1, neq_pp
      IF (eq_count(ieq) /= 0) THEN
        n = ntot * eq_count(ieq)
        d(1:n) =                                                               &
          RESHAPE(g_g_all                                                      &
                    (:,eq_el(eq_start(ieq):eq_start(ieq)+eq_count(ieq)-1,2)),  &
                  (/n/))
        CALL PetscSortRemoveDupsInt(n,d,p_ierr)
        ! remove a possible zero entry (from restraints).
        IF (d(1) /= 0) THEN
          low = 1
        ELSE
          low = 2
        END IF
        dnz(ieq) = COUNT(ieq_start <= d(low:n) .AND. d(low:n) <= ieq_finish)
        onz(ieq) = n - low + 1 - dnz(ieq) 
        ! == onz(ieq) = COUNT(d(low:n) < ieq_start .OR. ieq_finish < d(low:n))
        nnz_pp = nnz_pp + n - low + 1
      END IF
    END DO
    DEALLOCATE(d,eq_start,eq_count,eq_el)
    ! number of non-zeroes in the global matrix
    CALL MPI_Reduce(nnz_pp,nnz,1,MPIU_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    IF (numpe == 1) THEN
      WRITE(*,'(A,I0)')"Number of non-zeroes in the global matrix = ", nnz
    END IF
  END SUBROUTINE row_nnz

  SUBROUTINE p_create_vectors

    !/****if* petsc/p_create_vectors
    !*  NAME
    !*    SUBROUTINE: p_create_vectors
    !*  SYNOPSIS
    !*    Usage:      p_create_vectors
    !*  FUNCTION
    !*    Create the global RHS and solution vectors
    !*  ARGUMENTS
    !*    None
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

    PetscInt :: p_neq_pp

    p_neq_pp = neq_pp

    ! RHS vector.  For this particular order (create, set size, set from
    ! options) the allocation is done by set from options.
    CALL VecCreate(PETSC_COMM_WORLD,p_object%b,p_ierr)
    CALL VecSetSizes(p_object%b,p_neq_pp,PETSC_DECIDE,p_ierr)
    ! The default vector type set by VecSetFromOptions is VECSEQ for one
    ! process and VECMPI for more than one process.
    CALL VecSetFromOptions(p_object%b,p_ierr)

    ! Solution vector
    CALL VecDuplicate(p_object%b,p_object%x,p_ierr)
  END SUBROUTINE p_create_vectors
  
  SUBROUTINE p_destroy_vectors

    !/****if* petsc/p_destroy_vectors
    !*  NAME
    !*    SUBROUTINE: p_destroy_vectors
    !*  SYNOPSIS
    !*    Usage:      p_destroy_vectors
    !*  FUNCTION
    !*    Destroy the global RHS and solution vectors
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    31.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 31.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    CALL VecDestroy(p_object%x,p_ierr)
    CALL VecDestroy(p_object%b,p_ierr)
  END SUBROUTINE p_destroy_vectors
  
  SUBROUTINE p_create_ksps(error)

    !/****if* petsc/p_create_ksps
    !*  NAME
    !*    SUBROUTINE: p_create_ksps
    !*  SYNOPSIS
    !*    Usage:      p_create_ksps(error)
    !*  FUNCTION
    !*    Create the Krylov solvers and preconditioners
    !*  ARGUMENTS
    !*    INTENT(OUT)
    !*
    !*    error              : Logical
    !*                         error = .false. if no error occurred
    !*                         error = .true.  if an error occurred
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    28.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 13.06.2016, Mark Filipiak
    !*    Version 2, 15.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  KSP type, per-KSP tolerances (rtol, abstol, dtol, maxits), KSP
    !*  options, PC type, PC options are set in the xxx.petsc file.  Those
    !*  options are used to set up the preconditioned Krylov solver.  If there
    !*  are several KSP types to be chosen from, then each one will be
    !*  bracketed by -prefix_push and -prefix_pop.  For example
    !* 
    !*  -nsolvers 2
    !* 
    !*  -prefix_push solver_1_
    !*    -ksp_type cg
    !*  -prefix_pop
    !* 
    !*  -prefix_push solver_2_
    !*    -ksp_type gmres
    !*  -prefix_pop
    !* 
    !*  in the xxx.petsc file and 
    !* 
    !*  CALL KSPSetOptionsPrefix(p_object%ksp,"solver_1_",p_ierr)
    !* 
    !*  before KSPSetFromOptions before using the 'solver_1_' solver, i.e. CG,
    !*  and
    !* 
    !*  CALL KSPSetOptionsPrefix(p_object%ksp,"solver_2_",p_ierr)
    !* 
    !*  before KSPSetFromOptions before using the 'solver_2_' solver, i.e.
    !*  GMRES.  Thus you can switch from CG to GMRES (for example) during a
    !*  simulation.
    !* 
    !*  For the simple case of one solver, -nsolvers and the
    !*  -prefix_push/-prefix_pop pair can be omitted.  For example
    !* 
    !*  -ksp_type minres
    !* 
    !*  in the xxx.petsc file and no KSPSetOptionsPrefix used.
    !*/

    LOGICAL, INTENT(out) :: error

    PetscInt                     :: s
    CHARACTER(len=string_length) :: s_string
    PetscBool                    :: set
    PC                           :: p
    KSPType                      :: ksp_type
    PCType                       :: pc_type

    error = .FALSE.

    ! p_object%nsolvers is initialized to 1.  PetscOptionsGetInt does not change
    ! p_object%nsolvers if -nsolvers is not set in xxx.petsc
    CALL PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,            &
                            "-nsolvers",p_object%nsolvers,                     &
                            p_object%nsolvers_set,p_ierr)

    IF (p_object%nsolvers < 1) THEN
      IF (numpe == 1) THEN
        WRITE(error_unit,'(A)') "Error:  -nsolvers is less than 1"
      END IF
      error = .TRUE.
      RETURN
    END IF

    ALLOCATE(p_object%ksp(p_object%nsolvers),                                  &
             p_object%prefix(p_object%nsolvers))

    IF (.NOT. p_object%nsolvers_set) THEN
      ! Simple case, one solver, no prefix needed in the xxx.petsc file
      s = 1
      CALL KSPCreate(PETSC_COMM_WORLD,p_object%ksp(s),p_ierr)
      p_object%prefix(s) = "" ! for tests
      ! No KSPSetOptionsPrefix
      CALL KSPSetFromOptions(p_object%ksp(s),p_ierr)
    ELSE
      ! General case, one or more solvers, prefixes needed
      DO s = 1, p_object%nsolvers
        CALL KSPCreate(PETSC_COMM_WORLD,p_object%ksp(s),p_ierr)
        WRITE(s_string,'(I0)') s
        p_object%prefix(s) = pre_prefix//"_"//TRIM(s_string)//"_"
        CALL KSPSetOptionsPrefix(p_object%ksp(s),p_object%prefix(s),p_ierr)
        CALL KSPSetFromOptions(p_object%ksp(s),p_ierr)
      END DO
    END IF

    ! Fail if the KSP, tolerance, max iterations, or PC have not been set.
    ! PetscOptionsHasName(NULL,"","-ksp_type",...) will look for the option
    ! "-ksp_type", i.e. as if there were no prefix.
    DO s = 1, p_object%nsolvers
      CALL PetscOptionsHasName(PETSC_NULL_OBJECT,p_object%prefix(s),           &
                               "-ksp_type",set,p_ierr)
      IF (.NOT. set) THEN
        IF (numpe == 1) THEN
          WRITE(error_unit,'(A)')                                              &
            "Error:  -"//TRIM(p_object%prefix(s))//"ksp_type not set"
        END IF
        error = .TRUE.
      END IF
      CALL PetscOptionsHasName(PETSC_NULL_OBJECT,p_object%prefix(s),           &
                               "-ksp_rtol",set,p_ierr)
      IF (.NOT. set) THEN
        IF (numpe == 1) THEN
          WRITE(error_unit,'(A)')                                              &
            "Error:  -"//TRIM(p_object%prefix(s))//"ksp_rtol not set"
        END IF
        error = .TRUE.
      END IF
      CALL PetscOptionsHasName(PETSC_NULL_OBJECT,p_object%prefix(s),           &
                               "-ksp_max_it",set,p_ierr)
      IF (.NOT. set) THEN
        IF (numpe == 1) THEN
          WRITE(error_unit,'(A)')                                              &
            "Error:  -"//TRIM(p_object%prefix(s))//"ksp_max_it not set"
        END IF
        error = .TRUE.
      END IF
      CALL PetscOptionsHasName(PETSC_NULL_OBJECT,p_object%prefix(s),           &
                               "-pc_type",set,p_ierr)
      IF (.NOT. set) THEN
        IF (numpe == 1) THEN
          WRITE(error_unit,'(A)')                                              &
            "Error:  -"//TRIM(p_object%prefix(s))//"pc_type not set"
        END IF
        error = .TRUE.
      END IF
    END DO
    IF (error) THEN
      CALL p_destroy_ksps
      RETURN
    END IF
  END SUBROUTINE p_create_ksps
  
  SUBROUTINE p_destroy_ksps

    !/****if* petsc/p_destroy_ksps
    !*  NAME
    !*    SUBROUTINE: p_destroy_ksps
    !*  SYNOPSIS
    !*    Usage:      p_destroy_ksps
    !*  FUNCTION
    !*    Destroy the Krylov solver and preconditioner
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    31.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 02.06.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    PetscInt :: s

    DO s = 1, p_object%nsolvers
      CALL KSPDestroy(p_object%ksp(s),p_ierr)
    END DO
    DEALLOCATE(p_object%ksp,p_object%prefix)
  END SUBROUTINE p_destroy_ksps
  
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
    !*  PETSc 64-bit indices and 64-bit reals are used.
    !*/

    INTEGER, INTENT(in) :: ntot_max

    ! Arrays of PETSc types to be proof against changes in index and scalar
    ! sizes.
    ALLOCATE(p_object%rows(ntot_max),p_object%cols(ntot_max),                  &
             p_object%values(ntot_max*ntot_max))
  END SUBROUTINE p_create_workspace

  SUBROUTINE p_destroy_workspace

    !/****if* petsc/p_destroy_workspace
    !*  NAME
    !*    SUBROUTINE: p_destroy_workspace
    !*  SYNOPSIS
    !*    Usage:      p_destroy_workspace
    !*  FUNCTION
    !*      Destroys some workspace
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    31.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 31.05.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    DEALLOCATE(p_object%rows,p_object%cols,p_object%values)
  END SUBROUTINE p_destroy_workspace

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
    !*
    !*  MAT_IGNORE_ZERO_ENTRIES is not set:  for a series of non-linear solves
    !*  with the same adjacency matrix for the mesh, the matrix non-zero
    !*  structure will remain fixed.
    !*/

    INTEGER, INTENT(in) :: g(:)
    REAL,    INTENT(in) :: km(:,:)

    PetscInt :: p_ntot

    p_ntot = SIZE(g)
    p_object%rows(1:p_ntot) = g - 1
    p_object%cols(1:p_ntot) = p_object%rows(1:p_ntot)
    ! PETSc uses C array order, so normally a transpose would be needed, but
    ! MAT_ROW_ORIENTED is set to false in p_create_matrix to make MatSetValues
    ! work with Fortran-order arrays.
    p_object%values(1:p_ntot*p_ntot) = RESHAPE(km,(/p_ntot*p_ntot/))
    CALL MatSetValues(p_object%A,p_ntot,p_object%rows,p_ntot,p_object%cols,    &
                      p_object%values,ADD_VALUES,p_ierr)
  END SUBROUTINE p_add_element

  SUBROUTINE p_assemble

    !/****if* petsc/p_assemble
    !*  NAME
    !*    SUBROUTINE: p_assemble
    !*  SYNOPSIS
    !*    Usage:      p_assemble
    !*  FUNCTION
    !*    Assemble the global matrix.
    !*  ARGUMENTS
    !*    None
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

    CALL MatAssemblyBegin(p_object%A,MAT_FINAL_ASSEMBLY,p_ierr)
    CALL MatAssemblyEnd(p_object%A,MAT_FINAL_ASSEMBLY,p_ierr)

    ! All subsequent assemblies SHOULD not create any new entries: fail with an
    ! error message if they do.
    CALL MatSetOption(p_object%A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE,p_ierr)
  END SUBROUTINE p_assemble

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

    REAL, INTENT(in) :: r_pp(:)

    PetscScalar,POINTER :: varray(:)

    CALL VecGetArrayF90(p_object%b,varray,p_ierr)
    ! This is OK as long as PetscScalars are not smaller than ParaFEM reals.
    ! There should be a test for sizes of PetscScalars (and PetscInts) and
    ! ParaFEM reals (and indices).
    varray = r_pp
    CALL VecRestoreArrayF90(p_object%b,varray,p_ierr)
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

    REAL, INTENT(in) :: x_pp(:)

    PetscScalar,POINTER :: varray(:)

    CALL VecGetArrayF90(p_object%x,varray,p_ierr)
    ! This is OK as long as PetscScalars are not smaller than ParaFEM reals.
    ! There should be a test for sizes of PetscScalars (and PetscInts) and
    ! ParaFEM reals (and indices).
    varray = x_pp
    CALL VecRestoreArrayF90(p_object%x,varray,p_ierr)
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

    REAL, INTENT(out) :: x_pp(:)

    PetscScalar,POINTER :: varray(:)

    CALL VecGetArrayF90(p_object%x,varray,p_ierr)
    ! This is OK as long as ParaFEM reals are not smaller than PetscScalars.
    ! There should be a test for sizes of PetscScalars (and PetscInts) and
    ! ParaFEM reals (and indices).
    x_pp = varray
    CALL VecRestoreArrayF90(p_object%x,varray,p_ierr)
  END SUBROUTINE p_get_solution

  SUBROUTINE p_use_solver(solver,error)

    !/****if* petsc/p_use_solver
    !*  NAME
    !*    SUBROUTINE: p_use_solver
    !*  SYNOPSIS
    !*    Usage:      p_use_solver(solver,error)
    !*  FUNCTION
    !*    Choose one of the PETSc solvers specified by options and read in by
    !*    p_create_ksps (which is called by p_setup).
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    solver                : Integer
    !*                            The number of the solver to use.  This allows
    !*                            you to use different solvers at different
    !*                            stages of a calculation.
    !*    error                 : Logical
    !*                            error = .false. if no error occurred
    !*                            error = .true.  if an error occurred
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    02.06.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 13.06.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER, INTENT(in)  :: solver
    LOGICAL, INTENT(out) :: error

    error = .FALSE.

    IF (solver < 1 .OR. solver > p_object%nsolvers) THEN
      IF (numpe == 1) THEN
        WRITE(error_unit,'(A,I0,A,I0)')                                        &
          "Error: p_use_solver uses solver ",solver,                           &
          " but nsolvers is ",p_object%nsolvers
      END IF
      error = .TRUE.
      RETURN
    END IF

    p_object%solver = solver
  END SUBROUTINE p_use_solver

  SUBROUTINE p_solve(r_pp,x_pp,initial_guess_nonzero,reuse_preconditioner)

    !/****if* petsc/p_solve
    !*  NAME
    !*    SUBROUTINE: p_solve
    !*  SYNOPSIS
    !*    Usage:      p_solve(r_pp,x_pp,
    !*                        initial_guess_nonzero,reuse_preconditioner)
    !*  FUNCTION
    !*    Solve using PETSc, using the only solver, or the current solver
    !*    chosen using p_use_solver.  Tolerance and maximum number of
    !*    iterations are set in the PETSc control file (xxxx.petsc).
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
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
    !*    reuse_preconditioner  : Logical
    !*                            Normally the preconditioner is recalculated
    !*                            before each solve because in non-linear
    !*                            problems the matrix will change for each
    !*                            solve.  You may want to save time (usually
    !*                            at the cost of more linear solver iterations)
    !*                            by recalculating the preconditioner less
    !*                            frequently than every non-linear iteration.
    !*                            Set this flag to .true. to reuse the existing
    !*                            preconditioner (obviously the very first
    !*                            solve cannot have this flag set to .true.).
    !*
    !*                            You don't need to use this argument if you
    !*                            are doing several linear solves with the same
    !*                            matrix:  if the matrix hasn't changed, the
    !*                            preconditioner is not recalculated.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    17.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 13.06.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    REAL,    INTENT(in)           :: r_pp(:)
    REAL,    INTENT(inout)        :: x_pp(:)
    LOGICAL, INTENT(in), OPTIONAL :: initial_guess_nonzero
    LOGICAL, INTENT(in), OPTIONAL :: reuse_preconditioner

    KSP,POINTER :: ksp

    ! RHS vector
    CALL p_set_rhs(r_pp)

    ! Reduce typing
    ksp => p_object%ksp(p_object%solver)
    
    ! Solution vector 
    CALL KSPSetInitialGuessNonzero(ksp,PETSC_FALSE,p_ierr)
    IF (PRESENT(initial_guess_nonzero)) THEN
      IF (initial_guess_nonzero) THEN
        ! For non-linear solves, the previous solution is usually used as an
        ! initial guess: copy the ParaFEM solution vector to PETSc.
        CALL KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,p_ierr)
        CALL p_set_solution(x_pp)
      END IF
    END IF

    CALL KSPSetOperators(ksp,p_object%A,p_object%A,p_ierr)

    CALL KSPSetReusePreconditioner(ksp,PETSC_FALSE,p_ierr)
    IF (PRESENT(reuse_preconditioner)) THEN
      IF (reuse_preconditioner) THEN
        CALL KSPSetReusePreconditioner(ksp,PETSC_TRUE,p_ierr)
      END IF
    END IF

    CALL KSPSolve(ksp,p_object%b,p_object%x,p_ierr)

    CALL p_get_solution(x_pp)
  END SUBROUTINE p_solve

  SUBROUTINE p_print_info(unit)

    !/****if* petsc/p_print_info
    !*  NAME
    !*    SUBROUTINE: p_print_info
    !*  SYNOPSIS
    !*    Usage:      p_print_info(unit)
    !*  FUNCTION
    !*    Prints the convergence information
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    unit               : Integer
    !*                         File unit to print to.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    29.05.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 06.06.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    INTEGER, INTENT(in) :: unit

    Vec         :: w,y
    KSP,POINTER :: ksp
    PetscScalar :: s

    ! Reduce typing
    ksp => p_object%ksp(p_object%solver)
    
    ! Preconditioned residual norm tolerance
    CALL KSPGetTolerances(ksp,p_object%rtol,PETSC_NULL_REAL,                   &
                          PETSC_NULL_REAL,PETSC_NULL_INTEGER,p_ierr)

    ! True residual L2 norm
    CALL VecDuplicate(p_object%b,y,p_ierr)
    CALL VecDuplicate(y,w,p_ierr)
    ! y = Ax
    CALL MatMult(p_object%A,p_object%x,y,p_ierr)
    ! w = b - y == b - Ax
    s = -1.0 ! VecWAXPY needs a PetscScalar
    CALL VecWAXPY(w,s,y,p_object%b,p_ierr)
    ! ||b - Ax||
    CALL VecNorm(w,NORM_2,p_object%r_n2,p_ierr)
    CALL VecDestroy(w,p_ierr)
    CALL VecDestroy(y,p_ierr)

    ! L2 norm of RHS
    ! ||b||
    CALL VecNorm(p_object%b,NORM_2,p_object%b_n2,p_ierr)
    
    CALL KSPGetIterationNumber(ksp,p_object%its,p_ierr)
    CALL KSPGetConvergedReason(ksp,p_object%reason,p_ierr)
    p_object%description = p_describe_reason(p_object%reason)
  
    IF (numpe == 1)THEN
      IF (.NOT. p_object%nsolvers_set) THEN
        WRITE(unit,'(A)') "Solver"
      ELSE
        WRITE(unit,'(A,I0)') "Solver ",p_object%solver
      END IF
      WRITE(unit,'(A,I0,A)')                                                   &
        "The reason for convergence was ",p_object%reason," "                  &
        //TRIM(p_object%description)
      WRITE(unit,'(A,I0)')                                                     &
        "The number of iterations to convergence was ",p_object%its
      WRITE(unit,'(A,E17.7)')                                                  &
        "The preconditioned relative error tolerance was ",p_object%rtol
      WRITE(unit,'(A,E17.7)')                                                  &
        "The relative error ||b-Ax||/||b|| was           ",                    &
        p_object%r_n2/p_object%b_n2
    END IF
  END SUBROUTINE p_print_info

  FUNCTION p_describe_reason(reason)

    !/****if* petsc/p_describe_reason
    !*  NAME
    !*    FUNCTION: p_describe_reason
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
    
    ! The string constants in this routine have to be string_length or
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

  FUNCTION p_memory_use()

    !/****if* petsc/p_memory_use
    !*  NAME
    !*    FUNCTION: p_memory_use
    !*  SYNOPSIS
    !*    Usage:      p_memory_use()
    !*  FUNCTION
    !*    Return amount of memory (in GB) in use in the program.  Call
    !*    p_release_memory just before p_memory_use to get the amount of memory
    !*    really in use:  otherwise you get (approximately) the peak memory
    !*    use since the last p_release_memory call.
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    27.06.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 27.06.2016, Mark Filipiak
    !*    Version 2, 10.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  Linux only.
    !*/

    REAL    :: p_memory_use

    INTEGER :: proc,ios,start,finish
    CHARACTER(len=string_length) :: line

    INTEGER(int64) :: kbytes
    ! real32 matches MPI_REAL4
    REAL(real32)   :: VmRSS,sum_VmRSS

    kbytes = 0

    OPEN(newunit=proc,file="/proc/self/status",action='read')
    DO
      READ(proc,"(A)",iostat=ios) line
      IF (ios < 0 ) EXIT ! end of file

      ! test for a line like "VmRSS:       972 kB"
      IF (INDEX(line,"VmRSS:") == 1 .AND. INDEX(line,"kB") /= 0) THEN
        start  = LEN("VmRSS:") + 1
        finish = INDEX(line,"kB") - 1
        READ(line(start:finish),*) kbytes
        EXIT
      END IF
    END DO
    CLOSE(proc)

    VmRSS = kbytes / 1048576.0 ! VmRSS in GB
    ! MPI_REAL4 matches real32
    CALL MPI_Allreduce(VmRSS,sum_VmRSS,1,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,     &
                       m_ierr)
    p_memory_use = sum_VmRSS
  END FUNCTION p_memory_use
  
  FUNCTION p_memory_peak()

    !/****if* petsc/p_memory_peak
    !*  NAME
    !*    FUNCTION: p_memory_peak
    !*  SYNOPSIS
    !*    Usage:      p_memory_peak()
    !*  FUNCTION
    !*    Return maximum amount of memory (in GB) used so far in the program.
    !*    Linux only.
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    10.08.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 10.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    REAL    :: p_memory_peak

    INTEGER :: proc,ios,start,finish
    CHARACTER(len=string_length) :: line

    INTEGER(int64) :: kbytes
    ! real32 matches MPI_REAL4
    REAL(real32)   :: VmHWM,sum_VmHWM

    kbytes = 0

    OPEN(newunit=proc,file="/proc/self/status",action='read')
    DO
      READ(proc,"(A)",iostat=ios) line
      IF (ios < 0 ) EXIT ! end of file

      ! test for a line like "VmHWM:       972 kB"
      IF (INDEX(line,"VmHWM:") == 1 .AND. INDEX(line,"kB") /= 0) THEN
        start  = LEN("VmHWM:") + 1
        finish = INDEX(line,"kB") - 1
        READ(line(start:finish),*) kbytes
        EXIT
      END IF
    END DO
    CLOSE(proc)

    VmHWM = kbytes / 1048576.0 ! VmHWM in GB
    ! MPI_REAL4 matches real32
    CALL MPI_Allreduce(VmHWM,sum_VmHWM,1,MPI_REAL4,MPI_SUM,MPI_COMM_WORLD,     &
                       m_ierr)
    p_memory_peak = sum_VmHWM
  END FUNCTION p_memory_peak
  
  SUBROUTINE p_release_memory

    !/****if* petsc/p_release_memory
    !*  NAME
    !*    SUBROUTINE: p_release_memory
    !*  SYNOPSIS
    !*    Usage:      p_release_memory
    !*  FUNCTION

    !*    Release all freed memory back to the system.  Maybe not all if there
    !*    are gaps in the heap.  Linux only.  Does not work for the Cray
    !*    compiler with -h system_alloc.
    !*
    !*    Releasing memory reduces the peak memory requirement (perhaps the
    !*    system allocator is better than the tcmalloc and glibc allocators).
    !*  ARGUMENTS
    !*    None.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    08.08.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 10.08.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  For Cray, the default malloc is tcmalloc.  Other compilers (Gnu, Intel)
    !*  use the malloc in glibc.  The #ifdef _CRAYFTN test can only tell if the
    !*  Cray compiler is used, not if -h system_alloc is used.

    USE, INTRINSIC :: iso_c_binding, ONLY: c_int,c_size_t

    INTERFACE
      SUBROUTINE MallocExtension_ReleaseFreeMemory()                           &
                   BIND(c,name="MallocExtension_ReleaseFreeMemory")
      END SUBROUTINE MallocExtension_ReleaseFreeMemory
      FUNCTION malloc_trim(pad) BIND(c,name="malloc_trim")
        USE, INTRINSIC :: iso_c_binding, ONLY: c_int,c_size_t
        INTEGER(c_int)           :: malloc_trim
        INTEGER(c_size_t), VALUE :: pad
      END FUNCTION malloc_trim
    END INTERFACE
    
    INTEGER(c_int) :: i

#ifdef _CRAYFTN
    CALL MallocExtension_ReleaseFreeMemory() ! Cray
#else
    i = malloc_trim(0_c_size_t) ! others
#endif
  END SUBROUTINE p_release_memory
END MODULE parafem_petsc
