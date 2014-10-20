MODULE MP_INTERFACE

  !/****h* /mp_interface
  !*  NAME
  !*    MODULE: mp_interface
  !*  SYNOPSIS
  !*    Usage:      USE mp_interface
  !*  FUNCTION
  !*    Contains subroutines required to initialize MPI at the beginning
  !*    of the calling program and finalize MPI at the end of the calling 
  !*    program
  !*    
  !*    Subroutine             Purpose
  !*
  !*    FIND_PE_PROCS          Find number of processors and own rank
  !*    SHUTDOWN               Terminate MPI and STOP the program
  !*  AUTHOR
  !*    L. Margetts
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  USE mpi_wrapper
  USE precision
  USE global_variables
  
  IMPLICIT NONE

  INCLUDE "mpif.h"

  CONTAINS

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE FIND_PE_PROCS(num,numprocs)

  !/****f* mp_interface/find_pe_procs
  !*  NAME
  !*    SUBROUTINE: find_pe_procs
  !*  SYNOPSIS
  !*    Usage:      CALL find_pe_procs(num,numprocs)
  !*  FUNCTION
  !*    Initialise MPI
  !*    Get the rank of the processes and the total number of processes
  !*  INPUTS
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    num                    : Integer
  !*                           : Process number
  !*
  !*    numprocs               : Integer
  !*                           : Number of processes
  !*  AUTHOR
  !*    M. Pettipher
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: num, numprocs

    CALL MPI_INIT(ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_INIT: errcode = ',ier
    END IF

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,num,ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_COMM_RANK: errcode = ',ier
    END IF

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ier)
    IF (ier /= MPI_SUCCESS) THEN
      WRITE(*,'(A,A,I5)')'Error in MPI_COMM_SIZE: errcode = ',ier
    END IF

    num = num + 1

  END SUBROUTINE FIND_PE_PROCS

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE SHUTDOWN()

    IMPLICIT NONE

    CALL MPI_FINALIZE(ier)
    STOP "ParaFEM: shutdown: the program terminated successfully"
    
    RETURN

  END SUBROUTINE SHUTDOWN

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

END MODULE MP_INTERFACE
