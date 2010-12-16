MODULE MP_INTERFACE

  !/****h* /mp_interface
  !*  NAME
  !*    MODULE: mp_interface
  !*  SYNOPSIS
  !*    Usage:      USE mp_interface
  !*  FUNCTION
  !*    Contains subroutines required to initialize MPI at the beginning
  !*    of the calling program and finalize MPI at the end of the calling 
  !*    program. Dummy serial version.
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
  
  USE precision
  USE global_variables
  
  IMPLICIT NONE

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

    num      = 1
    numprocs = 1

  END SUBROUTINE FIND_PE_PROCS

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE SHUTDOWN()

    IMPLICIT NONE

    STOP
    
    RETURN

  END SUBROUTINE SHUTDOWN

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

END MODULE MP_INTERFACE
