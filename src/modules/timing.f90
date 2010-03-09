MODULE timing

  !/****h* /timing
  !*  NAME
  !*    MODULE: timing
  !*  SYNOPSIS
  !*    Usage:      USE timing
  !*  FUNCTION
  !*    Contains subroutines used for timing purposes. These require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    ELAP_TIME              Computes wall clock time
  !*  AUTHOR
  !*    L. Margetts
  !*    M. Pettipher
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  USE precision
  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  FUNCTION elap_time()
  REAL(iwp) elap_time

  !/****f* timing/elap_time
  !*  NAME
  !*    SUBROUTINE: elap_time
  !*  SYNOPSIS
  !*    Usage:      CALL elap_time()
  !*  FUNCTION
  !*    Returns the wall clock time
  !*  AUTHOR
  !*    M. Pettipher
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
   
 
  INTEGER :: count,count_rate
 
  CALL SYSTEM_CLOCK(count,count_rate)
  elap_time = REAL(count,iwp)/REAL(count_rate,iwp)
 
 
  END FUNCTION elap_time

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END MODULE timing
