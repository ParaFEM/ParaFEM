MODULE global_variables

  !/****h* modules/global_variables
  !*  NAME
  !*    MODULE: global_variables
  !*  SYNOPSIS
  !*    Usage:      USE global_variables
  !*  FUNCTION
  !*    Contains variables used in other modules and the main programs
  !*  AUTHOR
  !*    M. Pettipher   
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010 
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  USE precision

  IMPLICIT NONE

! The following may be used in the main program and gather_scatter

  INTEGER, PARAMETER   :: details = 6  ! Use standard output for details
  INTEGER              :: ntot, neq, nels_pp, neq_pp, numpe

! Variables for packing data to be broadcast. (utility and gather_scatter)

  INTEGER              :: position, bufsize, recbufsize, ier
  INTEGER, ALLOCATABLE :: tempbuf(:)

END MODULE global_variables
