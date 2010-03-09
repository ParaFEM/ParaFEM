MODULE precision

  !/****h* modules/precision
  !*  NAME
  !*    MODULE: precision
  !*  SYNOPSIS
  !*    Usage:      USE precision
  !*  FUNCTION
  !*    Sets IWP 
  !*  AUTHOR
  !*    M.A. Pettipher
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  INTEGER, PARAMETER :: iwp = SELECTED_REAL_KIND(15,300)

END MODULE precision
