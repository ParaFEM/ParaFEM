MODULE PCG

  !/****h* /pcg
  !*  NAME
  !*    MODULE: pcg
  !*  SYNOPSIS
  !*    Usage:      USE pcg
  !*  FUNCTION
  !*    Contains subroutines required by the preconditioned conjugate gradient
  !*    method. These subroutines are parallel and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    CHECON_PAR             Convergence test
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
  USE maths
  USE gather_scatter

  IMPLICIT NONE
  
  CONTAINS


END MODULE PCG
