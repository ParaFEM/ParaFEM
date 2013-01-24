MODULE BICG

  !/****h* /bicg
  !*  NAME
  !*    MODULE: bicg
  !*  SYNOPSIS
  !*    Usage:      USE bicg
  !*  FUNCTION
  !*    Contains subroutines required by BiCGSTAB(l)
  !*    These subroutines are parallel and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    FORM_S                 Forms the s vector in bicgstab(l)
  !*  AUTHOR
  !*    L. Margetts
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2013 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  Obsolete module. Moved content to new_library.f90
  !*/
  
  USE precision
  USE maths
  USE gather_scatter

  IMPLICIT NONE
  
  CONTAINS


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

END MODULE BICG
