MODULE ELEMENTS

  !/****h* /elements
  !*  NAME
  !*    MODULE: elements
  !*  SYNOPSIS
  !*    Usage:      USE elements
  !*  FUNCTION
  !*    Contains subroutines used for computing element level matrices
  !*    
  !*    Subroutine             Purpose
  !*
  !*    SHAPE_FUN              Computes the shape functions
  !*    SHAPE_DER              Compute the derivatives of the shape functions
  !*    BEEMAT                 Compute the B matrix
  !*    SAMPLE                 Returns local coords of the integrating points
  !*    ECMAT                  Returns the consistent mass matrix
  !*    DEEMAT                 Compute the D matrix
  !*    SHAPE_DERIVATIVES      Compute the derivatives of the shape functions
  !* 
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2004-2011 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  
  !*  Module is now obsolete and earmarked for deletion
  !*  LM has moved all subroutines to MODULE NEW_LIBRARY.F90
  !*
  !*  To update your code, replace "USE elements" with "USE new_library"
  !*
  !*  Warning: Argument order in DEEMAT is now (dee,e,v)
  !*  Warning: Argument order in BEEMAT is now (bee,deriv)
  !*/
  
  USE precision

  CONTAINS                                  


  
END MODULE ELEMENTS
