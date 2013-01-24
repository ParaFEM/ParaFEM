MODULE PLASTICITY

  !------------------------------------------------------------------------------
  !/****h* ParaFEM-Shared/plasticity
  !*  NAME
  !*    plasticity module
  !*  FUNCTION
  !*    Contains the following subroutines and functions needed to solve 
  !*    plasticity problems:
  !*
  !*    Subroutine             Purpose
  !*
  !*    INVAR                  Returns stress invariants SIGMA, DSBAR and 
  !*                           Lode angle THETA from current stresses held in 
  !*                           STRESS.
  !*
  !*    VMFLOW                 Returns the von Mises flow vector VMFL. STRESS 
  !*                           holds the second deviatoric invariant (2D only).
  !*
  !*    FMACAT                 Intermediate step.
  !*
  !*    FMRMAT                 Returns matrix RMAT from the von Mises flow 
  !*                           vector VMFL, invariant DSBAR, plastic multiplier
  !*                           DLAM and elastic matrix DEE.
  !*
  !*    FORMAA                 Returns modified matrix DAATD from the von Mises
  !*                           flow vector VMFL and matrix RMAT.
  !*
  !*    VMPL                   Returns the plastic stress-strain matrix PL for
  !*                           a von Mises material. STRESS holds the stresses
  !*                           and DEE holds the elastic stress-strain matrix.
  !*
  !*    FORM_TEMP              Forms temp matrix for consistent return.
  !*
  !*    MOCOUF                 Returns the Mohr-Coulomb failure function, F, 
  !*                           from the strength parameters PHI and C and 
  !*                           stress invariants SIGM, DSBAR and THETA.
  !*
  !*    MOCOUQ                 Returns the plastic potentail terms DQ1, DQ2
  !*                           and DQ3 for a Mohr-Coulomb material from
  !*                           dilation angle PSI (in degrees) and invariants
  !*                           DSBAR and THETA.
  !*
  !*  AUTHOR
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !****
  !*
  !* Module obsolete
  !* Content moved to new_library
  !*
  !*/

  USE precision
  
  CONTAINS

END MODULE PLASTICITY
