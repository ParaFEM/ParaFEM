c  ***********************************************************************
c  *                                                                     *
c  *                             Function cabs                           *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute and return the real*4 modulus of a complex*8 value
c
c  This routine accepts as input a complex*8 argument `c' (interpreted as a
c  pair of real*4 values: (real,imaginary) components), and computes the
c  modulus of the complex number.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      real function cabs( c )
      real c(*)

      cabs = sqrt( c(1)*c(1) + c(2)*c(2) )

      return
      end
