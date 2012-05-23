c  ***********************************************************************
c  *                                                                     *
c  *                            Function dcabs                           *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute and return the real*8 modulus of a complex*16 value
c
c  This routine accepts as input a complex*16 argument `c' (interpreted as a
c  pair of real*8 values: (real,imaginary) components), and computes the
c  modulus of the complex number.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      real*8 function dcabs( c )
      real*8 c(*), dsqrt

      dcabs = dsqrt( c(1)*c(1) + c(2)*c(2) )

      return
      end
