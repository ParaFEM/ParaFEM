c  ********************************************************************
c  *                                                                  *
c  *                       Function dpolyn                            *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to return the result of an (n-1)'th order polynomial f(x)
c
c  This function evaluates the (n-1)'th order polynomial
c
c     f(x) = a_1 + a_2*x + a_3*x^2 + . . . + a_n*x^{n-1}
c
c  by Horner's rule at the point x and returns the result.
c  Arguments to the function are;
c
c     x    the point at which the function is to be evaluated. (input)
c
c     a    real vector of length at least n which contains the coefficients
c          of the polynomial. (input)
c
c     n    the number of coefficients in the vector {a}. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      real*8 function dpolyn( x, a, n )
      implicit real*8 (a-h,o-z)
      dimension a(*)

      dpolyn = a(n)
      do 10 i = n-1, 1, -1
         dpolyn = a(i) + dpolyn*x
  10  continue

      return
      end
