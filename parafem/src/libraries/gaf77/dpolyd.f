c  *********************************************************************
c  *                                                                   *
c  *                          Function dpolyd                          *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to return the result of an (n-1)th order polynomial in Newton's
c            divided difference format evaluated at the point z.
c
c  This function accepts a vector of data locations {x_i}, i = 1, 2,..., n,
c  a vector of coefficients {b_i}, i = 1, 2, ..., n, and a point z, and
c  evaluates the polynomial
c
c   f_{n-1}(z) = b_1 + b_2(z - x_1) + b_3(z - x_1)(z - x_2) + ....
c                    + b_n(z - x_1)(z - x_2)...(z - x_{n-1})
c
c  The result is returned in the function's name. Arguments to the function
c  are as follows;
c
c     z   real value at which the polynomial is to be evaluated. (input)
c
c     x   real vector of length at least n which contains the data
c         locations {x_i}, i = 1, 2, ..., n. (input)
c
c     b   real vector of length at least n which contains the coefficients
c         of the Newton Divided Difference polynomial. (input)
c
c     n   the number of data points. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      real*8 function dpolyd( z, x, b, n )
      implicit real*8 (a-h,o-z)
      dimension x(*), b(*)

      dpolyd = b(n)
      do 10 i = n-1, 1, -1
         dpolyd = b(i) + (z - x(i))*dpolyd
  10  continue

      return
      end
