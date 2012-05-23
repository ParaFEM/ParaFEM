c  *********************************************************************
c  *                                                                   *
c  *                        Subroutine dntddm                          *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to compute the coefficients of an (n-1)th order polynomial
c            via Newton's Divided Difference algorithm.
c
c  This routine accepts a vector of data points (x_i, f(x_i)), for
c  i = 1, 2, ..., n, and computes the unique coefficients of the
c  interpolating polynomial
c
c   f_{n-1}(x) = b_1 + b_2(x - x_1) + b_3(x - x_1)(x - x_2) + ....
c                    + b_n(x - x_1)(x - x_2)...(x - x_{n-1})
c
c  which passes through all the data points. Arguments to the routine
c  are as follows;
c
c     x   real vector of length at least n which contains the data
c         locations {x_i}, i = 1, 2, ..., n. (input)
c
c     f   real vector of length at least n which on input contains the
c         observed function values at each x_i, i = 1, 2, ..., n. On
c         output, {f} will contain the set of coefficients {b_i}, i = 1,
c         2, ..., n. (input/output)
c
c     n   the number of data points to be fitted. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine dntddm( x, f, n )
      implicit real*8 (a-h,o-z)
      dimension x(*), f(*)

      do 10 i = 2, n
      do 10 j = n, i, -1
         f(j) = (f(j) - f(j-1))/(x(j) - x(j-i+1))
  10  continue

      return
      end
