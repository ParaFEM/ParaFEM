c  ********************************************************************
c  *                                                                  *
c  *                       Subroutine dtrisv                          *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 1.11
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to solve a tridiagonal system of equations given the LU
c            decomposition of the tridiagonal coefficient matrix
c
c  This routine takes the LU decomposition of a tridiagonal matrix and
c  computes the solution to the system of equations [A]{x} = {b}. The
c  LU decomposition is assumed to be stored in three vectors and L is
c  assumed to be unit-lower bidiagonal. See `trilu' for the decomposition
c  routine. The matrix is assumed to be stored in the form
c
c            _                                      _
c           |  d(1)  eu(2)   0    0     0   ...   0  |
c           |  el(2) d(2)  eu(3)  0     0   ...   0  |
c           |   0    el(3) d(3) eu(4)   0   ...   0  |
c           |   0     0    el(4) d(4) eu(5) ...   0  |
c           |   .     .        .    .     .       .  |
c           |   .     .           .    .     .    .  |
c           |   .     .              .    .    .  .  |
c           |   .     .                 .    .  eu(n)|
c           |   0     0    . . . . . .   el(n)  d(n) |
c            -                                      -
c  Arguments to the function are as follows;
c
c    eu    real vector of length at least n which on input contains the
c          super-diagonal elements of U with eu(1) = 0. (input)
c
c     d    real vector of length at least n which on input contains the
c          diagonal elements of U. This routine performs no checks to
c          ensure that elements of d are non-zero. This should have been
c          checked in the LU decomposition stage. (input)
c
c    el    real vector of length at least n which on input contains the
c          subdiagonal elements of L (diagonal elements of L are all 1's).
c          (input)
c
c     b    real vector of length at least n which on input contains the
c          RHS vector. On output b contains the solution
c
c     n    size of the tridiagonal matrix. (input)
c
c  REVISION HISTORY:
c  1.1	eliminated unused local variable `zero' (Dec 5/96)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine dtrisv( eu, d, el, b, n )
      implicit real*8 (a-h,o-z)
      dimension eu(*), d(*), el(*), b(*)

      do 10 k = 2, n
         b(k) = b(k) - el(k)*b(k-1)
  10  continue
      b(n) = b(n)/d(n)
      do 20 k = n-1, 1, -1
         b(k) = (b(k) - eu(k+1)*b(k+1))/d(k)
  20  continue

      return
      end
