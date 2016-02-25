c  ********************************************************************
c  *                                                                  *
c  *                          Function dtrilu                         *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to compute the LU decomposition of a tridiagonal matrix
c            and return its determinant.
c
c  This routine takes a tridiagonal matrix input as a set of three
c  vectors and computes its LU decomposition, where L is unit-lower
c  bidiagonl and U is upper bidiagonal. The determinant is computed
c  as the product of the diagonal terms of U and is set to zero if
c  any of the diagonal terms are zero (in which case the factorization
c  is discontinued). The matrix is assumed to be stored in the form
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
c          super-diagonal elements with eu(1) = 0. On ouput eu contains
c          the superdiagonal elements of U. (input/output)
c
c     d    real vector of length at least n which on input contains the
c          diagonal elements. On output, d contains the diagonal elements
c          of U. (input/output)
c
c    el    real vector of length at least n which on input contains the
c          subdiagonal elements. On output, el contains the subdiagonal
c          elements of L. (input/output)
c
c     n    size of the tridiagonal matrix. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      real*8 function dtrilu( eu, d, el, n )
      implicit real*8 (a-h,o-z)
      dimension eu(*), d(*), el(*)
      data zero/0.d0/

      dtrilu = d(1)
      do 10 k = 2, n
         if( d(k-1) .eq. zero ) return
         el(k) = el(k)/d(k-1)
         d(k)  = d(k) - el(k)*eu(k)
         dtrilu = dtrilu*d(k)
  10  continue

      return
      end
