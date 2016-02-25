c  *****************************************************************
c  *                                                               *
c  *                  Subroutine dsmslv                            *
c  *                                                               *
c  *****************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, Jan. 20, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  solve [A]{x} = {b} given [A-inverse] for [A] symmetric
c
c  Solves Ax = b for x given A-inverse stored in the lower triangle of A.
c  Arguments are as follows;
c
c     A    real array of size at least n x n which contains the elements
c          of [A-inverse] in its lower triangle. (input)
c
c    ia    leading dimension of A exactly as specified in the calling routine.
c          (input)
c
c     b    real vector of length at least n which contains the RHS elements.
c          (input)
c
c     x    real vector of length at least n which will contain the solution.
c          (output)
c
c     n    number of equations. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine dsmslv( A, ia, b, x, n )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*), b(*), x(*)
      data zero/0.d0/

      do 20 i = 1, n
         x(i) = zero
         do 10 j = 1,i
            x(i) = x(i) + A(i,j) * b(j)
  10     continue
         do 20 j = i+1, n
            x(i) = x(i) + A(j,i) * b(j)
  20  continue

      return
      end
