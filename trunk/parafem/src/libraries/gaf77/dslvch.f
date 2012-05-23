c  ***********************************************************************
c  *                                                                     *
c  *                          Subroutine dslvch                          *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the solution to the matrix equation [A]{x} = {b}
c           given the Cholesky LL' decomposition of A as produced by DCHOLL.
c
c  This routine accepts as input the Cholesky LL' decomposition of [A] stored
c  in place (ie in [A], see below) and the right-hand-side vector {b}, and
c  computes the solution {x} of [A]{x} = {b}. The solution is returned in the
c  vector {b}. Arguments to the routine are as follows;
c
c     A   real array of size at least n x n which on input contains the
c         Cholesky LL' decomposition of the matrix [A]. [L] is
c         assumed to be stored in the lower triangle of A. (input)
c
c    ia   integer giving the leading dimension of A exactly as specified
c         in the calling routine. (input)
c
c     n   integer giving the size of the matrix A. (input)
c
c     b   real vector of length at least n containing the right-hand-side
c         vector of [A]{x} = {b}. On output, b will contain the solution.
c         (input/output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine dslvch( A, ia, n, b )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*), b(*)
c					forward substitution (solve Ly = b)
      do 10 j = 1, n-1
         b(j) = b(j)/a(j,j)
         do 10 i = j+1, n
            b(i) = b(i) - A(i,j)*b(j)
  10  continue
      b(n) = b(n)/(A(n,n)*A(n,n))
c					back-sub (solve L'x = y)
      do 30 i = n-1, 1, -1
         do 20 j = i+1, n
            b(i) = b(i) - A(j,i)*b(j)
  20     continue
         b(i) = b(i)/A(i,i)
  30  continue

      return
      end
