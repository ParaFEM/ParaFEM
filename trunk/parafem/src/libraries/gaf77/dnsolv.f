c  ***********************************************************************
c  *                                                                     *
c  *                          Subroutine dnsolv                          *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the solution to the matrix equation [A]{x} = {b}
c           given the LU decomposition of A (Naive G.E. Version)
c
c  This routine accepts as input the LU decomposition of [A] stored in place
c  (ie in [A]) (it is assumed that L is unit lower triangular) and the
c  right-hand-side vector {b}, and computes the solution {x} of [A]{x} = {b}.
c  The solution is returned in the vector {b}.
c
c  This routine performs no checks of the matrix [U] to ensure that diagonal
c  values are not zero. If any are, a division by zero error will occur. It
c  is assumed that this check was performed in the LU decomposition stage.
c
c  Arguments to the routine are as follows;
c
c     A   real array of size at least n x n which on input contains the
c         LU decomposition of the row-permuted matrix [P][A]. [L] is
c         assumed to be stored in the lower sub-triangle of [A] and [U]
c         is stored in the upper triangle of [A]. (input)
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
      subroutine dnsolv( A, ia, n, b )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*), b(*)
c					forward substitution
      do 10 j = 1, n
      do 10 i = j+1, n
         b(i) = b(i) - A(i,j)*b(j)
  10  continue
c					backward substitution
      do 20 j = n, 1, -1
         b(j) = b(j)/A(j,j)
         do 20 i = 1, j-1
            b(i) = b(i) - A(i,j)*b(j)
  20  continue

      return
      end
