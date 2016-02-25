c  ***********************************************************************
c  *                                                                     *
c  *                          Subroutine psolv                           *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the solution to the matrix equation [A]{x} = {b}
c           given the LU decomposition of PA where P is a permutation matrix.
c
c  This routine accepts as input the LU decomposition of [P][A] stored in place
c  (ie in [A]) (it is assumed that L is unit lower triangular), the vector
c  of row permutations, and the right-hand-side vector {b}, and computes the
c  solution {x} of [A]{x} = {b}. The solution is returned in the vector {x}.
c  [P] is a permutation matrix that arose out of the partial pivoting LU
c  decomposition. It is represented by the integer vector {indx} of final
c  row indices (ie if indx(2) = 3, then row 3 has been swapped with row 2
c  during the LU decomposition -- the swap did not actually take place,
c  only indices were swapped).
c
c  This routine performs no checks of the matrix [U] to ensure that diagonal
c  values are not zero. If any are, a division by zero error will occur.
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
c  indx   integer vector of length at least n containing the row indices of
c         the permuted matrix [A], that is the index indx[2] = 3 implies
c         that the third row of [A] would occupy the 2nd row of [P][A]. (input)
c
c     b   real vector of length at least n containing the right-hand-side
c         vector of [A]{x} = {b}. (input)
c
c     x   real vector of length at least n which on output will contain the
c         solution vector (in the correct order). (output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine psolv( A, ia, n, indx, b, x )
      dimension A(ia,*), b(*), x(*)
      integer indx(*)
c					forward substitution
      i1 = indx(1)
      x(1) = b(i1)
      do 20 i = 2, n
         ii = indx(i)
         s = b(ii) - A(ii,1)*x(1)
         do 10 j = 2, i-1
            s = s - A(ii,j)*x(j)
  10     continue
         x(i) = s
  20  continue
c					backward substitution
      nn   = indx(n)
      x(n) = x(n)/A(nn,n)
      do 40 i = n-1, 1, -1
         ii = indx(i)
         s  = x(i) - A(ii,n)*x(n)
         do 30 j = i+1, n-1
            s = s - A(ii,j)*x(j)
  30     continue
         x(i) = s/A(ii,i)
  40  continue

      return
      end
