c  ***********************************************************************
c  *                                                                     *
c  *                          Subroutine dslvld                          *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the solution to the matrix equation [A]{x} = {b}
c           given the LDL' decomposition of A as produced by DFACLD.
c
c  This routine accepts as input the LDL' decomposition of [A] stored in place
c  (ie in [A], see below) and the right-hand-side vector {b}, and computes
c  the solution {x} of [A]{x} = {b}. The solution is returned in the
c  vector {b}. Arguments to the routine are as follows;
c
c     A   real array of size at least n x n which on input contains the
c         LDL' decomposition of the matrix [A]. [L] is
c         assumed to be stored in the lower sub-triangle of A and the
c         reciprocals of the elements of [D] in the diagonal of A, as
c         produced by DFACLD. (input)
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
      subroutine dslvld( A, ia, n, b )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*), b(*)
c					forward substitution (solve Ly = b)
      do 10 j = 1, n
      do 10 i = j+1, n
         b(i) = b(i) - A(i,j)*b(j)
  10  continue
c					back-sub (solve Dz = y and L'x = z)
      do 20 i = n, 1, -1
         b(i) = b(i)*A(i,i)
         do 20 j = i+1, n
            b(i) = b(i) - A(j,i)*b(j)
  20  continue

      return
      end
