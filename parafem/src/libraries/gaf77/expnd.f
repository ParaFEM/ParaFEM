c  *******************************************************************
c  *                                                                 *
c  *                     Subroutine Expnd                            *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, Jan. 20, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  converts an n x n matrix to an (n+k) x (n+k) matrix.
c
c  Converts an N x N matrix to an (N+K) x (N+K) matrix. It is assumed that
c  the original matrix is stored columnwise in a vector and that suffic-
c  ient room has been left for the larger matrix. Essentially this routine
c  creates an N x N submatrix within the larger matrix. Argument are as
c  follows;
c
c      A     the N x N matrix(vector) which needs to be put into a (N+K) x
c            (N+K) matrix(vector). Sufficient room must have been left in
c            A in the calling routine to allow it to grow.
c
c      N     the original dimension of A
c
c      K     the increase in size of A (to N+K)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      Subroutine expnd( A, N, K )
      dimension A(*)

      nn = N*N
      nk = nn + (N-1)*K

      do 20 i = 1,N-1
         do 10 j = 1,N
            A(nk) = A(nn)
            nk = nk - 1
            nn = nn - 1
  10     continue
         nk = nk - K
  20  continue

      return
      end
