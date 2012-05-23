c  **********************************************************************
c  *                                                                    *
c  *                          Function odmax                            *
c  *                                                                    *
c  **********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, 1990
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  returns the maximum off-diagonal element of a square matrix
c
c  Returns the maximum off-diagonal value of a square matrix.
c  Arguments are as follows;
c
c    A   real array of size N x N containing the matrix to consider.
c        (input)
c
c    N   the size of the matrix A. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c--------------------------------------------------------------------------
      real function odmax( A, n )
      real A(n,*)
c                                      initialize odmax
      odmax = A(2,1)
c                                      search for larger off-diagonal values
      do 10 i = 3, n
         if( odmax .lt. A(i,1) ) odmax = A(i,1)
  10  continue

      do 30 j = 2, n-1
         do 20 i = 1, j-1
            if( odmax .lt. A(i,j) ) odmax = A(i,j)
  20     continue
         do 30 i = j+1, n
            if( odmax .lt. A(i,j) ) odmax = A(i,j)
  30  continue

      do 40 i = 1, n-1
         if( odmax .lt. A(i,n) ) odmax = A(i,n)
  40  continue

      return
      end
