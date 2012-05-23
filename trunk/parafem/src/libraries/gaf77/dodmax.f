C  **********************************************************************
C  *                                                                    *
C  *                          Function dodmax                           *
C  *                                                                    *
C  **********************************************************************
C  Double Precision Version 1.01
C  Written by Gordon A. Fenton, Princeton, 1990
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  returns the maximum off-diagonal element of a square matrix
C
C  Returns the maximum off-diagonal value of a square matrix.
C  Arguments are as follows;
C
C    A   real array of size N x N containing the matrix to consider.
C        (input)
C
C    N   the size of the matrix A. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
C--------------------------------------------------------------------------
      real*8 function dodmax( A, n )
      real*8 A(n,*)
c                                      initialize dodmax
      dodmax = A(2,1)
c                                      search for larger off-diagonal values
      do 10 i = 3, n
         if( dodmax .lt. A(i,1) ) dodmax = A(i,1)
  10  continue

      do 30 j = 2, n-1
         do 20 i = 1, j-1
            if( dodmax .lt. A(i,j) ) dodmax = A(i,j)
  20     continue
         do 30 i = j+1, n
            if( dodmax .lt. A(i,j) ) dodmax = A(i,j)
  30  continue

      do 40 i = 1, n-1
         if( dodmax .lt. A(i,n) ) dodmax = A(i,n)
  40  continue

      return
      end
