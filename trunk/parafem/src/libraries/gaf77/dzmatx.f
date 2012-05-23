c   ********************************************************************
c   *                                                                  *
c   *                      Subroutine dzmatx                           *
c   *                                                                  *
c   ********************************************************************
c   Double Precision Version 1.01
c   Written by Gordon A. Fenton, Princeton, Jan. 20, 1988.
c  Latest Update: Jun 9, 1999
c
c   PURPOSE zeroes the extra rows and columns created by DEXPND.
c
c   Zeroes the new locations created by DEXPND (the last K rows and columns)
c   Arguments are as follows;
c
c     A     real array of size at least NK x NK which is the matrix to be
c           partially zeroed. (input/output)
c
c    NK     row and column dimension of A as produced by DEXPND. (input)
c
c     N     row and column dimension of the original (non-zero) matrix
c           A prior to modification by DEXPND. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      Subroutine dzmatx( A, NK, N )
      implicit real*8 (a-h,o-z)
      dimension A(NK,*)
      data zero/0.d0/
c
      do 10 j = 1,N
      do 10 i = N+1, NK
         A(i,j) = zero
         A(j,i) = zero
  10  continue
c
      do 20 j = N+1, NK
         A(j,j) = zero
         do 20 i = N+1, j-1
            A(i,j) = zero
            A(j,i) = zero
  20  continue
c
      return
      end
