c  ***********************************************************************
c  *                                                                     *
c  *                            Function gselu                           *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  LU decompose a matrix A and return its determinant using naive
c           Gaussian Elimination.
c
c  This function reduces a matrix A into its (Doolittle) LU decomposition
c  `in place',
c
c      [L][U] = [A]    ===>   [A] = [L\U]
c
c  No scaling or pivoting is performed during the decomposition.
c  As a by-product of the decomposition, the determinant is returned in
c  the function name. If the matrix [A] is found to be algorithmically
c  singular, then the determinant is set to zero and the function exited.
c  The matrix [L] is assumed to be unit lower triangular (ie., having 1's
c  on the diagonal) -- thus only the subdiagonal terms of [L] are stored
c  in [A]. Arguments to the routine are as follows;
c
c     A    real array of size at least n x n containing on input the
c          matrix of coefficients. On ouput, A will contain its own LU
c          decomposition. (input/output)
c
c    ia    leading dimension of A exactly as specified in the calling routine.
c          (input)
c
c     n    integer giving the number of equations in the problem (ie, the
c          number of rows in A). (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      real function gselu( A, ia, n )
      dimension A(ia,*)
      data zero/0.0/, one/1.0/
c					forward reduce A
      gselu = one
      do 20 k = 1, n - 1
         gselu = gselu*A(k,k)
         if( A(k,k) .eq. zero ) return
c						compute L(i,k)
         kp1 = k + 1
         do 10 i = kp1, n
            A(i,k) = A(i,k)/A(k,k)
  10     continue
c						reduce the rest of A
         do 20 j = kp1, n
         do 20 i = kp1, n
            A(i,j) = A(i,j) - A(i,k)*A(k,j)
  20  continue
c					include the last diagonal element
      gselu = gselu*A(n,n)

      return
      end
