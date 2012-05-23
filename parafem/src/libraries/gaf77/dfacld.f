c  ***********************************************************************
c  *                                                                     *
c  *                            Function dfacld                          *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 2.11
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the LDL' factorization of a symmetric matrix A and
c           return its determinant.
c
c  This function reduces a matrix A into its LDL' decomposition (where
c  prime indicates transpose) in place. No scaling or pivoting is performed
c  during the decomposition. If the matrix A is found to be algorithmically
c  singular, then the determinant is set to zero and returned. The matrix [L]
c  is assumed to be unit lower triangular (ie., having 1's on the diagonal)
c  and [D] is a diagonal matrix.
c  Arguments to the routine are as follows;
c
c     A    real array of size at least n x n containing on input the
c          matrix of coefficients in its upper triangle at least. On
c          ouput, A will contain L written in the lower sub-triangle
c          and the RECIPROCALS of D written in the diagonal. (input/output)
c
c    ia    leading dimension of A exactly as specified in the calling routine.
c          (input)
c
c     n    integer giving the number of equations in the problem (ie, the
c          number of rows in A). (input)
c
c  REVISION HISTORY:
c  2.1	made function itself real*8
c  2.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      real*8 function dfacld( A, ia, n )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*)
      data zero/0.d0/, one/1.d0/
c					find elements of L
      dfacld = one
      do 40 i = 1, n
         do 20 j = 1, i-1
            s = A(j,i)
            do 10 k = j-1, 1, -1
               s = s - A(i,k)*A(j,k)
  10        continue
            A(i,j) = s
  20     continue
c					find diagonal terms of D
         s = A(i,i)
         do 30 k = i-1, 1, -1
            t = A(i,k)
            A(i,k) = t*A(k,k)
            s = s - t*A(i,k)
  30     continue
c						update determinant
         dfacld = s*dfacld
         if( s .eq. zero ) return
c						and put 1/D in A
         A(i,i) = one/s
  40  continue

      return
      end
