c  *****************************************************************
c  *                                                               *
c  *                  Subroutine Factr                             *
c  *                                                               *
c  *****************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, Jan. 20, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  factorizes a real symmetric submatrix A(n1:n2) into L*D*L'
c
c  Factorize a real positive definite symmetric matrix A into a lower
c  triangular and diagonal matrix where A = L*D*L'. L is a lower
c  triangular unit matrix (ones on the diagonal) and D is a diagonal
c  matrix.  This method is believed superior to the Cholesky factor-
c  ization method because the square roots have been eliminated. The
c  number of other operations is unchanged. The upper triangle of A is
c  left unchanged. On output the lower triangle contains the elements of
c  the matrix L and the diagonal contains the reciprocals of the elements
c  of D.
c  Variables are described as follows;
c
c      A    an square real symmetric matrix with its values stored in the
c           upper triangle at least. On output, the factorized values are
c           overwritten on the lower triangular part of A. Note that the
c           original diagonal elements of A are lost so if A is desired
c           later, the diagonal values should be stored seperately.
c
c      IA   input integer denoting the column size of A as specified in the
c           calling routine
c
c      n1,n2  starting and ending diagonal locations defining the submatrix
c             to be factorized, ie. if the entire matrix is to be factorized,
c             then n1 = 1 and n2 = n
c
c      ierr   flag set to zero if factorization successfull and to -1 if A
c             is algorithmically singular
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      Subroutine factr( A, IA, n1, n2, ierr )
      dimension A(ia,*)
      data zero/0.0/, one/1.0/

      ierr = -1
      do 40 i = n1,n2
         do 20 j = n1, i-1
            x = A(j,i)
            do 10 k = j-1, n1, -1
               x = x - A(i,k) * A(j,k)
  10        continue
            A(i,j) = x
  20     continue

         x = A(i,i)
         do 30 k = i-1, n1, -1
            y = A(i,k)
            A(i,k) = y * A(k,k)
            x = x - y * A(i,k)
  30     continue
         if ( x .eq. zero ) return
         A(i,i) = one/x
  40  continue
c
      ierr = 0
      return
      end
