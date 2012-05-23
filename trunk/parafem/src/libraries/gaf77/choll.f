c  *****************************************************************
c  *                                                               *
c  *                      Function choll                           *
c  *                                                               *
c  *****************************************************************
c  Single Precision Version 2.01, Matrix storage
c  Written by Gordon A. Fenton, TUNS, Feb. 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE to compute the Cholesky decomposition L according to L*L' = A
c          and return the determinant of A
c
c  Computes the Cholesky decomposition L of a symmetric positive definite
c  real matrix A = L*L' where the prime indicates transpose. L is a real
c  non-singular lower-triangular matrix. This method involves N divisions,
c  N square-roots, and (0.16*N**3 + 0.5*n**2) multiplications. If A is
c  discovered to be non-positive definite or singular, the determinant is
c  set to 0 and the function returned.
c  The variable description is as follows;
c
c    A    input positive definite symmetric matrix with necessary values
c         stored in the lower triangle (at least). On output, A will contain
c         the Cholesky decomposition in its lower triangle (including the
c         diagonal). The upper triangle of A is not touched.
c
c    IA   input integer containing the column dimension of A (and thus L) as
c         specified in the calling routine
c
c    N    input integer containing the order of matrix A.
c
c  REVISION HISTORY:
c  2.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c ===========================================================================
      real function choll( a, ia, n )
      dimension a(ia,*)
      data zero/0.0/, one/1.0/

      det   = one
      choll = zero
      do 20 k = 1, n
c						check diagonal element
         if( a(k,k) .le. zero ) return
c						find L and determinant
         det = det*a(k,k)
         a(k,k) = sqrt( a(k,k) )
         do 10 j = k+1, n
            a(j,k) = a(j,k)/a(k,k)
  10     continue

         do 20 j = k+1, n
         do 20 i = j, n
            a(i,j) = a(i,j) - a(i,k)*a(j,k)
  20  continue

      choll = det
      return
      end
