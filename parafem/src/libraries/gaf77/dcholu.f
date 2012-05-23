c  *****************************************************************
c  *                                                               *
c  *                     Function dcholu                           *
c  *                                                               *
c  *****************************************************************
c  Double precision version 2.01, Matrix storage
c  Written by G. A. Fenton, TUNS, Feb. 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE to compute the Cholesky decomposition U from U'*U = A
c          and return the determinant of A
c
c  Computes the Cholesky decomposition U of a symmetric positive definite
c  real matrix A = U'*U where the prime indicates transpose. U is a real
c  non-singular upper-triangular matrix and U' is its tranpose. Thus this
c  routine returns the upper triangular matrix U overwritten in A. The
c  method involves N divisions, N square-roots, and (0.16*N**3 + 0.5*n**2)
c  multiplications. If A is discovered to be non-positive definite (or
c  singular), then the determinant is set to zero and the function returned.
c  Arguments to the routine are as follows;
c
c    A    input positive definite symmetric matrix with necessary values
c         stored in the upper triangle (at least). On output, the upper
c         triangle of A will be overwritten with U. The lower triangle
c         of A is not touched. (input/output)
c
c    IA   input integer containing the column dimension of A as
c         specified in the calling routine
c
c    N    input integer containing the order of matrix A.
c
c  REVISION HISTORY:
c  2.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c =========================================================================
      real*8 function dcholu( a, ia, n )
      implicit real*8 (a-h,o-z)
      dimension a(ia,*)
      data zero/0.d0/, one/1.d0/
      sqrt(x) = dsqrt(x)

      det    = one
      dcholu = zero
      do 10 k = 1, n
c						check diagonal elements
         if( a(k,k) .le. zero ) return
c						update determinant and find U
         det    = det*a(k,k)
         a(k,k) = sqrt( a(k,k) )
         do 10 j = k+1, n
            a(k,j) = a(k,j)/a(k,k)
            do 10 i = k+1, j
               a(i,j) = a(i,j) - a(k,i)*a(k,j)
  10  continue

      dcholu = det
      return
      end
