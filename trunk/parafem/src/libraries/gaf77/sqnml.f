c  *******************************************************************
c  *                                                                 *
c  *                          Function sqnml                         *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to solve the general least squares problem via the normal
c           equations [A^t][A]{c} = [A^t]{b} and return the residual
c           sum of squared errors.
c
c  This routine accepts an n x m matrix [A] (with n > m-1) and computes
c  the set of coefficients {c} which minimize the sum of squared errors
c  (using the 2-norm);
c
c                             2
c           || {b} - [A]{c} ||
c
c  This is accomplished by forming [W] = [A^t][A], {r} = [A^t]{b} and
c  solving
c
c           [W]{c} = {r}
c
c  using an LDL' factorization/solver ([W] is symmetric and positive
c  definite). Thus this routine calls the function FACLD(A,ia,n) and
c  the subroutine SLVLD(A,ia,n,b). If the matrix [W] is found to be
c  singular, then the sum of squared errors is set to -1 and returned.
c  Arguments to this routine are as follows;
c
c    A    real array of size at least n x m containing the elements of
c         the array [A]. (input)
c
c   ia    leading dimension of the array A exactly as given in the dimension
c         statement of the calling routine. (input)
c
c    W    real array of size at least m x m which is used to store [A^t][A]
c         and then its LDL' factorization. (output)
c
c   iw    leading dimension of the array W exactly as given in the dimension
c         statement of the calling routine. (input)
c
c    b    real vector of length at least n which contains the RHS. (input)
c
c    c    real vector of length at least m which will contain the solution.
c         (output)
c
c    n    integer column dimension of the matrix [A] (n > m-1). (input)
c
c    m    row dimension of the matrix [A]. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      real function sqnml( A, ia, W, iw, b, c, n, m )
      dimension A(ia,*), W(iw,*), b(*), c(*)
      data zero/0.0/, one/1.0/

c					form [W] = [A^t][A] and {r} = [A^t]{b}
      do 30 j = 1, m
         do 20 i = 1, j
            s = A(1,i)*A(1,j)
            do 10 k = 2, n
               s = s + A(k,i)*A(k,j)
  10        continue
            W(i,j) = s
            W(j,i) = s
  20     continue
         c(j) = A(1,j)*b(1)
         do 30 i = 2, n
            c(j) = c(j) + A(i,j)*b(i)
  30  continue
c					factorize [W] = [L][D][L^t]
      det = facld( W, iw, m )
      if( det .eq. zero ) then
         sqnml = -one
         return
      endif
c					solve [W]{c} = {r} for {c}
      call slvld( W, iw, m, c )
c					compute residual sum of squared errors
      sqnml = zero
      do 50 i = 1, n
         s = b(i)
         do 40 k = 1, m
            s = s - c(k)*A(i,k)
  40     continue
         sqnml = sqnml + s*s
  50  continue

      return
      end
