c  **********************************************************************
c  *                                                                    *
c  *                          Function qrhim                            *
c  *                                                                    *
c  **********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to solve the general least squares problem MIN|| b - Ac || for
c            the vector {c} given the QR decomposition of [A] stored in [A]
c            and return the sum of squared errors. This routine performs
c            interative improvement.
c
c  This routine takes an n x m  (n > m) matrix [A] which has been decomposed
c  into an upper triangular matrix [R] and an orthogonal transformation [Q]
c  stored directly in [A] in factorized form (Householder's method),
c  and solves the general least squares problem which is to minimize
c  the sum of squared errors (using the 2-norm)
c
c                       2                           2
c     || {b} - [A]{c} ||   = || [Q^t]{b} - [R]{c} ||
c
c  for the vector of model coefficients {c}. For values of nit > 0, the
c  routine proceeds to improve on the solution by solving the residual
c  problem iteratively. The residual is calculated as
c
c       {r} = {b} - [A]{c}
c
c  whereupon the least squares problem
c
c     || {r} - [A]{z} ||   = || [Q^t]{r} - [R]{z} ||
c
c  is solved for {z}. The improved solution is then {c} = {c} + {z}.
c  Arguments to the routine are as follows;
c
c    A    real array of size at least n x m which contains the elements of
c         the matrix [A]. (input)
c
c   ia    integer giving the leading dimension of A exactly as specified in
c         the dimension statement of the calling routine. (input)
c
c    Q    real array of size at least n x m which contains the elements of
c         the upper triangular matrix [R] stored in its upper triangle and
c         the matrix [Q] stored in factorized form in its lower sub-triangle.
c         See qrhdec. (input)
c
c   iq    integer giving the leading dimension of Q exactly as specified in
c         the dimension statement of the calling routine. (input)
c
c    b    real vector of length at least n which contains the RHS. (input)
c
c    c    real vector of length at least m which on output will contain
c         the solution vector. (output)
c
c    r    temporary real vector of length at least n.
c
c    n    column dimension of the matrix [A]. (input)
c
c    m    row dimension of the matrix [A]. (input)
c
c  nit    the number of improvement interations to perform. If nit = 0, then
c         just the initial solution is returned. (input)
c
c  The function returns the residual sum of squared errors. If the system
c  is found to be singular, the function returns -1.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      function qrhim( A, ia, Q, iq, b, c, r, n, m, nit )
      dimension A(ia,*), Q(iq,*), b(*), r(*), c(*)
      real*8 sum, dble
      data zero/0.0/, one/1.0/, two/2.0/, tol/1.1e-6/
c					check diagonal elements of R
      qrhim = -one
      do 20 i = 1, m
         if( Q(i,i) .eq. zero ) return
  20  continue
c					initialize first solution
      do 10 i = 1, m
         c(i) = zero
  10  continue
c					obtain initial and subsequent solutions
      sqerr = -two
      nn = max0( 0, nit )
      do 100 iter = 0, nn
c						compute residual
         qrhim = zero
         do 60 i = 1, n
            sum = dble(b(i)) - dble(A(i,1))*dble(c(1))
            do 50 j = 2, m
               sum = sum - dble(A(i,j))*dble(c(j))
  50        continue
            r(i) = sum
            qrhim = qrhim + sum*sum
  60     continue
c						any point in continuing?

         if( abs( sqerr - qrhim ) .lt. tol ) return
         sqerr = qrhim
c						update RHS with [Q^t]{r}
         do 80 k = 1, m
            sn = r(k)
            sd = one
            do 70 j = k+1, n
               sn = sn + Q(j,k)*r(j)
               sd = sd + Q(j,k)*Q(j,k)
  70        continue
            s = -two*sn/sd
            r(k) = r(k) + s
            do 80 j = k+1, n
               r(j) = r(j) + s*Q(j,k)
  80     continue
c						solve [R}{c} = [Q^t]{r} for {c}
         do 90 j = m, 1, -1
            z    = r(j)/Q(j,j)
            c(j) = c(j) + z
            do 90 i = 1, j-1
               r(i) = r(i) - Q(i,j)*z
  90     continue
 100  continue
c						compute final error
      qrhim = zero
      do 120 i = 1, n
         sum = dble(b(i)) - dble(A(i,1))*dble(c(1))
         do 110 j = 2, m
            sum = sum - dble(A(i,j))*dble(c(j))
 110     continue
         qrhim = qrhim + sum*sum
 120  continue

      return
      end
