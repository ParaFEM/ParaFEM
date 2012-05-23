c  **********************************************************************
c  *                                                                    *
c  *                           Function qrhsq                           *
c  *                                                                    *
c  **********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to solve the general least squares problem MIN|| b - Ac || for
c            the vector {c} given the QR decomposition of [A] stored in [A]
c            and return the sum of squared errors. (Householder version)
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
c  for the vector of model coefficients {c}. Arguments to the routine
c  are as follows;
c
c    A    real array of size at least n x m which contains the elements of
c         the upper triangular matrix [R] stored in its upper triangle and
c         the matrix [Q] stored in factorized form in its lower sub-triangle.
c         See qrhdec. (input)
c
c   ia    integer giving the leading dimension of A exactly as specified in
c         the dimension statement of the calling routine. (input)
c
c    b    real vector of length at least n which on input contains the
c         RHS. On output, b will contain the solution {c} in its first m
c         elements. (input/output)
c
c     n    column dimension of the matrix [A]. (input)
c
c     m    row dimension of the matrix [A]. (input)
c
c  ierr    integer error flag which is set to zero on normal execution, and
c          to -1 if the problem is found to be singular. (output)
c
c  The function returns the residual sum of squared errors. If the system
c  is found to be singular, the function returns -1.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      function qrhsq( A, ia, b, n, m )
      dimension A(ia,*), b(*)
      data zero/0.0/, one/1.0/, two/2.0/
c					update RHS with [Q^t]{b}
      do 20 k = 1, m
         sn = b(k)
         sd = one
         do 10 j = k+1, n
            sn = sn + A(j,k)*b(j)
            sd = sd + A(j,k)*A(j,k)
  10     continue
         s = -two*sn/sd
         b(k) = b(k) + s
         do 20 j = k+1, n
            b(j) = b(j) + s*A(j,k)
  20  continue
c					solve [R}{c} = [Q^t]{b} for {c} --->{b}
      qrhsq = -one
      do 30 j = m, 1, -1
         if( A(j,j) .eq. zero ) return
         b(j) = b(j)/A(j,j)
         do 30 i = 1, j-1
            b(i) = b(i) - A(i,j)*b(j)
  30  continue
c					obtain sum of squared errors
      qrhsq = b(m+1)*b(m+1)
      do 40 k = m+2, n
         qrhsq = qrhsq + b(k)*b(k)
  40  continue

      return
      end
