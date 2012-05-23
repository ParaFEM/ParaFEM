c  **********************************************************************
c  *                                                                    *
c  *                           Function dqrgsq                          *
c  *                                                                    *
c  **********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to solve the general least squares problem MIN|| b - Ac || for
c            the vector {c} using the QR decomposition of [A] and return the
c            residual sum of squared errors. (Given's version)
c
c  This routine takes the [Q][R] decomposition of a matrix [A] (see qrgdec)
c  and solves the general least squares problem which is to minimize the sum
c  of squared errors (using the 2-norm)
c
c                       2                           2
c     || {b} - [A]{c} ||   = || [Q^t]{b} - [R]{c} ||
c
c  for the vector of model coefficients {c}. Arguments to the routine
c  are as follows;
c
c    R    real array of size at least m x m which contains the elements of
c         the upper triangular matrix [R] = [Q^t][A] stored in its upper
c         triangle. (input)
c
c   ir    integer giving the leading dimension of R exactly as specified in
c         the dimension statement of the calling routine. (input)
c
c    Q    real array of size at least n x n which contains
c         the elements of the orthogonal transformation matrix required to
c         rotate [A] into [R], ie [Q^t][A] = [R]. (input)
c
c   iq    leading dimension of Q exactly as specified in the calling routine.
c         (input)
c
c    b    real vector of length at least n which contains the RHS. (input)
c
c    c    real vector of length at least n which on output will contain
c         the solution {c} in its first m elements. (output)
c
c    n    column dimension of the matrix [Q], n > m-1. (input)
c
c    m    row dimension of the matrix [R]. (input)
c
c ierr    integer error flag which is set to zero on normal execution, and
c         to -1 if the problem is found to be singular. (output)
c
c  The function returns the residual sum of squared errors. If the system
c  is found to be singular, the function returns -1.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      real*8 function dqrgsq( R, ir, Q, iq, b, c, n, m )
      implicit real*8 (a-h,o-z)
      dimension R(ir,*), Q(iq,*), b(*), c(*)
      data zero/0.d0/, one/1.d0/
c					update RHS with [Q^t]{b}
      do 10 j = 1, n
         c(j) = Q(1,j)*b(1)
         do 10 i = 2, n
            c(j) = c(j) + Q(i,j)*b(i)
  10  continue
c					solve [R}{c} = [Q^t]{b} for {c}
      dqrgsq = -one
      do 20 j = m, 1, -1
         if( R(j,j) .eq. zero ) return
         c(j) = c(j)/R(j,j)
         do 20 i = 1, j-1
            c(i) = c(i) - R(i,j)*c(j)
  20  continue
c					obtain sum of squared errors
      dqrgsq = c(m+1)*c(m+1)
      do 30 k = m+2, n
         dqrgsq = dqrgsq + c(k)*c(k)
  30  continue

      return
      end
