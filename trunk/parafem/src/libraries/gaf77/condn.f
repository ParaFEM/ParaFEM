c  *********************************************************************
c  *                                                                   *
c  *                        subroutine condn                           *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, Apr. 15, 1993
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  computes the 1-norm and infinite norm condition numbers of the
c           supplied matrix
c
c  This routine computes the condition number of a supplied matrix using both
c  the 1-norm and the infinite-norms of the matrix and its inverse. As a
c  byproduct of the process, the inverse of A is also returned. Iterative
c  improvement is employed to maximize the accuracy of the result.
c  Arguments to the routine are as follows;
c
c     A      real array of size at least N x N which contains the matrix for
c            which the condition numbers are required. (input)
c
c     AINV   real array of size at least N x N which on output will contain
c            the inverse matrix of the array A. (output)
c
c     ATMP   real array of size at least N x N which is used for workspace.
c
c     ia     leading dimension of the arrays A and AINV (and ATMP) as
c            specified in the calling routine. (input)
c
c     n      the size of the matrix A. (input)
c
c     BTMP   real vector of length at least 3*N which is used for workspace.
c
c     ITMP   integer vector of length at least N which is used for workspace.
c
c     C1     real value which on output contains the 1-norm condition number
c            of the matrix A. If A is found to be algorithmically singular,
c            C1 is set to the largest floating point number possible. (output)
c
c     CINF   real value which on output contains the infinite-norm condition
c            number of the matrix A. If A is found to be algorithmically
c            singular, CINF is set to the largest floating point number
c            possible (see "big"). (output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine condn( A, AINV, ATMP, ia, n, BTMP, ITMP, C1, CINF )
      dimension A(ia,*), AINV(ia,*), ATMP(ia,*), BTMP(*)
      integer ITMP(*)
      data zero/0.0/, one/1.0/, big/1.7014117e+38/
c					obtain 1 and infinite norms of A
      anorm = zero
      cnorm = zero
      do 20 i = 1, n
         btmp(i) = zero
         s = abs(A(i,1))
         t = abs(A(1,i))
         ATMP(i,1) = A(i,1)
         do 10 j = 2, n
            s = s + abs(A(i,j))
            t = t + abs(A(j,i))
            ATMP(i,j) = A(i,j)
  10     continue
         if( s .gt. anorm ) anorm = s
         if( t .gt. cnorm ) cnorm = t
  20  continue
c					LU decompose [A] ==> ATMP

      call spplu( ATMP, ia, n, itmp, btmp(n+1), det )
      if( det .eq. zero ) then
         C1   = big
         CINF = big
         return
      endif
c					now find inverse matrix
      do 30 j = 1, n
         btmp(j) = one
         call psolv( ATMP, ia, n, itmp, btmp, AINV(1,j) )
         call imprv( A, ATMP, ia, n, itmp, AINV(1,j), btmp, btmp(n+1),
     >               btmp(2*n+1), 1 )
         btmp(j) = zero
  30  continue
c					obtain infinite norm of A-inv
      bnorm = zero
      dnorm = zero
      do 50 i = 1, n
         s = abs(AINV(i,1))
         t = abs(AINV(1,i))
         do 40 j = 2, n
            s = s + abs(AINV(i,j))
            t = t + abs(AINV(j,i))
  40     continue
         if( s .gt. bnorm ) bnorm = s
         if( t .gt. dnorm ) dnorm = t
  50  continue

      C1   = cnorm*dnorm
      CINF = anorm*bnorm

      return
      end
