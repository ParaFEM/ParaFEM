c  ********************************************************************
c  *                                                                  *
c  *                        Subroutine deigen                         *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 2.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to compute the eigenvalues (and eigenvectors) of a real
c            symmetric matrix
c
c  This routine takes an n x n real symmetric matrix [A] and estimates
c  all of its eigenvalues and, optionally, its eigenvectors. The procedure
c  involves first tridiagonalizing [A] via Householder transformations,
c  then estimating the eigenvalues by iteratively performing Given's
c  rotations (with implicit shifts) on the resulting tridiagonal matrix.
c  The algorithm follows that proposed in "Numerical Recipes in C" by
c  Press etal, pages 367-381. Arguments to the routine are as follows;
c
c      A     real symmetric matrix of size at least n x n which contains
c            the elements of the matrix for which the eigenvalues are
c            desired. The matrix A is overwritten by the associated
c            eigenvectors [Q] on output (see LQ) such that [Q^t][A][Q] = [D]
c            where [D] is a diagonal matrix consisting of the eigenvectors.
c            (input/output)
c
c     ia     leading dimension of A exactly as specified in the calling
c            routine. (input)
c
c      d     real vector of length at least n which on output will contain
c            the set of eigenvalues for the problem. These eigenvalues will
c            not necessarily be ordered by magnitude. (output)
c
c      e     temporary vector of length at least n which is used to store
c            the sub-diagonal terms.
c
c      n     size of the matrix A. (input)
c
c     LQ     logical flag which is true if the eigenvectors are desired. Note
c            that if LQ is false, then A will contain no useful information
c            on return. The computation of the eigenvectors takes significantly
c            longer than the eigenvalues alone (O(n^3) vs. O(n)). (input)
c
c  ORDER     logical flag which is true if the eigenvalues (and corresponding
c            eigenvectors) are to be ordered in descending magnitude. (input)
c
c   ierr     integer flag denoting the status of the run. ierr > 0 if
c            the eigenvalues are successfully computed and -1 if convergence
c            is not achieved within 30 iterations for a particular eigenvalue.
c            If successful, ierr reports the total number of iterations taken
c            in finding all n eigenvalues. (output)
c
c NOTE: the parameter MAXIT controls the maximum number of iterations allowed
c       in the search for each eigenvalue.
c
c  REVISION HISTORY:
c  2.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine deigen( A, ia, d, e, n, LQ, ORDER, ierr )
      implicit real*8 (a-h,o-z)
      parameter (MAXIT = 30)
      dimension A(ia,*), d(*), e(*)
      logical LQ, ORDER
      data zero/0.d0/, one/1.d0/, two/2.d0/
c					real*8 conversion stuff
      abs(x)      = dabs(x)
      sqrt(x)     = dsqrt(x)
      sign(x1,x2) = dsign(x1,x2)

      ierr = 0
c					reduce [A] to tridiagonal form
c					using Householder transformations
      do 70 i = n, 3, -1
         l = i - 1
         s = abs( A(i,1) )
         do 10 k = 2, i-2
            s = s + abs( A(i,k) )
  10     continue
c					do we need to transform this column?
         if( s .eq. zero ) then
            e(i) = A(i,l)
            d(i) = zero
         else
c						yes, perform transformation ...
            s = s + abs( A(i,l ) )
c						compute h = 0.5*V^t*V
            h = zero
            do 20 k = 1, l
               A(i,k) = A(i,k)/s
               h      = h + A(i,k)*A(i,k)
  20        continue
            g      = -sign( sqrt(h), A(i,l) )
            e(i)   = s*g
            h      = h - g*A(i,l)
            A(i,l) = A(i,l) - g
            d(i)   = h
c						store V/h in i'th col. of A
            f      = zero
            do 50 j = 1, l
               A(j,i) = A(i,j)/h
c						compute A*V/h
               g      = zero
               do 30 k = 1, j
                  g = g + A(j,k)*A(i,k)
  30           continue
               do 40 k = j+1, l
                  g = g + A(k,j)*A(i,k)
  40           continue
c						put A*V/h in unused e(j)
               e(j) = g/h
               f    = f + e(j)*A(i,j)
  50        continue
c						hh = V^t*A*V/(V^t*V)
            hh = f/(h+h)
            do 60 j = 1, l
               f    = A(i,j)
               e(j) = e(j) - hh*f
               g    = e(j)
               do 60 k = 1, j
                  A(j,k) = A(j,k) - f*e(k) - g*A(i,k)
  60        continue
         endif
  70  continue
c					set remaining values
      d(1) = zero
      d(2) = zero
      e(1) = zero
      e(2) = A(2,1)
c					accumulate transformation [Q]
      if( LQ ) then
         do 100 i = 1, n
            l = i - 1
            if( d(i) .ne. zero ) then
               do 90 j = 1, l
                  g = zero
                  do 80 k = 1, l
                     g = g + A(i,k)*A(k,j)
  80              continue
                  do 90 k = 1, l
                     A(k,j) = A(k,j) - g*A(k,i)
  90           continue
            endif
            d(i)   = A(i,i)
            A(i,i) = one
            do 100 j = 1, l
               A(j,i) = zero
               A(i,j) = zero
 100     continue
      else
         do 110 i = 1, n
            d(i) = A(i,i)
 110     continue
      endif
c					now iterate for eigenvalues
c						renumber e(i) = e(1), e(2),...
      do 120 i = 2, n
         e(i-1) = e(i)
 120  continue
      e(n) = zero
c					look for l'th eigenvalue
      do 180 l = 1, n
c						maximum MAXIT iterations
         do 170 iter = 1, MAXIT
            ierr = ierr + 1
c						look for zero sub-diagonal
c						to split the matrix
            do 130 m = l, n-1
               dd = abs(d(m)) + abs(d(m+1))
               if( (dd + abs(e(m))) .eq. dd ) go to 140
 130        continue
 140        if( m .ne. l ) then
c						form shift
               g = (d(l+1) - d(l))/(two*e(l))
               r = sqrt( one + g*g )
               g = d(m) - d(l) + e(l)/(g + sign(r,g))

c						apply Given rotations and chase
c						zeroes to restore tridiag. form
               s = one
               c = one
               p = zero
               do 160 i = m-1, l, -1
                  f = s*e(i)
                  b = c*e(i)
                  if( abs( f ) .gt. abs( g ) ) then
                     c = g/f
                     r = sqrt( one + c*c )
                     e(i+1) = f*r
                     s = one/r
                     c = c*s
                  else
                     s = f/g
                     r = sqrt( one + s*s )
                     e(i+1) = g*r
                     c = one/r
                     s = s*c
                  endif
                  g = d(i+1) - p
                  r = (d(i) - g)*s + two*c*b
                  p = s*r
                  d(i+1) = g + p
                  g = c*r - b
c						update eigenvectors
                  if( LQ ) then
                     do 150 k = 1, n
                        f = A(k,i+1)
                        A(k,i+1) = s*A(k,i) + c*f
                        A(k,i) = c*A(k,i) - s*f
 150                 continue
                  endif
 160           continue
               d(l) = d(l) - p
               e(l) = g
               e(m) = zero
            else
c						go back and find next eigenval
               go to 180
            endif
 170     continue
c						if we get here, we've failed
         ierr = -1
         return
 180  continue
c					all done, now order if desired
      if ( ORDER ) then
         do 210 ii = 2, n
            i = ii - 1
            k = i
            p = d(i)
c					find largest eigenvalue
            do 190 j = ii, n
               if( d(j) .gt. p ) then
                  k = j
                  p = d(j)
               endif
 190        continue
c					swap
            if( k .ne. i ) then
               d(k) = d(i)
               d(i) = p
               if( LQ ) then
                  do 200 j = 1, n
                     p      = A(j,i)
                     A(j,i) = A(j,k)
                     A(j,k) = p
 200              continue
               endif
            endif
 210     continue
      endif

      return
      end
