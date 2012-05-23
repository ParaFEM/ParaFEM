c  ***********************************************************************
c  *                                                                     *
c  *                           Function dpnorm                           *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  return the p-norm of a vector
c
c  This routine accepts as input a vector of length n and computes the
c  associated p-norm. If
c
c     a) p = 1, then norm = |x(1)| + |x(2)| + ... + |x(n)|
c
c     b) p = 2, then norm = sqrt[ x(1)*x(1) + x(2)*x(2) + ... + x(n)*x(n) ]
c
c     c) p = other, then norm = max_i |x(i)|    (infinite norm)
c
c  all of which are `equivalent' norms associated with the `size' of the
c  vector {x}. Arguments to the function are;
c
c    X    real vector of length at least n. (input)
c
c    n    integer length of the vector X. (input)
c
c    p    integer flag denoting the type of norm to take (see above). (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c------------------------------------------------------------------------------
      real*8 function dpnorm( X, n, p )
      implicit real*8 (a-h,o-z)
      dimension X(*)
      integer p
      abs(xx)  = dabs(xx)
      sqrt(xx) = dsqrt(xx)
c					switch to norm type p = 1, 2, infinite
c						1-norm
      if( p .eq. 1 ) then
         dpnorm = abs( x(1) )
         do 10 i = 2, n
            dpnorm = dpnorm + abs( x(i) )
  10     continue
c						2-norm
      else if( p .eq. 2 ) then
         dpnorm = x(1)*x(1)
         do 20 i = 2, n
            dpnorm = dpnorm + x(i)*x(i)
  20     continue
         dpnorm = sqrt(dpnorm)
c						infinite-norm
      else
         dpnorm = abs( x(1) )
         do 30 i = 2, n
            s = abs( x(i) )
            if( s .gt. dpnorm ) dpnorm = s
  30     continue
      endif

      return
      end
