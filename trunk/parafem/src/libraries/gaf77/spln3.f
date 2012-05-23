c  ********************************************************************
c  *                                                                  *
c  *                         subroutine spln3                         *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.11
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to produce a set of cubic splines interpolating n data points.
c
c  This routine takes a set of n data points {x_i, f(x_i)}, i = 1, 2, ..., n,
c  and computes the coefficients of the (n-1) interpolating cubic splines;
c
c    S_j(x) = a_j + b_j*(x - x_j) + c_j*(x - x_j)^2 + d_j*(x - x_j)^3
c
c  for j = 1, 2, ..., n-1. Each S_j(x) is valid on the interval [x_j, x_{j+1}].
c  The splines produced are either natural or clamped, depending on the value
c  of LCLAMP. See also clamp1 and clampn. Arguments to the routine are as
c  follows;
c
c     x    real vector of length at least n which contains the locations
c          of the known function values f(x_i). (input)
c
c     a    real vector of length at least n which contains the observed
c          function values f(x_i), i = 1, 2, ..., n. Note that the vector
c          `a' is also the first coefficient of the set of cubic splines.
c          (input)
c
c     b    real vector of length at least (n-1) which on output will contain
c          the second coefficient of the set of cubic splines. (output)
c
c     c    real vector of length at least (n-1) which on output will contain
c          the third coefficient of the set of cubic splines. (output)
c
c     d    real vector of length at least (n-1) which on output will contain
c          the fourth coefficient of the set of cubic splines. (output)
c
c     n    the number of data points for which the fit is required. (input)
c
c LCLAMP   logical flag which is true if the cubic spline is to be `clamped'
c          at the end points with first derivatives as specified by `clamp1'
c          and `clamp2'. If LCLAMP is false, then a natural spline is used
c          (having zero second derivative at the end points). (input)
c
c clamp1   if LCLAMP is true, then this is a real value giving the derivative
c          of the spline at the first data point. If LCLAMP is false, this
c          is ignored. (input)
c
c clamp2   if LCLAMP is true, then this is a real value giving the derivative
c          of the spline at the last data point. If LCLAMP is false, this
c          is ignored. (input)
c
c   ierr   integer error flag: ierr = 0 if all goes well, ierr = -1 if the
c          system is singular (usually means that two or more of the data
c          points are identical). (output)
c
c  REVISION HISTORY:
c  1.1	eliminated unused variable `one' (Dec 5/96)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine spln3( x, a, b, c, d, n, LCLAMP, clamp1, clamp2, ierr)
      dimension x(*), a(*), b(*), c(*), d(*)
      logical LCLAMP
      data zero/0.0/, half/0.5/, onept5/1.5/, two/2.0/
      data three/3.0/, forpt5/4.5/

      ierr = -1
      if( n .lt. 2 ) return
c						set up h (store in b)
      do 10 i = 1, n-1
         b(i) = x(i+1) - x(i)
         if( b(i) .eq. zero ) return
  10  continue
c						set up tridiagonal system
c							put diagonals in d
c							put RHS in c
      if( LCLAMP ) then
         d(1)   = two*b(2) + onept5*b(1)
         d(n-2) = two*b(n-2) + onept5*b(n-1)
         c(2)   =    (three*(a(3) - a(2))/b(2))
     >            - (forpt5*(a(2) - a(1))/b(1))
     >            + onept5*clamp1
         c(n-1) =  (forpt5*(a(n)   - a(n-1))/b(n-1))
     >            - (three*(a(n-1) - a(n-2))/b(n-2))
     >            - onept5*clamp2
      else
         d(1)   = two*(b(1) + b(2))
         d(n-2) = two*(b(n-2) + b(n-1))
         c(2)   = three*((a(3) - a(2))/b(2) - (a(2) - a(1))/b(1))
         c(n-1) = three*((a(n)   - a(n-1))/b(n-1)
     >                -  (a(n-1) - a(n-2))/b(n-2))
      endif
      do 20 i = 3, n-2
         c(i)   = three*((a(i+1) - a(i))/b(i) - (a(i) - a(i-1))/b(i-1))
         d(i-1) = two*(b(i) + b(i-1))
  20  continue
c						solve the tridiagonal system
      if( d(1) .eq. zero ) return
      do 30 i = 2, n-2
         r = b(i)/d(i-1)
         d(i) = d(i) - r*b(i)
         if( d(i) .eq. zero ) return
         c(i+1) = c(i+1) - r*c(i)
  30  continue
      c(n-1) = c(n-1)/d(n-2)
      do 40 i = n-2, 2, -1
         c(i) = (c(i) - b(i)*c(i+1))/d(i-1)
  40  continue
c						solve for remaining coeff's
      if( LCLAMP ) then
         c(1) = onept5*((a(2) - a(1))/b(1) - clamp1)/b(1) - half*c(2)
         cn   = onept5*(clamp2-(a(n)-a(n-1))/b(n-1))/b(n-1)-half*c(n-1)
      else
         c(1) = zero
         cn   = zero
      endif

      do 50 i = 1, n-2
         d(i) = (c(i+1) - c(i))/(three*b(i))
         b(i) = (a(i+1) - a(i))/b(i) - b(i)*(two*c(i) + c(i+1))/three
  50  continue
      d(n-1) = (cn - c(n-1))/(three*b(n-1))
      b(n-1) = (a(n) - a(n-1))/b(n-1) - b(n-1)*(two*c(n-1) + cn)/three

c						all done
      ierr = 0
      return
      end
