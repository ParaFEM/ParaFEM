c  ********************************************************************
c  *                                                                  *
c  *                         function dspint                          *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 1.11
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  returns the definite integral of a function through a cubic spline
c           approximation.
c
c  This routine takes a set of n data points {x_i, f(x_i)}, i = 1, 2, ..., n,
c  and computes the coefficients of the (n-1) interpolating cubic splines;
c
c   f(x) ~ S_j(x) = a_j + b_j*(x - x_j) + c_j*(x - x_j)^2 + d_j*(x - x_j)^3
c
c  for j = 1, 2, ..., n-1. Each S_j(x) is valid on the interval [x_j, x_{j+1}]
c  and S_j(x) is an approximation to f(x). The definite integral
c
c                r
c          I = int f(x) dx
c                q
c
c  is then calculated as the sum
c
c              n-1
c          I ~ sum a_i*h_i + (b_i*h_i**2)/2 + (c_i*h_i**3)/3 + (d_i*h_i**4)/4
c              i=1
c
c  in which h_i = (x_{i+1} - x_i)
c  The splines produced are either natural or clamped, depending on the value
c  of LCLAMP. See also clamp1 and clampn. Arguments to the routine are as
c  follows;
c
c     x    real vector of length at least n which contains the locations
c          of the known function values f(x_i). (input)
c
c     f    real vector of length at least n which contains the observed
c          function values f(x_i), i = 1, 2, ..., n. (input)
c
c     b    temporary real vector of length at least 4*(n-1) used for workspace.
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
c  1.1	eliminated unused local variable `one' (Dec 5/96)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      real*8 function dspint( x, f, b, n, LCLAMP, clamp1, clamp2, ierr)
      implicit real*8 (a-h,o-z)
      dimension x(*), f(*), b(*)
      logical LCLAMP
      data zero/0.d0/, half/0.5d0/, onept5/1.5d0/, two/2.d0/
      data three/3.d0/, four/4.d0/, forpt5/4.5d0/

      ierr = -1
      dspint = zero
      if( n .lt. 2 ) return
c						set up h (store in b)
      do 10 i = 1, n-1
         b(3*n-3+i) = x(i+1) - x(i)
         if( b(3*n-3+i) .eq. zero ) return
  10  continue
c						set up tridiagonal system
c							put diagonals in d
c							put RHS in c
      if( LCLAMP ) then
         b(2*n-1) = two*b(3*n-1) + onept5*b(3*n-2)
         b(3*n-4) = two*b(4*n-5) + onept5*b(4*n-4)
         b(n+1)   = (three* (f(3) - f(2))/b(3*n-1))
     >            - (forpt5*(f(2) - f(1))/b(3*n-2))
     >            +  onept5*clamp1
         b(2*n-2) = (forpt5*(f(n)  - f(n-1))/b(4*n-4))
     >            - (three*(f(n-1) - f(n-2))/b(4*n-5))
     >            -  onept5*clamp2
      else
         b(2*n-1) = two*(b(3*n-2) + b(3*n-1))
         b(3*n-4) = two*(b(4*n-5) + b(4*n-4))
         b(n+1)   = three*((f(3)-f(2))/b(3*n-1) - (f(2)-f(1))/b(3*n-2))
         b(2*n-2) = three*((f(n)-f(n-1))/b(4*n-4)
     >                -  (f(n-1)-f(n-2))/b(4*n-5))
      endif
      do 20 i = 3, n-2
         b(n-1+i)   = three*((f(i+1)-f(i))/b(3*n-3+i)
     >                     - (f(i)-f(i-1))/b(3*n-4+i))
         b(2*n-3+i) = two*(b(3*n-3+i) + b(3*n-4+i))
  20  continue
c						solve the tridiagonal system
      if( b(2*n-1) .eq. zero ) return
      do 30 i = 2, n-2
         r = b(3*n-3+i)/b(2*n-3+i)
         b(2*n-2+i) = b(2*n-2+i) - r*b(3*n-3+i)
         if( b(2*n-2+i) .eq. zero ) return
         b(n+i) = b(n+i) - r*b(n-1+i)
  30  continue
      b(n+n-2) = b(n+n-2)/b(3*n-4)
      do 40 i = n-2, 2, -1
         b(n-1+i) = (b(n-1+i) - b(3*n-3+i)*b(n+i))/b(2*n-3+i)
  40  continue
      do 50 i = 2, n-2
         b(n-1+i) = b(n-1+i)/three
  50  continue
c						solve for remaining coeff's
      if( LCLAMP ) then
         b(n) = half*((f(2)-f(1))/b(3*n-2)-clamp1)/b(3*n-2)-half*b(n+1)
         cn   = half*(clamp2-(f(n)-f(n-1))/b(4*n-4))/b(4*n-4)
     >         -half*b(n+n-2)
      else
         b(n) = zero
         cn   = zero
      endif

      do 60 i = 1, n-2
         b(2*n-2+i) = (b(n+i) - b(n-1+i))/(four*b(3*n-3+i))
         b(i) = half*((f(i+1) - f(i))/b(3*n-3+i)
     >        - b(3*n-3+i)*(two*b(n-1+i) + b(n+i)))
  60  continue
      b(3*n-3) = (cn - b(n+n-2))/(four*b(4*n-4))
      b(n-1) = half*((f(n) - f(n-1))/b(4*n-4)
     >       - b(4*n-4)*(two*b(n+n-2) + cn))

c						compute integral

      dspint = b(3*n-2)*(f(1) + b(3*n-2)*(b(1) + b(3*n-2)*(b(n)
     >       + b(3*n-2)*b(2*n-1))))
      do 70 i = 2, n-1
         dspint = dspint + b(3*n-3+i)*(f(i)+b(3*n-3+i)*(b(i)
     >          + b(3*n-3+i)*(b(n-1+i)+b(3*n-3+i)*b(2*n-2+i))))
  70  continue

      ierr = 0
      return
      end
