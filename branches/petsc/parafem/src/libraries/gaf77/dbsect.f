c  ***********************************************************************
c  *                                                                     *
c  *                         Function dbsect                             *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Jan. 2, 1991
c
c  PURPOSE  to obtain the roots of a user-supplied function via the Bisection
c	    algorithm.
c
c  This function takes as input the name of the user-supplied function, the
c  bracketing values, and the tolerance, and estimates the enclosed root
c  via the bisection algorithm. Note that roots very close
c  to zero will not necessarily be obtained with a relative accuracy less
c  than the specified tolerance but will have an absolute accuracy of
c  at least tol*(machine epsilon). For example, if tol = 1.0e-4, then
c  this function will accept ~2.e-20 as a root where the true root is
c  zero. Conversely, a true root equal to 2.e-20 may be represented as
c  zero (in double precision, machine epsilon = 2.e-16).
c
c  The bisection algorithm is employed, involving successively halving the
c  interval, to locate the root. Thus if a root does occur between the
c  two starting values, this method is sure to succeed if the sign of f(x)
c  at the two starting values differs. If the sign of f(x) at the two starting
c  values is the same, then convergence to a unique root, or any root at
c  all, is not guaranteed. In fact the method will fail and report that the
c  number of iterations is -1 if the function does not cross the axis in
c  either the left or right half interval after any subdivision.
c
c  Arguments to the routine are as follows
c
c    f    user-supplied external function. It is assumed that calls to
c         this function take the form f(x), where x is the value at which
c         the function value is desired. (input)
c
c  xl, xu real values which bracket the root. (input)
c
c   tol   real value which gives the maximum acceptable relative error. (input)
c
c maxit   integer giving the maximum number of iterations allowed. (input)
c
c   nit   integer giving the number of iterations actually performed. If
c         nit is greater than maxit on return, then convergence did not
c         occur to the desired error tolerance. Generally this can be
c         fixed by increasing either the tolerance or the value of maxit.
c         If nit has value -1 on return, then the function never crossed the
c         axis during the search. In this case the boundary points probably
c         need to be changed. (output)
c-----------------------------------------------------------------------------
      real*8 function dbsect( f, xl, xu, tol, maxit, nit )
      implicit real*8 (a-h,o-z)
      external f
      data zero/0.d0/, half/0.5d0/, deps/2.22045d-16/
      abs(y) = dabs(y)

      x1     = xl
      x2     = xu
      xo     = x1
      g1     = f(x1)
      g2     = f(x2)
c					begin interations
      do 10 nit = 1, maxit
c					new estimate of root
         dbsect = half*(x1 + x2)
c					relative error estimate
         if( abs(dbsect) .lt. deps ) then
            e = abs( (dbsect - xo)/deps )
         else
            e = abs( (dbsect - xo)/dbsect )
         endif
c					check if error < tolerance
         if( e .lt. tol ) return
c					which side is the root on? (minimize
c					the number of calls to f(x))
         g0 = f(dbsect)
         if( g0 .eq. zero ) return
         if( g1*g0 .lt. zero ) then
            x2 = dbsect
            g2 = g0
         elseif( g2*g0 .lt. zero ) then
            x1 = dbsect
            g1 = g0
         else
            go to 20
         endif
         xo = dbsect

  10  continue
c					convergence failed
      return
c				line doesn't cross axis in either interval!
  20  nit = -1
      return
      end
