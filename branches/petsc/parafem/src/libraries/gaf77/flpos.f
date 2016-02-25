c  ***********************************************************************
c  *                                                                     *
c  *                          Function flpos                             *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Jan. 2, 1991
c  Latest Update: Dec 5, 1996
c
c  PURPOSE  to obtain the roots of a user-supplied function via the False
c	    Position algorithm.
c
c  This function takes as input the name of the user-supplied function, the
c  bracketing values, and the tolerance, and estimates the enclosed root
c  via the false-position algorithm. Note that roots very close
c  to zero will not necessarily be obtained with a relative accuracy less
c  than the specified tolerance but will have an absolute accuracy of
c  at least tol*(machine epsilon). For example, if tol = 1.0e-4, then
c  this function will accept ~1.2e-11 as a root where the true root is
c  zero. Conversely, a true root equal to 1.2e-11 may be represented as
c  zero (in single precision, machine epsilon = 1.2e-07).
c
c  The false-position algorithm is employed to locate the root. Thus if a
c  root does occur between the two starting values, this method is sure to
c  succeed if the sign of f(x) at the two starting values differs (and if
c  MAXIT is large enough). If the sign of f(x) at the two starting values
c  is the same, then convergence to a unique root, or any root at all, is
c  not guaranteed. In this case if a root is not found, the false-position
c  algorithm reports that the number of iterations is -1.
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
c         axis during the search and/or a zero slope was encountered.
c         In this case the boundary points probably need to be changed.
c         (output)
c
c  REVISION HISTORY:
c  1.1	eliminated unused variable `half' (Dec 5/96)
c-----------------------------------------------------------------------------
      real function flpos( f, xl, xu, tol, maxit, nit )
      external f
      data zero/0.0/, eps/1.19209e-07/

      x1     = xl
      x2     = xu
      xo     = x1
      g1     = f(x1)
      g2     = f(x2)
c					begin interations
      do 10 nit = 1, maxit
c					zero slope implies convergence failure
         if( g2 .eq. g1 ) go to 20
c					new estimate of root
         flpos = (g2*x1 - g1*x2)/(g2 - g1)
c					relative error estimate
         if( abs(flpos) .lt. eps ) then
            e = abs( (flpos - xo)/eps )
         else
            e = abs( (flpos - xo)/flpos )
         endif
c					return if error < tolerance
         if( e .lt. tol ) return
c					which side is the root on? (minimize
c					the number of calls to f(x))
         g0 = f(flpos)
         if( g0 .eq. zero ) return
         if( g1*g0 .lt. zero ) then
            x2 = flpos
            g2 = g0
         elseif( g2*g0 .lt. zero ) then
            x1 = flpos
            g1 = g0
         else
            go to 20
         endif
         xo = flpos

  10  continue
c					convergence failed
      return
c					no axis crossing within interval
  20  nit = -1
      return
      end
