c  ***********************************************************************
c  *                                                                     *
c  *                          Function newtr                             *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Feb. 5, 1991
c
c  PURPOSE  to obtain a root of a user-supplied function via the Newton-
c	    Raphson algorithm.
c
c  This function takes as input the name of the user-supplied function, 
c  the name of the user-supplied derivative function, a starting value,
c  and the tolerance, and estimates the closest root via the Newton-Raphson
c  algorithm. Note that roots very close to zero will not necessarily be
c  obtained with a relative accuracy less than the specified tolerance but
c  will have an absolute accuracy of at least tol*(machine epsilon). For
c  example, if tol = 1.0e-4, then a true root equal to zero may be
c  represented as 1.2e-11 (in single precision, machine epsilon = 1.2e-07).
c
c  The routine checks for the occurrence of a zero slope. If encountered,
c  the `zero-slope encountered' flag is set (nit = -1) and newtr returns the
c  point at which the zero slope occurred. If convergence is not attained, the
c  value of nit returned will be greater than maxit. This can often be solved
c  by increasing maxit and/or tol.
c
c  Arguments to the routine are as follows
c
c     f   user-supplied external function. It is assumed that calls to
c         this function take the form f(x), where x is the value at which
c         the function value is desired. (input)
c
c    df   user-supplied external function. df(x) returns the derivative of
c         the function f(x) at the point x. (input)
c
c    xs   real starting value from which the search starts. (input)
c
c   tol   real value which gives the maximum acceptable relative error. (input)
c
c maxit   integer giving the maximum number of iterations allowed. (input)
c
c   nit   integer giving the number of iterations actually performed. If
c         nit is greater than maxit on return, then convergence did not
c         occur to the desired error tolerance. Generally this can be
c         fixed by increasing either the tolerance or the value of maxit.
c         If nit has value -1 on return, then a zero slope was encountered
c         during the search and convergence was impossible. This can sometimes
c         be fixed by choosing a different starting point. (output)
c
c  NOTE: be sure to declare newtr as a real function in the calling routine.
c-----------------------------------------------------------------------------
      real*8 function newtr( f, df, xs, tol, maxit, nit )
      external f, df
      data zero/0.0/, eps/1.19209e-07/
c					first guess
      xo     = xs
      newtr = xo
c					begin iterations
      do 10 nit = 1, maxit
c					check to see if G(x) --> 0 (converged)
         fx = f(xo)
         if( fx .eq. zero ) return
c					zero slope implies convergence failure
         sl = df(xo)
         if( sl .eq. zero ) go to 20
c					new estimate of root
         newtr = xo - fx/sl
c					relative error estimate
         if( abs(newtr) .lt. eps ) then
            e = abs((newtr - xo)/eps)
         else
            e = abs((newtr - xo)/newtr)
         endif
c					check if error < tolerance
         if( e .lt. tol ) return
c					reset `old' value
         xo = newtr
c					go back for next estimate
  10  continue
c					convergence failed within maxit
      return
c					zero slope encountered
  20  nit = -1
      return
      end
