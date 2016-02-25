c  ***********************************************************************
c  *                                                                     *
c  *                         Function dsecnt                             *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Feb. 5, 1991
c
c  PURPOSE  to obtain a root of a user-supplied function via the Secant
c	    algorithm.
c
c  This function takes as input the name of the user-supplied function, 
c  two starting values (from which to estimate the initial slope),
c  and the tolerance, and estimates the closest root via the Secant
c  algorithm. Note that roots very close to zero will not necessarily be
c  obtained with a relative accuracy less than the specified tolerance but
c  will have an absolute accuracy of at least tol*(machine epsilon). For
c  example, if tol = 1.0e-4, then a true root equal zero to may be
c  represented as 2.0e-20 (in double precision, machine epsilon = 2.e-16).
c
c  The routine checks for the occurrence of a zero slope. If encountered,
c  the `zero-slope encountered' flag is set (nit = -1) and dsecnt returns the
c  last point defining the zero slope. If convergence is not attained, the
c  value of nit returned will be greater than maxit. This can often be solved
c  by increasing maxit and/or tol.
c
c  Arguments to the routine are as follows
c
c     f   user-supplied external function. It is assumed that calls to
c         this function take the form f(x), where x is the value at which
c         the function value is desired. (input)
c
c xa,xb   real starting values from which the search starts. (input)
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
c         be fixed by choosing different starting points. (output)
c-----------------------------------------------------------------------------
      real*8 function dsecnt( f, xa, xb, tol, maxit, nit )
      implicit real*8 (a-h,o-z)
      external f
      data zero/0.d0/, deps/2.22045d-16/
      abs(y) = dabs(y)
c					initialize
      x0     = xa
      x1     = xb
      u0     = f(x0)
      dsecnt = x1
c					begin iterations
      do 10 nit = 1, maxit
c					check to see if G(x) --> 0 (converged)
         u1 = f(x1)
         if( u1 .eq. zero ) return
c					zero slope implies convergence failure
         if( u0 .eq. u1 ) go to 20
c					new estimate of root
         dsecnt = x1 - u1*(x0-x1)/(u0-u1)
c					relative error estimate
         if( abs(dsecnt) .lt. deps ) then
            e = abs((dsecnt - x1)/deps)
         else
            e = abs((dsecnt - x1)/dsecnt)
         endif
c					check if error < tolerance
         if( e .lt. tol ) return
c					save `old' values
         u0 = u1
         x0 = x1
         x1 = dsecnt
c					go back for next estimate
  10  continue
c					convergence failed
      return
c					zero slope encountered
  20  nit = -1
      return
      end
