c  ***********************************************************************
c  *                                                                     *
c  *                          Function modnr                             *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Feb. 5, 1991
c
c  PURPOSE  to obtain a root of a user-supplied function via the modified
c	    Newton-Raphson algorithm.
c
c  This function takes as input the name of the user-supplied function, 
c  the name of the user-supplied first derivative function, the name of
c  the user supplied second derivative function, a starting value,
c  and the tolerance, and estimates the closest root via the modified
c  Newton-Raphson algorithm. Note that roots very close to zero will not
c  necessarily be obtained with a relative accuracy less than the specified
c  tolerance but will have an absolute accuracy of at least
c  tol*(machine epsilon). For example, if tol = 1.0e-4, then a true root
c  equal to 2.0e-20 may be represented as zero in single precision.
c
c  This routine checks for the occurrence of a zero slope. If it occurs
c  while G(x) is not equal to zero, then convergence will never occur (the
c  iteration keeps returning to xo). In this case the `convergence failed'
c  flag is set (nit < 0) and the value of xo is returned (the point at
c  which the zero slope occurred).
c
c  Arguments to the routine are as follows
c
c     G   user-supplied external function. It is assumed that calls to
c         this function take the form G(x), where x is the value at which
c         the function value is desired. (input)
c
c    DG   user-supplied external function. DG(x) returns the first
c         derivative of the function G(x) at the point x. (input)
c
c   DDG   user-supplied external function. DDG(x) returns the second
c         derivative of the function G(x) at the point x. (input)
c
c    xs   real starting value from which the search starts. (input)
c
c   tol   real value which gives the maximum acceptable relative error. (input)
c
c maxit   integer giving the maximum number of iterations allowed. (input)
c
c   nit   integer giving the number of iterations actually performed. If
c         nit is less than zero on return, then convergence failed; (output)
c          = -1 convergence failed in maxit iterations
c          = -2 convergence failed due to zero slope
c          = -3 convergence failed due to zero denominator
c
c-----------------------------------------------------------------------------
      real function modnr( G, DG, DDG, xs, tol, maxit, nit )
      external G, DG, DDG
      data zero/0.0/, eps/0.119209E-06/
c					first guess
      xo = xs
      modnr = xo
c					begin iterations
      do 10 nit = 1, maxit
c					check to see if G(x) --> 0 (converged)
         fx = G(xo)
         if( fx .eq. zero ) return
c					obtain the slope at xo and check
         sl = DG(xo)
         if( sl .eq. zero ) go to 20
c					check denominator
         den = sl*sl - fx*DDG(xo)
         if( den .eq. zero ) go to 30
c					new estimate of root
         modnr = xo - fx*sl/den
c					relative error estimate
         if( modnr .eq. zero ) then
            e = abs((modnr - xo)/eps)
         else
            e = abs((modnr - xo)/modnr)
         endif
c					check if error < tolerance
         if( e .lt. tol ) return
c					reset `old' value
         xo = modnr
c					go back for next estimate
  10  continue
c					convergence failed
      nit = -1
      return
  20  nit = -2
      return
  30  nit = -3
      return
      end
