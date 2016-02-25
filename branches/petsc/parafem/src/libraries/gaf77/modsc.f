c  ***********************************************************************
c  *                                                                     *
c  *                          Function modsc                             *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Feb. 5, 1991
c
c  PURPOSE  to obtain a root of a user-supplied function via the modified
c	    Secant algorithm.
c
c  This function takes as input the name of the user-supplied function,
c  the name of the user-supplied first derivative function, two starting
c  values, and the tolerance, and estimates the closest root via the
c  modified Secant algorithm. Note that roots very close to zero will not
c  necessarily be obtained with a relative accuracy less than the specified
c  tolerance but will have an absolute accuracy of at least
c  tol*(machine epsilon). For example, if tol = 1.0e-4, then a true root
c  equal to 2.0e-20 may be represented as zero in single precision.
c
c  The routine checks for the occurrence of a zero slope. If encountered,
c  the `convergence failed' flag is set (nit < 0) and modsc returns the
c  point at which the zero slope occurred.
c
c  Arguments to the routine are as follows
c
c     G   user-supplied external function. It is assumed that calls to
c         this function take the form G(x), where x is the value at which
c         the function value is desired. (input)
c
c    DG   user-supplied external function. DG(x) returns the first
c         derivative of G(x) at the point x. (input)
c
c xa,xb   real starting values from which the search starts. (input)
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
      real function modsc( G, DG, xa, xb, tol, maxit, nit )
      external G, DG
      data zero/0.0/, eps/0.119209E-06/
c					initialize
      x0     = xa
      x1     = xb
      modsc  = x0
      nit    = 0
c					check DG(x0) = 0?
      sl = DG(x0)
      if( sl .eq. zero ) go to 20
      fx = G(x0)
      if( fx .eq. zero ) return
      u0    = fx/sl
      modsc = x1
c					begin iterations
      do 10 nit = 1, maxit
c					check to see if G(x) --> 0 (converged)
         fx = G(x1)
         if( fx .eq. zero ) return
c					zero slope implies convergence failure
         sl = DG(x1)
         if( sl .eq. zero ) go to 20
         u1 = fx/sl
c					zero denominator also implies failure
         den = u0 - u1
         if( den .eq. zero ) go to 30
c					new estimate of root
         modsc = x1 - u1*(x0-x1)/den
c					relative error estimate
         if( modsc .eq. zero ) then
            e = abs( (modsc - x1)/eps )
         else
            e = abs( (modsc - x1)/modsc )
         endif
c					check if error < tolerance
         if( e .lt. tol ) return
c					save `old' values
         u0 = u1
         x0 = x1
         x1 = modsc
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
