c  ***********************************************************************
c  *                                                                     *
c  *                          Function fxpnt                             *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Jan. 2, 1991
c  Latest Update: Dec 5, 1996
c
c  PURPOSE  to obtain the roots of a user-supplied function via the Fixed-
c	    Point Iteration algorithm
c
c  This function takes as input the name of the user-supplied function, the
c  starting value, and the tolerance, and estimates the closest root
c  via the Fixed-Point Iteration algorithm. Note that roots very close
c  to zero will not necessarily be obtained with a relative accuracy less
c  than the specified tolerance but will have an absolute accuracy of
c  at least tol*(machine epsilon). For example, if tol = 1.0e-4, then
c  a true root equal to 2.0e-20 may be represented as zero in single
c  precision. Arguments to the routine are as follows
c
c    G    user-supplied external function. It is assumed that calls to
c         this function take the form G(x), where x is the value at which
c         the function value is desired. (input)
c
c   xs    real value giving the point at which the search should start. (input)
c
c   tol   real value which gives the maximum acceptable relative error. (input)
c
c maxit   integer giving the maximum number of iterations allowed. (input)
c
c   nit   integer giving the number of iterations actually performed. If
c         nit is less than 0 on return, then convergence failed. (output)
c
c  REVISION HISTORY:
c  1.1	eliminated unused variable `half' (Dec 5/96)
c-----------------------------------------------------------------------------
      real function fxpnt( G, xs, tol, maxit, nit )
      external G
      data zero/0.0/, eps/0.119209E-06/

      xo = xs
c					begin interations
      do 10 nit = 1, maxit
c					new estimate of root
         fxpnt = G(xo)
c					relative error estimate
         if( fxpnt .eq. zero ) then
            e = abs( (fxpnt - xo)/eps )
         else
            e = abs( (fxpnt - xo)/fxpnt )
         endif
c					check if error < tolerance
         if( e .lt. tol ) return
c					set old `x'
         xo = fxpnt

  10  continue
c					convergence failed`
      nit = -1
      return
      end
