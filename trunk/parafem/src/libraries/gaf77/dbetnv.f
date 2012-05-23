c  ********************************************************************
c  *                                                                  *
c  *                         function dbetnv                          *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 3.0
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Dec 1, 1998
c
c  PURPOSE  returns the inverse Beta distribution function for given
c           probability p and parameters alpha and beta.
c
c  This routine employs the Newton-Raphson root finding technique to
c  obtain the inverse Beta distribution. The Beta distribution is defined
c  by the incomplete Beta function
c
c                G(a+b)     x
c     F(x) =   ---------  INT (1-t)**(b-1) * t**(a-1) dt
c              G(a)*G(b)    0
c
c  where G(.) is the Gamma function (a! if `a' is an integer).
c  This routine returns x given p = F(x). Relative error tolerance is
c  fixed at 1.e-10.
c  Note that the modified Newton-Raphson algorithm, designed for multiple
c  roots, does not work very well on this function.
c  Arguments to the routine are as follows;
c
c     p   real value giving the probability at which the corresponding x value
c         is desired. (input)
c
c     a   the alpha parameter of the Beta distribution. (input)
c
c     b   the beta parameter of the Beta distribution. (input)
c
c  ierr   integer error flag which is 0 if all goes well, -1 if convergence
c         not achieved.
c
c  REVISION HISTORY:
c  2.0	added some checks of the input parameters (May 10/94)
c  3.0	tossed the innacurate Abramowitz and Stegun approximation in favour
c	of a Newton-Raphson approach (Dec 1/98)
c---------------------------------------------------------------------------
      real*8 function dbetnv( p, a, b, ierr )
      parameter (MAXIT = 40)
      implicit real*8 (a-h,o-z)
      data tol/1.d-10/
      data zero/0.d0/, half/0.5d0/, one/1.d0/
      data pt001/0.001d0/
      data pt999/0.999d0/, ptz99/0.00002475d0/
c					function statements
      beta(z1,z2,z3) = dbeta(z1,z2,z3)
      gamln(z) = dgamln(z)
      float(i) = dble(i)
      exp(z)   = dexp(z)
      abs(z)   = dabs(z)
c-------------------------------- start executable statements ------------
      ierr = 0
c					eliminate the obvious
      if( p .le. zero ) then
         dbetnv = zero
         return
      elseif( p .ge. one ) then
         dbetnv = one
         return
      endif
c					try Newton-Raphson first
c						start at midpoint
      dbetnv = half
c						now refine using Newton-Raphson
      bg = exp( gamln(a+b) - gamln(a) - gamln(b) )
      a0 = a - one
      b0 = b - one
      xo = dbetnv
      do 20 i = 1, MAXIT
         od  = one - dbetnv
         da0 = dbetnv**a0
         db0 = od**b0
c						function and derivatives
         g   = beta(dbetnv,a,b) - p
         gp  = bg*da0*db0
c						Newton-Raphson update
         dbetnv = dbetnv - g/gp
c						ensure we stay in bounds
         if( dbetnv .ge. one ) then
            dbetnv = pt999 + ptz99*float(i)
         elseif( dbetnv .le. zero ) then
            dbetnv = pt001 - ptz99*float(i)
         endif

         if( abs((dbetnv - xo)/dbetnv) .lt. tol ) return
         xo = dbetnv

  20  continue
c					NR Convergence failed, try Bisection
      xo = zero
      x1 = zero
      x2 = one
      g1 = -p
      g2 = one - p
      do 30 j = 1, 2*MAXIT
         dbetnv = half*(x1 + x2)
         e = abs( (dbetnv - xo)/dbetnv )
         if( e .lt. tol ) return
         g0 = beta(dbetnv,a,b) - p
         if( g0 .eq. zero ) return
         if( g1*g0 .lt. zero ) then
            x2 = dbetnv
            g2 = g0
         elseif( g2*g0 .lt. zero ) then
            x1 = dbetnv
            g1 = g0
         else
            go to 40
         endif
         xo = dbetnv
  30  continue
c					give up, convergence failed
  40  ierr = -1

      return
      end
