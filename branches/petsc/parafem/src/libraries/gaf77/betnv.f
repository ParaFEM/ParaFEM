c  ********************************************************************
c  *                                                                  *
c  *                          function betnv                          *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 3.0
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
      real*4 function betnv( p, a, b, ierr )
      parameter (MAXIT = 100)
      data tol/1.e-06/
      data zero/0.0/, half/0.5/, one/1.0/
      data pt001/0.001/
      data pt999/0.999/, ptz99/0.0000099/
c-------------------------------- start executable statements ------------
      ierr = 0
c					eliminate the obvious
      if( p .le. zero ) then
         betnv = zero
         return
      elseif( p .ge. one ) then
         betnv = one
         return
      endif
c					start at midpoint
      betnv = half
c					now refine using Newton-Raphson

      bg = exp( gamln(a+b) - gamln(a) - gamln(b) )
      a0 = a - one
      b0 = b - one
      xo = betnv
      do 20 i = 1, MAXIT
         od  = one - betnv
         da0 = betnv**a0
         db0 = od**b0
c					function and derivatives
         g   = beta(betnv,a,b) - p
         gp  = bg*da0*db0

c					Newton-Raphson update
         betnv = betnv - g/gp
c					ensure we stay in bounds
         if( betnv .ge. one ) then
            betnv = pt999 + ptz99*float(i)
         elseif( betnv .le. zero ) then
            betnv = pt001 - ptz99*float(i)
         endif

         if( abs((betnv - xo)/betnv) .lt. tol ) return
         xo = betnv

  20  continue
c					NR Convergence failed, try Bisection
      xo = zero
      x1 = zero
      x2 = one
      g1 = -p
      g2 = one - p
      do 30 j = 1, 2*MAXIT
         betnv = half*(x1 + x2)
         e = abs( (betnv - xo)/betnv )
         if( e .lt. tol ) return
         g0 = beta(betnv,a,b) - p
         if( g0 .eq. zero ) return
         if( g1*g0 .lt. zero ) then
            x2 = betnv
            g2 = g0
         elseif( g2*g0 .lt. zero ) then
            x1 = betnv
            g1 = g0
         else
            go to 40
         endif
         xo = betnv
  30  continue
c					convergence not achieved
  40  ierr = -1

      return
      end
