c  *******************************************************************
c  *                                                                 *
c  *                     Function dlsfr3                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.31
c  Written by Gordon A. Fenton, TUNS, Mar. 15, 1995
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points in a separable 3-D
c           fractional Gaussian noise field
c
c  Returns the covariance between two points in the field separated by the
c  lag vector {X,Y,Z}. The covariance function is a separable fractional
c  Gaussian noise, given by
c
c         B(X,Y,Z) = var*r(X)*r(Y)*r(Z)
c
c  for two points separated by distance {X,Y,Z}, where r(X), r(Y), and r(Z)
c  are the correlation functions given by
c
c               1               2H       2H           2H
c     r(X) = -------- [ |X + pb|  -  2|X|  +  |X - pb|  ]		(1)
c                 2H
c             2*pb
c
c
c               1               2G       2G           2G
c     r(Y) = -------- [ |Y + pb|  -  2|Y|  +  |Y - pb|  ]		(2)
c                 2G
c             2*pb
c
c
c               1               2F       2F           2F
c     r(Z) = -------- [ |Z + pb|  -  2|Z|  +  |Z - pb|  ]		(2)
c                 2F
c             2*pb
c
c  and where `pb' is the length over which the fractional Brownian motion is
c  averaged in order to make this, the derivative process, exist. Normally
c  `pb' is selected to be quite small (of the order of the size of the
c  discretization interval). The parameters `F', `G', and `H' are the
c  self-similar parameters, or Hurst exponents, (0.5 < F,G,H < 1) with
c  F = G = H = 1 giving a total correlated field and F = G = H = 0.5 giving
c  a white noise field (for sufficiently small pb). Finally var is the
c  point variance of the process.
c
c  The parameters var, pb, H, G, and F are brought in via the common
c  block /dparam/.
c
c  If var < 0, then this function computes the variance of a local average
c  of the random field, |var|*V(X,Y,Z), averaged over a domain having side
c  dimensions {X,Y,Z}. For a separable covariance function, the variance
c  function is also separable,
c
c     V(X,Y,Z) = V(X)*V(Y)*V(Z)
c
c  where (X,Y,Z) gives the dimensions of the averaging region. The 1-D
c  variance function corresponding to Eq. (1) is given by
c
c                        r        r      r            r
c                |X + pb| - 2( |X| + |pb| ) + |X - pb|
c       V(X) = -----------------------------------------		(3)
c                      2                  2H
c                     X * (2H+1)*(2H+2)*pb
c
c  where r = (2*H+2). Notice that as X goes to zero the above becomes 0/0, but
c  its limiting value is 1.0. A similar relationship holds for the variance
c  function corresponding to (2).
c
c  Parameters to this process are brought in through the common block
c  /dparam/ and are described as follows;
c
c     var	the point variance of the fGn process.
c
c      pb	the averaging length discussed above. This is assumed common
c		for all three directions.
c
c       H	the Hurst exponent governing the process in the X direction.
c		In this implementation, 0.5 < H < 1; values of H near 1 yield
c		a correlation structure which remains very high (and thus a
c		covariance matrix which may be nearly singular). Values of H
c		near 0.5 yield a band- limited white noise process.
c
c       G	the Hurst exponent governing the process in the Y direction.
c		In this implementation, 0.5 < G < 1; values of G near 1 yield
c		a correlation structure which remains very high (and thus a
c		covariance matrix which may be nearly singular). Values of G
c		near 0.5 yield a band- limited white noise process.
c
c       F	the Hurst exponent governing the process in the Z direction.
c		In this implementation, 0.5 < F < 1; values of F near 1 yield
c		a correlation structure which remains very high (and thus a
c		covariance matrix which may be nearly singular). Values of F
c		near 0.5 yield a band- limited white noise process.
c
c  Arguments to this function are just the X, Y, and Z elements of the
c  lag vector between points in the field for which the covariance is
c  desired (or the dimensions of the physical averaging region, if var < 0).
c
c  REVISION HISTORY:
c  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (May 1/00)
c  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Apr 9/01)
c  1.3	reversed default - now return covariances if var > 0. (Apr 11/01)
c  1.31	revised above documentation to reflect revision 1.3 (May 9/01)
c---------------------------------------------------------------------- 
      real*8 function dlsfr3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, G, F
      data zero/0.d0/, one/1.d0/, two/2.d0/
      abs(q)  = dabs(q)

      if( var .lt. zero ) then			! return variance function
c							X-direction var func
         if( X .eq. zero ) then
            v1 = one
         else
            e0 = two*H
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(X+pb)**e2-two*(abs(X)**e2+pb**e2)+abs(X-pb)**e2
            v1 = t1/(X*X*e1*e2*pb**e0)
         endif
c							Y-direction var func
         if( Y .eq. zero ) then
            v2 = one
         else
            e0 = two*G
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(Y+pb)**e2-two*(abs(Y)**e2+pb**e2)+abs(Y-pb)**e2
            v2 = t1/(Y*Y*e1*e2*pb**e0)
         endif
c							Z-direction var func
         if( Z .eq. zero ) then
            v3 = one
         else
            e0 = two*F
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(Z+pb)**e2-two*(abs(Z)**e2+pb**e2)+abs(Z-pb)**e2
            v3 = t1/(Z*Z*e1*e2*pb**e0)
         endif
c							final variance function
         dlsfr3 = -var*v1*v2*v3
      else					! return covariance
         eH = two*H
         eG = two*G
         eF = two*F
         pH = two*(pb**eH)
         pG = two*(pb**eG)
         pF = two*(pb**eF)
         tH = abs(X+pb)**eH-two*(abs(X)**eH)+abs(X-pb)**eH
         tG = abs(Y+pb)**eG-two*(abs(Y)**eG)+abs(Y-pb)**eG
         tF = abs(Y+pb)**eF-two*(abs(Y)**eF)+abs(Y-pb)**eF
         dlsfr2 = var*tH*tG*tF/(pH*pG*pF)
      endif

      return
      end
