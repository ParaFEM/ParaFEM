c  *******************************************************************
c  *                                                                 *
c  *                     Function dlsfr2                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.51
c  Written by Gordon A. Fenton, TUNS, Sept. 19, 1994
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points in a separable 2-D
c           fractional Gaussian noise random field.
c
c  Returns the covariance between two points, separated by the lag vector
c  {X,Y}, in a random field which is a separable 2-D fractional Gaussian
c  noise. The covariance function for such a field is separable, and given by
c
c         B(X,Y) = var*r(X)*r(Y)
c
c  for two points separated by distance {X,Y}, where r(X) and r(Y) are the
c  directional correlation functions,
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
c  and where `pb' is the length over which the fractional Brownian motion is
c  averaged in order to make this, the derivative process, exist. Normally
c  `pb' is selected to be quite small (of the order of the size of the
c  discretization interval). The parameters `G' and `H' are the
c  self-similar parameters, or Hurst exponents, (0.5 < G,H < 1) with
c  G = H = 1 giving a total correlated field and G = H = 0.5 giving
c  a white noise field (for sufficiently small pb). Finally var is the
c  point variance of the process.
c
c  The parameters var, H, and G are brought in via the common
c  block /dparam/.
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(X,Y), averaged over the domain X x Y.
c  If the covariance function is separable, then so too is the variance
c  function,
c
c     V(X,Y) = V(X)*V(Y)
c
c  where (X,Y) gives the dimensions of the averaging region. The 1-D
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
c		for both directions.
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
c      da	dummy placeholder provided to make the common block consistent
c		amongst a variety of variance functions.
c
c  Arguments to this function are just the X and Y components of the
c  lag vector separating the two points (or the dimensions of the averaging
c  region, if var < 0).
c
c  REVISION HISTORY:
c  1.1	used the same notation in all the fGn process variance functions
c  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (May 1/00)
c  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Mar 27/01)
c  1.4	properly handled sign on var for covariances (Apr 5/01)
c  1.5	reversed default - now return covariances if var > 0. (Apr 11/01)
c  1.51	revised above documentation to reflect revision 1.5 (May 9/01)
c---------------------------------------------------------------------- 
      real*8 function dlsfr2( X, Y )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, G, da
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
c							final variance func
         dlsfr2 = -var*v1*v2
      else					! var < 0, return covariance
         e0 = two*H
         e1 = two*G
         e2 = two*(pb**e0)
         e3 = two*(pb**e1)
         t1 = abs(X+pb)**e0-two*(abs(X)**e0)+abs(X-pb)**e0
         t2 = abs(Y+pb)**e1-two*(abs(Y)**e1)+abs(Y-pb)**e1
         dlsfr2 = var*t1*t2/(e2*e3)
      endif

      return
      end
