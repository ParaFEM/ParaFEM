c  *******************************************************************
c  *                                                                 *
c  *                     Function dlafs3                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.4
c  Written by Gordon A. Fenton, TUNS, Mar. 15, 1995
c  Latest Update: Jan 11, 2005
c
c  PURPOSE  returns the covariance between two points in a partially
c           separable 3-D fractional Gaussian noise field.
c
c  Returns the covariance between points in a 3-D fractional Gaussian noise
c  process which is isotropic in the X-Y plane. The covariance function is
c  partially separable, and is given by
c
c         B(X,Y,Z) = var*r(X,Y)*r(Z)
c
c  where r(Z) is the correlation function in the z-direction,
c
c               1               2G       2G           2G
c     r(Z) = -------- [ |Z + pb|  -  2|Z|  +  |Z - pb|  ]		(1)
c                 2G
c             2*pb
c
c  while, in the X-Y plane, the isotropic correlation function r(X,Y) is
c  given by
c
c                     1               2H       2H           2H
c  r(X,Y) = r(s) = -------- [ |s + pb|  -  2|s|  +  |s - pb|  ]		(2)
c                       2H
c                   2*pb
c
c  where s = sqrt(X*X + Y*Y), and `pb' is the length over which the fractional
c  Brownian motion is averaged in order to make this, the derivative process,
c  exist. Normally `pb' is selected to be quite small (of the order of the
c  size of the discretization interval). The parameters `G' and `H' are the
c  self-similar parameters, or Hurst exponents, (0.5 < G,H < 1) with
c  G = H = 1 giving a total correlated field and G = H = 0.5 giving a white
c  noise field (for arbitrarily small pb).
c
c  The parameters var, pb, H, and G are brought in through the common block
c  /dparam/. X, Y, and Z are the function arguments and are the elements
c  of the lag vector separating the two points in the field for which the
c  covariance is desired.
c
c  If var < 0, then this function computes the variance of a local average
c  of the random field, averaged over a domain having side dimensions {X,Y,Z}.
c  In this case, the variance function (which is also separable)
c
c     V(X,Y,Z) = V(X,Y)*V(Z)
c
c  is used, and the function returns |var|*V(X,Y,Z). The 1-D
c  variance function corresponding to Eq. (1) is given by
c
c                        r        r      r            r
c                |Z + pb| - 2( |Z| + |pb| ) + |Z - pb|
c       V(Z) = -----------------------------------------		(3)
c                      2                  2H
c                     Z * (2H+1)*(2H+2)*pb
c
c  where r = (2*G+2). Notice that as Z goes to zero the above becomes 0/0, but
c  its limiting value is 1.0.
c
c  The variance function which yields a close approximation to Eq. (2)c  is given by
c
c                      r          r        r              r
c              |D + pb|  -  2( |D|  +  |pb|  )  + |D - pb|
c   V(X,Y) = --------------------------------------------------		(4)
c                          2                  2G
c                         D * (2G+1)*(2G+2)*pb
c
c  where r = (2H+2) and D = |X| + |Y| (the 1-norm of {X,Y}). The actual
c  covariance function corresponding to (4) is equal to (2) along the
c  X,Y coordinate axes, and remains close to (2) for G > 0.8. For G < 0.7,
c  the covariance in the diagonal direction is too high, but since the
c  covariance in such fields dies off very rapidly, the error is
c  negligible.
c  
c  Parameters to this process are brought in through the common block
c  /dparam/ and are described as follows;
c
c     var	the point variance of the fGn process.
c
c      pb	the averaging length discussed above.
c
c       H	the Hurst exponent governing the process in the X,Y plane.
c		In this implementation, 0.5 < H < 1; values of H near 1 yield
c		a correlation structure which remains very high (and thus a
c		covariance matrix which may be nearly singular). Values of H
c		near 0.5 yield a band- limited white noise process.
c
c       G	the Hurst exponent governing the process in the Z direction.
c		In this implementation, 0.5 < G < 1; values of G near 1 yield
c		a correlation structure which remains very high (and thus a
c		covariance matrix which may be nearly singular). Values of G
c		near 0.5 yield a band- limited white noise process.
c
c      da	dummy placeholder provided so that the common block /dparam/
c		has the same form for a variety of variance functions.
c
c  Arguments to this function are just the X, Y, and Z dimensions of the
c  lag distance (or the physical averaging region).
c
c  REVISION HISTORY:
c  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (Apr 29/00)
c  1.11	minor changes to documentation above (Jul 14/00)
c  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Apr 9/01)
c  1.3	reversed default - now return covariances if var > 0 (Apr 11/01)
c  1.31	revised above documentation to reflect revision 1.3 (May 9/01)
c  1.4	made X-Y the isotropic plane (Z is usually vertical) (Jan 11/05)
c---------------------------------------------------------------------- 
      real*8 function dlafs3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, G, da
      data zero/0.d0/, one/1.d0/, two/2.d0/
      abs(q)  = dabs(q)
      sqrt(q) = dsqrt(q)

      if( var .lt. zero ) then			! return variance function
c						Z-direction variance function
         if( Z .eq. zero ) then
            v1 = one
         else
            e0 = two*G
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(Z+pb)**e2 - two*(abs(Z)**e2+pb**e2) + abs(Z-pb)**e2
            v1 = t1/(Z*Z*e1*e2*(pb**e0))
         endif
c						X-Y plane variance function
         if( (X .eq. zero) .and. (Y .eq. zero) ) then
            v2 = one
         else
            D  = abs(Y) + abs(X)
            e0 = two*H
            e1 = e0 + one
            e2 = e0 + two
            t1 = (D+pb)**e2 - two*(D**e2 + pb**e2) + abs(D-pb)**e2
            v2 = t1/(D*D*e1*e2*pb**e0)
         endif
         dlafs3 = -var*v1*v2
      else					! return covariance
c						Z-direction covariance
         if( Z .eq. zero ) then
            v1 = one
         else
            e0 = two*G
            t1 = abs(Z+pb)**e0 - two*(abs(Z)**e0) + abs(Z-pb)**e0
            v1 = t1/(two*(pb**e0))
         endif
c						Y-direction covariance
         if( (Y .eq. zero) .and. (X .eq. zero) ) then
            v2 = one
         else
            e0 = two*H
            D  = sqrt(Y*Y + X*X)
            t1 = (D+pb)**e0 - two*(D**e0) + abs(D-pb)**e0
            v2 = t1/(two*(pb**e0))
         endif
         dlafs3 = var*v1*v2
      endif

      return
      end
