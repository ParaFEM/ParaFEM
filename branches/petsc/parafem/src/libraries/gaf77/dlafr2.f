c  *******************************************************************
c  *                                                                 *
c  *                     Function dlafr2                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.41
c  Written by Gordon A. Fenton, TUNS, Apr. 27, 1994
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points in an isotropic 2-D
c           fractional noise field.
c
c  This function computes and returns the covariance between two points,
c  separated by the distance T, in a 2-D fractional Gaussian noise
c  (fGn) process. The radial covariance function of the fractional
c  noise process is given by
c
c               var              2H       2H           2H
c      B(s) = -------- [ |s + pb|  -  2|s|  +  |s - pb|  ]		(1)
c                  2H
c              2*pb
c
c  where `pb' is the length over which the fractional Brownian motion
c  is averaged in order to make this, the derivative process, exist,
c  and s is the Euclidean length of the lag vector {X,Y}, s = sqrt(X*X+Y*Y).
c  Normally `pb' is selected to be quite small (of the order of the
c  size of the discretization interval).
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(D), averaged over the domain `D' (see below).
c  The variance function which yields a close approximation to the above
c  radial covariance function is given by
c
c                      r          r        r              r
c              |D + pb|  -  2( |D|  +  |pb|  )  + |D - pb|
c      V(D) = --------------------------------------------------	(2)
c                          2                  2H
c                         D * (2H+1)*(2H+2)*pb
c
c  where r = (2H+2) and D = |X| + |Y| (the 1-norm of {X,Y}, where now X and
c  Y are the dimensions of the rectangular averaging region). The actual
c  covariance function corresponding to (2) is equal to (1) along the
c  coordinate axes, and remains close to (1) for H > 0.8. For H < 0.7,
c  the covariance in the diagonal direction is too high, but since the
c  covariance in such fields dies off very rapidly, the error is
c  believed negligible.
c  
c  Parameters to this process are brought in through the common block
c  /dparam/ and are described as follows;
c
c     var	the point variance of the fGn process.
c
c     pb	the averaging length discussed above.
c
c     H		the Hurst exponent. In this implementation, 0.5 < H < 1;
c		values of H near 1 yield a correlation structure which
c		remains very high (and thus a covariance matrix which may
c		be nearly singular). Values of H near 0.5 yield a band-
c		limited white noise process.
c
c     da, db	dummy placeholders provided so that the common block /dparam/
c		has the same form for a variety of variance functions.
c
c  Arguments to this function are just the components, X and Y, of the lag
c  vector separating the two points (or the dimensions of the
c  physical averaging region, if var < 0).
c
c  REVISION HISTORY:
c  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (Apr 29/00)
c  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Mar 27/01)
c  1.3	properly handled sign on var for covariances (Apr 5/01)
c  1.4	reversed default - now return covariances if var > 0 (Apr 11/01)
c  1.41	revised above documentation to reflect revision 1.4 (May 9/01)
c---------------------------------------------------------------------- 
      real*8 function dlafr2( X, Y )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, da, db
      data zero/0.d0/, one/1.d0/, two/2.d0/
      abs(q)  = dabs(q)
      sqrt(q) = dsqrt(q)
c
      e0 = two*H
      if( var .lt. zero ) then			! return variance function
         if( (X. eq. zero) .and. (Y .eq. zero) ) then
            dlafr2 = -var
         else
            D  = abs(X) + abs(Y)
            e1 = e0 + one
            e2 = e0 + two
            t1 = (D+pb)**e2 - two*(D**e2 + pb**e2) + abs(D-pb)**e2
            dlafr2 = -var*t1/(D*D*e1*e2*pb**e0)
         endif
      else					! var < 0, return covariance
         D = sqrt(X*X + Y*Y)
         t1 = (D+pb)**e0 - two*(D**e0) + abs(D-pb)**e0
         dlafr2 = var*t1/(two*(pb**e0))
      endif

      return
      end
