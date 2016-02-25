C  *******************************************************************
c  *                                                                 *
c  *                     Real*8 Function dlafr1                      *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.41
c  Written by Gordon A. Fenton, Sept. 1989.
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points in a 1-D fractional
c           Gaussian noise process.
c
c  This function computes and returns the covariance between two points,
c  separated by the distance T, in a 1-D fractional Gaussian noise
c  (fGn) process. The continuous covariance function of the fGn process
c  is approximated by (for lag T)
c
c              var              2H       2H           2H
c     B(T) = -------- [ |T + pb|  -  2|T|  +  |T - pb|  ]		(1)
c                 2H
c             2*pb
c
c  where var is the point variance, H is the self-similarity parameter
c  as defined by Mandelbrot, also called the Hurst exponent, and pb is the
c  length over which the fractional Brownian noise is averaged in order to
c  define the derivative process which is this fGn. Typically pb should be
c  taken as equal to the minimum distance between field points in the
c  simulated process. var, H, and pb are brought into this routine
c  through the common block DPARAM (other parameters in this common are
c  ignored).
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(T), averaged over the distance T.
c  The variance function for the fGn process is given by
c
c                        r        r      r            r
c                |T + pb| - 2( |T| + |pb| ) + |T - pb|
c       V(T) = -----------------------------------------		(4)
c                      2                  2H
c                     T * (2H+1)*(2H+2)*pb
c
c  where r = (2*H+2) and T is the averaging length. Notice that as T goes to
c  zero the above becomes 0/0. However, the limiting value is 1.
c  When var < 0, this function returns the variance of the local average,
c  that is it returns |var|*V(T).
c
c  REVISION HISTORY:
c  1.1	used the same notation in all the fGn process variance functions
c	(Mar 15/95)
c  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (Apr 29/00)
c  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
c  1.4	reversed default - now return covariances if var > 0 (Apr 11/01)
c  1.41	revised above documentation to reflect revision 1.4 (May 9/01)
c =======================================================================
      real*8 function dlafr1( X )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, H, da, db
      data zero/0.d0/, one/1.d0/, two/2.d0/
      abs(x) = dabs(x)

      e0 = two*H
      if( var .lt. zero ) then			! return variance function
         if( X .eq. zero ) then
            dlafr1 = -var
         else
            e1 = e0 + one
            e2 = e0 + two
            t1 = abs(X+pb)**e2-two*(abs(X)**e2+pb**e2)+abs(X-pb)**e2
            dlafr1 = -var*t1/(X*X*e1*e2*pb**e0)
         endif
      else					! return covariance
           t1 = abs(X+pb)**e0 - two*(abs(X)**e0) + abs(X-pb)**e0
           dlafr1 = var*t1/(two*(pb**e0))
      endif

      return
      end
