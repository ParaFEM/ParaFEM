C  *******************************************************************
c  *                                                                 *
c  *                     Real*8 Function dlace1                      *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.31
c  Written by Gordon A. Fenton, 1989
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points in a 1-D exponentially
c           damped oscillatory noise process.
c
c  Returns the covariance between two points separated by distance T
c  in an exponentiallly damped oscillatory noise process.
c  The covariance function for such a process is given by
c
c       B(T) = var*cos(dpb*T)*exp( -2*|T|/dthx )
c
c  where var is the point variance, dthx is the scale of fluctuation, and
c  dpb is the frequency of the oscillation, and where var, dpb, and dthx
c  are brought into this routine through the common block DPARAM.
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(T), averaged over the distance T.
c  The variance function V(T) of the process averaged over a certain
c  distance T is calculated according to (see E.H. Vanmarcke, "Random Fields:
c  Analysis and Synthesis", MIT Press, 1984, Chapter 5).
c
c                         T
c       V(T) = (2/T) * INT  (1-s/T)*r(s) ds,
c                         0
c
c  where r(s) is the correlation function, r(s) = B(s)/B(0), so that
c  dlace1(T) = var * V(T).
c
c  REVISION HISTORY:
c  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (Apr 29/00)
c  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Feb 27/01)
c  1.3	reversed default - now return covariances if var > 0 (Apr 11/01)
c  1.31	revised above documentation to reflect revision 1.3 (May 9/01)
c =======================================================================
      real*8 function dlace1( T )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, two/2.d0/
      abs(x) = dabs(x)
      cos(x) = dcos(x)
      sin(x) = dsin(x)
      exp(x) = dexp(x)

      if( var .lt. zero ) then			! return variance function
         if( T .eq. zero ) then
            dlace1 = -var
         else
            ps = two/dthx
            s  = ps*ps + dpb*dpb
            r  = ps*ps - dpb*dpb
            a  = -ps*abs(T)
            b  = a*s + r
            c  = dpb*abs(T)
            d  = -two*var/(T*T*s*s)
            dlace1  = d*(exp(a)*(r*cos(c) - two*ps*dpb*sin(c)) - b)
         endif
      else					! return covariance
         if( dthx .eq. zero ) then
            if( T .eq. zero ) then
               dlace1 = var
            else
               dlace1 = zero
            endif
         else
            dlace1 = var*cos(dpb*T)*exp( -two*abs(T)/dthx )
         endif
      endif

      return
      end
