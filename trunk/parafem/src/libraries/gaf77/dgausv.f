c  *******************************************************************
c  *                                                                 *
c  *                    Function dgausv                              *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.2
c  Written by Gordon A. Fenton, Princeton, Mar. 19, 1988.
c  Latest Update: Sep 21, 2001
c
c  PURPOSE  return a normally distributed N(0,var) random realization.
c
c  Returns a normally distributed, zero mean, random variable.
c  The argument "var" is the variance of the random number.
c  Ensure that the machine specific random number generator "randu" is
c  initiated in the calling routine with a suitable seed. See `iseed'.
c  Arguments to the function are;
c
c    var  real value prescribing the variance of the random variate. (input)
c
c
c  Notes:
c   1) `randu' is in the libGAFsim library.
c   2) since randu returns random numbers in the range (0,1), excluding
c      the endpoints, we need not guard against taking the log of zero
c      below.
c
c  REVISION HISTORY:
c  1.1	now calls libGAFsim's (real*4) randu generator. (Mar 5/99)
c  1.2	return 'dgausv' rather than 'gausv' (Sep 21/01)
c---------------------------------------------------------------------------
      real*8 function dgausv( var )
      implicit real*8 (a-h,o-z)
      real*4 randu
      logical getnxt
      save a, r, getnxt
      data twopi/6.2831853071795864769d0/
      data two/2.d0/
      data getnxt/.true./
      sqrt(z) = dsqrt(z)
      alog(z) = dlog(z)
      cos(z)  = dcos(z)
      sin(z)  = dsin(z)

      if( getnxt ) then
         a = twopi*randu(0)
         b = randu(0)
         r = sqrt(-two*var*alog(b))
         dgausv = r*cos(a)
         getnxt = .false.
      else
         dgausv = r*sin(a)
         getnxt = .true.
      endif

      return
      end
