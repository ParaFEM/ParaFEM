c  *******************************************************************
c  *                                                                 *
c  *                    Function gausv                               *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 3.2
c  Written by Gordon A. Fenton, TUNS, Mar. 21, 1993.
c  Latest Update: Dec. 5, 1996
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
c   1) since randu returns random numbers in the range (0,1), excluding
c      the endpoints, we need not guard against taking the log of zero
c      below.
c
c  REVISION HISTORY:
c  3.1	now using new randu function (RAN2 from Numerical Recipes, 2nd Ed)
c	(Oct 14/96)
c  3.2	eliminated unused local variable `zero' (Dec 5/96)
c---------------------------------------------------------------------------
      real function gausv( var )
      real var, a, r
      logical getnxt
      save a, r
      data twopi/6.2831853071795864769/
      data two/2.0/
      data getnxt/.true./

      if( getnxt ) then
         a = twopi*randu(0)
         r = sqrt(-two*var*alog(randu(0)))
         gausv = r*cos(a)
         getnxt = .false.
      else
         gausv = r*sin(a)
         getnxt = .true.
      endif

      return
      end
