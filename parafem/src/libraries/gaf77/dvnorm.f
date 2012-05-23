c  *******************************************************************
c  *                                                                 *
c  *                    subroutine dvnorm                            *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.21
c  Written by Gordon A. Fenton, Princeton, 1989.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  return a vector of N(0,1) random distributed realizations
c
c  Returns a normally distributed, zero mean, unit variance random vector
c  in the argument `g'. `n' is the length of the desired vector.
c  Ensure that the random number generator is initialized prior to calling
c  this routine.
c
c  REVISION HISTORY:
c  1.1	now uses libGAFsim's randu random number generator (real*4) (Feb 26/98)
c  1.2	eliminates the checks that randu = 0 (since it won't) (Mar 7/99)
c  1.21	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      subroutine dvnorm( g, n )
      implicit real*8 (a-h,o-z)
      dimension g(*)
      real*4 randu
      data twopi/6.2831853071795864769d0/
      data two/2.d0/
      cos(z)  = dcos(z)
      sin(z)  = dsin(z)
      alog(z) = dlog(z)
      sqrt(z) = dsqrt(z)
c
      nn = n/2
      nn = 2*nn
      do 10 i = 1, nn, 2
         a = twopi*randu(0)
         b = randu(0)
         r = sqrt( -two*alog(b) )
         g(i)   = r*cos(a)
         g(i+1) = r*sin(a)
  10  continue

      if( n .eq. nn ) return
c                                     for n odd, set the last value
      a = twopi*randu(0)
      b = randu(0)
      r = sqrt( -two*alog(b) )
      g(n) = r*cos(a)

      return
      end
