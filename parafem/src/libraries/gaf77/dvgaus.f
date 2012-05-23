c  *******************************************************************
c  *                                                                 *
c  *                    subroutine dvgaus                            *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.11
c  Written by Gordon A. Fenton, Princeton, 1989.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  return a vector of N(0,v) independent randomly distributed
c           realizations
c
c  Generates a vector of independent normally distributed random variables
c  having zero mean and prescribed variance. Arguments to the routine are
c  as follows
c
c    G     real vector of length at least N which on output will contain
c          the desired vector of realizations. (output)
c
c    S     real vector of length at least N which contains the prescribed
c          standard deviations for each element of G. (input)
c
c    N     the number of elements in the vector G. (input)
c
c  Notes:
c    1) Ensure that the random number generator, randu (see libGAFsim), is
c	initialized prior to calling this routine. See iseed in libGAFsim.
c    2) since randu returns random numbers in the range (0,1), excluding
c       the endpoints, we need not guard against taking the log of zero
c       below.
c
c  REVISION HISTORY:
c  1.1	now calls libGAFsim's (real*4) randu generator. (Mar 5/99)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      subroutine dvgaus( g, s, n )
      implicit real*8 (a-h,o-z)
      dimension g(*), s(*)
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
         r = sqrt(-two*alog(b))
         g(i)   = r*s(i)*cos(a)
         g(i+1) = r*s(i+1)*sin(a)
  10  continue

      if( n .eq. nn ) return
c                                     for n odd, set the last value
      a = twopi*randu(0)
      b = randu(0)
      r = sqrt(-two*alog(b))
      g(n) = r*s(n)*cos(a)

      return
      end
