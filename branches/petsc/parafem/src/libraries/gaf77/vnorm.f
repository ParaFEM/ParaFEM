c  *******************************************************************
c  *                                                                 *
c  *                    subroutine vnorm                             *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 2.21
c  Written by Gordon A. Fenton, TUNS, Mar. 21, 1993.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  return a vector of N(0,1) random distributed realizations
c
c  Returns a normally distributed, zero mean, unit variance random vector
c  in the argument `g'. `n' is the length of the desired vector.
c
c  Notes:
c   1) Ensure that the random number generator is initialized prior to calling
c      this routine.
c   2) since randu returns random numbers in the range (0,1), excluding
c      the endpoints, we need not guard against taking the log of zero
c      below.
c
c  REVISION HISTORY:
c  2.1	now using new randu function (RAN2 from Numerical Recipes, 2nd Ed)
c	(Oct 14/96)
c  2.2	eliminated unused local variables `zero' and `big' (Dec 5/96)
c  2.21	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      subroutine vnorm( g, n )
      real g(*)
      data twopi/6.2831853071795864769/
      data two/2.0/
c
      nn = n/2
      nn = 2*nn
      do 10 i = 1, nn, 2
         a = twopi*randu(0)
         r = sqrt(-two*alog(randu(0)))
         g(i)   = r*cos(a)
         g(i+1) = r*sin(a)
  10  continue

      if( n .eq. nn ) return
c                                     for n odd, set the last value
      a = twopi*randu(0)
      r = sqrt(-two*alog(randu(0)))
      g(n) = r*cos(a)

      return
      end
