c  *********************************************************************
c  *                                                                   *
c  *                         function grcdf                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Dec 5, 1996
c
c  PURPOSE  computes the cumulative distribution function of an n-variate
c           jointly normal random vector having a correlation matrix with
c           a special multiplicative form.
c
c  This routine estimates the multi-dimensional definite integral of the
c  multi-variate normal density function, f(x1,x2,...,xn), having a correlation
c  matrix R with the special form (derived from a vector {g} of length n)
c
c      R_{ij} = g_i*g_j + (1 - g_i*g_i)*d_{ij}
c
c  where d_{ij} is the Chronecker delta and g_i is some constant between 0 and
c  1. The multi-dimensional integral reduces to a one-dimensional integral for
c  such a correlation matrix which is evaluated by fitting a cubic spline at
c  a sequence of (NS + 1) points lying in the range [-14, 14].
c  Arguments to this routine are as follows;
c
c    B     real vector of length at least N containing the upper bounds of the
c          integration. B is overwritten by B(i)/sqrt(1-G(i)**2).
c          (input/output)
c
c    G     real vector of length at least N containing the multiplicative
c          elements of the correlation matrix R (see above). G is overwritten
c          by G(i)/sqrt(1-G(i)**2). (input/output)
c
c    N     integer giving the number of elements in the random vector for
c          which the cumulative distribution is desired. (input)
c
c    ierr  integer flag which is set to zero if all goes well and to -1 if
c          the integral cannot be calculated (usually means that a cubic
c          spline cannot be fitted to the data - one or more data points
c          identical?). (output)
c
c  REVISION HISTORY:
c  1.1	allocated exactly enough memory for temp vector T (Dec 5/96)
c-----------------------------------------------------------------------------
      real function grcdf( B, G, N, ierr )
      dimension B(*), G(*), T(644)
      data zero/0.0/, half/0.50/, one/1.0/, two/2.0/
      data four/4.0/, five/5.0/, seven/7.0/, ten/10.0/
      data fourtn/14.0/, pt2/0.2/, pt05/0.05/
      data rt2pi/2.5066282746310005024/
c					normalize B and G
      do 10 i = 1, N
         r    = sqrt(one - G(i)*G(i))
         B(i) = B(i)/r
         G(i) = G(i)/r
  10  continue
c					set integration points
      T(1) = -fourtn
      T(2) = -ten
      T(3) = -seven
      T(4) = -five
      T(5) = -four
      do 20 i = 5, 14
         T(i) = -four + float(i-4)*pt2
  20  continue
      do 30 i = 15, 94
         T(i) = -two + float(i-14)*pt05
  30  continue
      do 40 i = 95, 104
         T(i) = two + float(i-94)*pt2
  40  continue
      T(105) = five
      T(106) = seven
      T(107) = ten
      T(108) = fourtn
c					set integrand values
      do 50 i = 1, 108
         T(i+108) = exp(-half*T(i)*T(i))/rt2pi
         do 50 j = 1, N
            T(i+108) = T(i+108)*phi(B(i)-G(i)*T(i))
  50  continue
c					compute integral

      grcdf = spint(T(1), T(109), T(217), 108, .true., zero, zero, ierr)

      return
      end
