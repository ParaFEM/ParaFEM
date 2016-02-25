c  *********************************************************************
c  *                                                                   *
c  *                       real*8 function derf                        *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Mar  7, 1999
c  Latest Update: Mar  7, 1999
c
c  PURPOSE  returns the error function
c
c  DESCRIPTION
c  This function returns the error function defined by
c
c                       x
c     erf(x) = (2/pi) int exp( -z^2 ) dz
c                       0
c
c  It does so by employing the cumulative standard normal distribution.
c  See dphi.f
c
c  ARGUMENTS
c
c	x	real value giving the upper bound on the integrand. If x
c		is negative, -erf(-x) is returned. (input)
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      real*8 function derf(x)
      implicit real*8 (a-h,o-z)
      data rttwo/1.4142135623730950488d0/
      data zero/0.d0/, one/1.d0/, two/2.d0/

      derf = two*dphi(rttwo*x) - one

      return
      end
