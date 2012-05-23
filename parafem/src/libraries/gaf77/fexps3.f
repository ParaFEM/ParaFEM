c  ********************************************************************
c  *                                                                  *
c  *                     real function fexps3                         *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c  Latest Update: Dec. 5, 1996
c
c  PURPOSE   returns the power spectral density of a 3-D random process having
c            an exponential covariance function.
c
c
c   Computes the power spectral density at the frequency having components
c   (w1,w2,w3). The random process considered has a simple exponential
c   covariance structure given by;
c
c      B(r1,r2,r3) = PA * exp( -2*sqrt{(r1/thx)^2 + (r2/thy)^2 + (r3/thz)^2} )
c
c   where PA is the point variance of the process and (thx, thy, thz) are the
c   x-, y-, and z-direction scales of fluctuation. These parameters are
c   brought in through the common block /PARAM/.
c   The power spectra for this covariance function is given by the
c   Wiener-Khinchine relationship to be
c
c                                          PA*thx*thy*thz
c       G(w1,w2,w3) = --------------------------------------------------------
c                     {pi*[1 + (w1*thx/2)^2 + (w2*thy/2)^2 + (w3*thz/2)^2]}**2
c
c
c  REVISION HISTORY:
c  1.1	eliminated unused local variable `onept5' (Dec 5/96)
c -----------------------------------------------------------------------
      real function fexps3( w1, w2, w3 )
      common/param/ pa, pb, thx, thy, thz
      data pi/3.1415926535897932384/
      data half/0.5/, one/1.0/

      t1 = half*thx*w1
      t2 = half*thy*w2
      t3 = half*thz*w3
      tt = pi*(one + t1*t1 + t2*t2 + t3*t3)

      fexps3 = (pa*thx*thy*thz)/(tt*tt)

      return
      end
