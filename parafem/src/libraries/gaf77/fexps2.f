c  ********************************************************************
c  *                                                                  *
c  *                     real function fexps2                         *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c
c  PURPOSE   returns the power spectral density of a 2-D random process having
c            an exponential covariance function.
c
c   Computes the power spectral density at the frequency having components
c   (w1,w2). The random process considered has a simple exponential covariance
c   structure given by;
c
c            B(r1,r2) = PA * exp( -2*sqrt{(r1/thx)^2 + (r2/thy)^2} )
c
c   where PA is the point variance of the process and (thx, thy) are the
c   x- and y-direction scales of fluctuation. These parameters are brought
c   in through the common block /PARAM/.
c   The power spectra for this covariance function is given by the
c   Wiener-Khinchine relationship to be
c
c                                         PA*thx*thy
c               G(w1,w2) = -------------------------------------------
c                          2*pi*(1 + (thx*w1/2)^2 + (thy*w2/2)^2)**1.5
c
c -----------------------------------------------------------------------
      real function fexps2( w1, w2 )
      common/param/ pa, pb, thx, thy, thz
      data pi/3.1415926535897932384/
      data half/0.5/, one/1.0/, onept5/1.5/

      t1 = half*thx*w1
      t2 = half*thy*w2

      fexps2 = (half*pa*thx*thy)/(pi*(one + t1*t1 + t2*t2)**onept5)

      return
      end
