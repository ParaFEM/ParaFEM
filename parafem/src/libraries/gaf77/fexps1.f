c  ********************************************************************
c  *                                                                  *
c  *                     real function fexps1                         *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c
c  PURPOSE   returns the power spectral density function of a 1-D random
c            process having an exponential covariance function.
c
c  This function computes the power spectral density at the
c  frequency/wavenumber given by the argument `w'. The random process
c  is one having a simple exponential covariance function of the form;
c
c               B(r) = PA * exp( -2*|r|/thx )
c
c  where PA is the point variance of the process and thx is the
c  scale of fluctuation. These parameters are brought in through the
c  common block /PARAM/.
c  The (one-sided) power spectra for this covariance function is given
c  through the Wiener-Khinchine relationship to be
c
c                          PA*thx
c               G(w) = ----------------
c                      pi*(1 + (thx*W/2)^2)
c
c -----------------------------------------------------------------------
      real function fexps1( w )
      common/param/ pa, pb, thx, thy, thz
      data pi/3.1415926535897932384/
      data half/0.5/, one/1.0/

      tt = thx*w*half
      fexps1 = (pa*thx)/(pi*(one + tt*tt))

      return
      end
