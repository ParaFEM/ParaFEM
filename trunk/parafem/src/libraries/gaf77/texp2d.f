c  *********************************************************************
c  *                                                                   *
c  *                     real function texp2d                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c
c  PURPOSE  returns the equivalent 1-d spectral density function corresponding
c           to the 2-d Gauss-Markov process (exponential correlation function)
c           required for the 2-D Turning Bands algorithm.
c
c  This routine takes one argument, frequency, and returns the spectral
c  density value of the 1-d process which is equivalent to the 2-d
c  process having an isotropic exponential correlation structure;
c
c             B(r) = exp(-2*|r|/thx)
c
c  The equivalent 1-d spectral density function (one-sided) is given by
c
c                            2*pa*w*thx^2
c             G_1(w) = -----------------------
c                       (4 + (w*thx)^2)**(1.5)
c
c  Parameters to this function are passed via
c  the common block PARAM: var is the point variance of the process, thx
c  and thy are the directional scales of fluctuation. At this time only the
c  isotropic field is considered, and so all calculations are performed
c  using only thx (assuming thx = thy).
c-----------------------------------------------------------------------------
      real function texp2d( w )
      common/param/ var, pb, thx, thy, thz
      data two/2.0/, four/4.0/

      wt     = w*thx
      texp2d = two*var*wt*thx*((four + wt*wt)**(-1.5))

      return
      end
