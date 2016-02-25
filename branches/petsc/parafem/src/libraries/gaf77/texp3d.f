c  *********************************************************************
c  *                                                                   *
c  *                     real function texp3d                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c
c  PURPOSE  returns the equivalent 1-D spectral density function corresponding
c           to the 3-D Gauss-Markov process (exponential correlation function)
c           required for the 3-D Turning Bands algorithm.
c
c  This routine takes one argument, frequency, and returns the spectral
c  density value of the 1-D process which is equivalent to the 3-D
c  process having an isotropic exponential correlation structure;
c
c             B(r) = exp(-2*|r|/tht)
c
c  The equivalent 1-D spectral density function (one-sided) is given by
c
c                            4*pa*tht                4 - (w*tht)^2
c             G_1(w) = ------------------- x [ 1 -  ---------------- ]
c                       pi*(4 + (w*tht)^2)           4 + (w*tht)^2
c
c  Parameters to this function are passed via
c  the common block PARAM: pa is the point variance of the process, thx
c  and thy are the directional scales of fluctuation. At this time only the
c  isotropic field is considered, and so all calculations are performed
c  using only thx (assuming thx = thy).
c-----------------------------------------------------------------------------
      real function texp3d( w )
      common/param/ pa, pb, thx, thy, thz
      data one/1.0/, four/4.0/, pi/3.1415926535897932384/

      wt     = w*thx
      ww     = wt*wt
      fww    = four + ww
      texp3d = four*pa*thx*(one - ((four - ww)/(fww)))/(pi*fww)

      return
      end
