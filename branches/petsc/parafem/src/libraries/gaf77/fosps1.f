C  *******************************************************************
c  *                                                                 *
c  *                       Real Function fosps1                      *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, 1989
c
c  PURPOSE  returns the power spectral density of a 1-D exponentially
c           damped oscillatory noise.
c
c  Returns the one-sided power spectral density function of an exponentiallly
c  damped oscillatory noise at the given frequency w. The continuous
c  covariance structure for such a process is given by
c
c       B(r) = pa*cos(pb*r)*exp( -2*|r|/tx )
c
c  where pa is the point variance, tx is the scale of fluctuation, and
c  pb is the frequency of the oscillation.
c  pa, pb, and tx are brought into this routine through the common block
c  PARAM. The corresponding one-sided spectral density function G(w) is
c
c                            _                                               _
c                   2*pa    |          1                           1          |
c           G(w) = ------ * | ---------------------  +  --------------------- |
c                   pi*tx   | (2/tx)^2 + (pb - w)^2     (2/tx)^2 + (pb + w)^2 |
c                           -                                                -
c
c  NOTES: 
c      1) parameters of the function are brought over from the calling
c         routine through the common /PARAM/ which has a standard form
c         for all libGAFsim variance external functions. This form is
c         not essential.
c      2) this routine is designed to be used in conjunction with FFT1G
c         which calls it to calculate Fourier coefficient variances.
c =======================================================================
      real function fosps1( w )
      common/param/ pa, pb, tx, ty, tz
      data one/1.0/, two/2.0/
      data pi/3.1415926535897932384/

      a  = two/tx
      aa = a*a
      b  = pb - w
      c  = pb + w

      fosps1 = pa*a*((one/(aa+b*b)) + (one/(aa+c*c)))/pi

      return
      end

