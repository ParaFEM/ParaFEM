c  *********************************************************************
c  *                                                                   *
c  *                         subroutine xregrs                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Jul  4, 1997
c  Latest Update: Aug 1, 1997
c
c  PURPOSE  computes regression coefficients given the data sums
c
c  DESCRIPTION
c  Given the sums of the data and the independent variable and their
c  squares, this routine computes the linear regression coefficients
c  of the predictor
c
c	y = a + b*x
c
c  where x is the independent variate. In addition, this routine tests
c  to see if the slope, b, is significantly different than zero at the
c  significance level 'siglvl'. This information is returned as a ratio
c  (bsig) of the test statistic to the critical statistic. If this ratio
c  is greater than 1.0, then we reject the hypothesis that the slope is
c  zero.
c
c  The regression or prediction error variance at any specific value of the
c  independent variate can be calculated via
c
c	s^2 = r0 + r1*x + r2*x^2
c
c  where 'x' is the independent variate and (r0,r1,r2) are returned by this
c  routine. The logical flag 'lpred' controls whether regression or prediction
c  error variance coefficients are returned.
c
c  ARGUMENTS
c
c      sn	real value containing the number of observations in the
c		regression. (input)
c
c      sx	real value containing the sum of the x's (independent
c		variates). (input)
c
c     sx2	real value containing the sum of the squared x's. (input)
c
c      sy	real value containing the sum of the y's (observations taken
c		at each x location). (input)
c
c     sy2	real value containing the sum of the squared y's. (input)
c
c     syx	real value containing the sum of the products (x*y). (input)
c
c	a	regression intercept. (output)
c
c	b	regression slope. (output)
c
c   r0,r1	three real values containing the coefficients to the
c      r2	regression error expression (see above). (output)
c
c  siglvl	significance level at which the test of the coefficient `b' is
c		performed. (input)
c
c    bsig	ratio of test statistic to critical statistic, for significance
c		level `siglvl'. Values of bsig greater than 1.0 imply that the
c		slope is significantly different than zero. (output)
c
c   lpred	logical flag which is true if prediction error variances are
c		to be returned rather than the regression error variances. If
c		lpred is true, then the values of r0, r1, and r2 reflect
c		the prediction error variance (which is larger than
c		the regression error variance). (input)
c
c  REVISION HISTORY:
c  1.1	added flag to control type of error variance to return (Aug 1/97)
c-------------------------------------------------------------------------
      subroutine xregrs(sn,sx,sx2,sy,sy2,syx,a,b,r0,r1,r2,
     >                  siglvl,bsig,lpred)
      implicit real*8 (a-h,o-z)
      logical lpred
      data zero/0.d0/, half/0.5d0/, one/1.d0/, two/2.d0/, big/1.d+10/
      abs(z)     = dabs(z)
      int(z)     = idint(z)
      sqrt(z)    = dsqrt(z)
      studi(z,i) = dstudi(z,i)

      barx = sx/sn
      bary = sy/sn
      Sxx  = sx2 - sn*barx*barx
      Syy  = sy2 - sn*bary*bary
      Sxy  = syx - sn*barx*bary
c						compute a and b
      b = Sxy/Sxx
      a = bary - b*barx
c						compute regression variances
      n = int( half + sn )
      if( n .gt. 2 ) then
         ve = (Syy - b*Sxy)/(sn-two)
      else
         ve = zero
      endif
      if( lpred ) then
         r0 = ve*(one + sx2/(sn*Sxx))		! prediction error
      else
         r0 = ve*sx2/(sn*Sxx)			! regression error
      endif
      r1 = -two*ve*barx/Sxx
      r2 = ve/Sxx
c						test significance of slope (b)

      Ttst = sqrt( ((sn - two)*Sxy*Sxy)/(Sxx*Syy - Sxy*Sxy) )
      asig = one - half*siglvl
      Tcrt = studi(asig,n-2)
      if( Tcrt .eq. zero ) then
         bsig = big
      else
         bsig = abs(Ttst/Tcrt)
      endif
c						all done
      return
      end
