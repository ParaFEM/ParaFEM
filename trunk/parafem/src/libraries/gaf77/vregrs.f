c  *********************************************************************
c  *                                                                   *
c  *                         subroutine vregrs                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.3
c  Written by Gordon A. Fenton, TUNS, Tue May 13 19:25:15 1997
c  Latest Update: Aug 1, 1997
c
c  PURPOSE  performs a simple linear regression on a set of data
c
c  This program fits a simple straight line to a set of (x,y) data. In
c  addition, the slope coefficient is tested to see if it is significantly
c  different than zero (at significance level `siglvl'). The trend determined
c  for the data is optionally removed from the data before returning. This
c  routine fits the line
c
c		y = a + b*x						(1)
c
c  where x is the independent variate. This line is fit to the (x,y) data
c  via least squares regression.
c
c  Arguments to this routine are as follows;
c
c	x	real vector of length at least n which contains the
c		independent abscissa values. (input)
c
c	y	real vector of length at least n which contains the dependent
c		ordinate data to be fitted. (input)
c
c	n	the number of data points. (input)
c
c	a	regression intercept. (output)
c
c	b	regression slope. (output)
c
c    bsig	ratio of test statistic to critical statistic, for significance
c		level `siglvl'. Values of bsig greater than 1.0 imply that the
c		slope is significantly different than zero. (output)
c
c  siglvl	significance level at which the test of the coefficient `b' is
c		performed. (input)
c
c    avea	real value which will be used to compute the average regression
c		intercept on the nit'th call to this routine. Between calls to
c		this routine, avea should not be changed. (input/output)
c
c    aveb	real value which will be used to compute the average regression
c		slope on the nit'th call to this routine. Between calls to
c		this routine, aveb should not be changed. (input/output)
c
c    avet	real value which will be used to compute the average slope
c		significance ratio on the nit'th call to this routine. Between
c		calls to this routine, avet should not be changed.
c		(input/output)
c
c   r0,r1	three real values containing running sums of the regression
c      r2	variance estimates used to compute a confidence interval on
c		the regression given by avea and aveb. These sums are over
c		the nit calls to this routine and so these parameters should
c		not be changed between calls to this routine. (input/output)
c
c   globa	real value which will contain the regression intercept obtained
c		by regressing on all data over the nit calls to this routine.
c		It is only computed when nen = nit. (output)
c
c   globb	real value which will contain the regression slope obtained
c		by regressing on all data over the nit calls to this routine.
c		It is only computed when nen = nit. (output)
c
c   g0,g1	three real values containing the global regression variance
c      g2	estimates used to compute a confidence interval on the
c		regression given by globa and globb. (output)
c
c    gsig	ratio of slope test statistic to critical statistic, for
c		significance level `siglvl'. Values of gsig greater than 1.0
c		imply that the globally regressed slope is significantly
c		different than zero. (output)
c
c    nsig	integer containing a running sum of the number of times the
c		slope was determined to be significant (bsig > 1) over the
c		nit calls to this routine.
c		(input/output)
c
c     nen	entry number - this is the number of times that this routine
c		has been called. When nen = nit, then the average regression
c		over all calls and the global regression are computed prior
c		to returning. (input)
c
c     nit	the total number of realizations over which the calculations
c		are to be performed. See nen for details. (input)
c
c   lpred	logical flag which is true if prediction error variances are
c		to be returned rather than the regression error variances. If
c		lpred is true, then the values of r0, r1, r2, and g0, g1, g2,
c		reflect the prediction error variance (which is larger than
c		the regression error variance). (input)
c
c   istat	unit number to which error messages are to be logged. (input)
c
c  verbos	logical flag which is true if error messages can be logged to
c		standard output. (input)
c
c  REVISION HISTORY:
c  1.1	modified to be driven by RGRS1D, fewer arguments, added CI stuff
c	(Jun 19/97)
c  1.2	added global regression over ensemble, moved some stuff to xregrs
c	(Jul 4/97)
c  1.3	added flag controlling computation of prediction vs regression error
c	(Aug 1/97)
c-------------------------------------------------------------------------
      subroutine vregrs(x,y,n,a,b,
     >                  bsig,siglvl,avea,aveb,avet,r0,r1,r2,
     >                  globa,globb,g0,g1,g2,gsig,nsig,nen,nit,
     >                  lpred, istat,verbos)
      implicit real*8 (a-h,o-z)
      dimension x(*), y(*)
      logical verbos, lpred
      save gsn, gsx, gsx2, gsy, gsy2, gsyx
      data one/1.d0/
      float(i) = dble(i)
c						formats
   1  format(5a)
c						compute regression
      if( n .lt. 2 ) then
         if( verbos ) write(6,1)'Cannot compute regression with less tha
     >n 2 data points.'
         write(istat,1)'Cannot compute regression with less than 2 data 
     >points.'
         return
      endif

      sn  = float(n)
      sx  = x(1)
      sx2 = x(1)*x(1)
      sy  = y(1)
      sy2 = y(1)*y(1)
      syx = x(1)*y(1)
      do 10 i = 2, n
         sx  = sx  + x(i)
         sx2 = sx2 + x(i)*x(i)
         sy  = sy  + y(i)
         sy2 = sy2 + y(i)*y(i)
         syx = syx + x(i)*y(i)
  10  continue
c						individual calcs
      call xregrs(sn,sx,sx2,sy,sy2,syx,a,b,s0,s1,s2,
     >            siglvl,bsig,lpred)
c						ensemble calcs
      if( nen .le. 1 ) then
         gsn  = sn
         gsx  = sx
         gsx2 = sx2
         gsy  = sy
         gsy2 = sy2
         gsyx = syx
         avea = a
         aveb = b
         avet = bsig
         r0   = s0
         r1   = s1
         r2   = s2
         if( bsig .gt. one ) then
            nsig = 1
         else
            nsig = 0
         endif
      else
         gsn  = gsn  + sn
         gsx  = gsx  + sx
         gsx2 = gsx2 + sx2
         gsy  = gsy  + sy
         gsy2 = gsy2 + sy2
         gsyx = gsyx + syx
         avea = avea + a
         aveb = aveb + b
         avet = avet + bsig
         r0   = r0 + s0
         r1   = r1 + s1
         r2   = r2 + s2
         if( bsig .gt. one ) nsig = nsig + 1
      endif
      if( (nen .eq. nit) .and. (nit .gt. 1) ) then
         call xregrs(gsn,gsx,gsx2,gsy,gsy2,gsyx,globa,globb,
     >               g0,g1,g2,siglvl,gsig,lpred)
         avea = avea/float(nit)
         aveb = aveb/float(nit)
         avet = avet/float(nit)
         snn  = float(nit*nit)
         r0   = r0/snn
         r1   = r1/snn
         r2   = r2/snn
      elseif( nit .eq. 1 ) then
         globa = a
         globb = b
         g0    = s0
         g1    = s1
         g2    = s2
         gsig  = bsig
      endif
c						all done
      return
      end
