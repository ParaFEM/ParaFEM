c  *********************************************************************
c  *                                                                   *
c  *                        real*4 function adfit                      *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.5
c  Written by Gordon A. Fenton, TUNS, Jun  9, 1999
c  Latest Update: Jul 30, 2006
c
c  PURPOSE  computes an Anderson-Darling goodness-of-fit test statistic
c
c  DESCRIPTION
c  This function computes and returns the Anderson-Darling goodness-of-fit
c  test statistic for 8 cases;
c
c	ntype= 1) normal distribution with mean and variance estimated
c	ntype=-1) normal distribution with mean and variance known
c	ntype= 2) exponential distribution with mean estimated
c	ntype=-2) exponential distribution with mean known
c	ntype= 3) Weibull distribution with lambda and beta estimated
c	ntype=-3) Weibull distribution with lambda and beta known
c	ntype= 4) lognormal distribution with mean and variance estimated
c	ntype=-4) lognormal distribution with mean and variance known
c
c  The test on the lognormal fit is simply a test on the normal fit of
c  the log-data. It is assumed, in this case, that p1 and p2 are the
c  mean and standard deviation of the log-data.
c
c  See page 392 of Law and Kelton's ``Simulation Modeling and Analysis,''
c  Second Edition, McGraw-Hill, NY, NY, 1991 for details on the algorithm
c  and critical Anderson-Darling values.
c
c  ARGUMENTS
c
c	x	real vector of length at least n containing the original
c		data (or log-data if the distribution type is lognormal,
c		using ntype = 1). (input)
c
c	z	real vector of length at least n used as workspace.
c
c	n	integer giving the number of data values in x. (input)
c
c   ntype	integer code denoting the distribution which is to be tested
c		against the data;
c		=  1 	normal distribution with parameters estimated
c		= -1 	normal distribution with parameters known
c		=  2 	exponential distribution with parameter estimated
c		= -2 	exponential distribution with parameter known
c		=  3 	Weibull distribution with parameters estimated
c		= -3 	Weibull distribution with parameters known
c		=  4 	lognormal distribution with parameters estimated
c		= -4 	lognormal distribution with parameters known
c		(input)
c
c  reject	logical flag which is true if the null hypothesis (that the
c		data follows the selected distribution) is to be rejected
c		at a significance value of alpha. (ouput)
c
c   alpha	significance value at which to conduct the test (see
c		`reject' above). This must be one of 0.1, 0.05, 0.025,
c		or 0.01. If not, a warning message is issued and the
c		nearest of the above values is selected to perform the
c		test. (input)
c
c      p1	real value containing the first known (or estimated) parameter
c		of the null distribution. In the case of the
c		0) distribution with all parameters known, this is not used
c		1) normal distribution, this is the (sample) mean
c		2) exponential distribution, with F(x) = 1 - exp{-lambda*x},
c		   this is (the estimate of) lambda.
c		3) Weibull distribution, with F(x) = 1 - exp{-(lambda*x)^b}
c		   this is (the estimate of) `lambda'.
c		4) lognormal distribution, this is the (sample) mean of the
c		   log-data
c		(input)
c
c      p2	real value containing the second known (or estimated)
c		parameter of the null distribution. In the case of the
c		0) distribution with all parameters known, this is not used
c		1) normal distribution, this is the (sample) standard deviation
c		2) exponential distribution, with F(x) = 1 - exp{-lambda*x},
c		   this is not used
c		3) Weibull distribution, with F(x) = 1 - exp{-(lambda*x)^b}
c		   this is (the estimate of) `b'.
c		4) lognormal distribution, this is the (sample) standard
c		   deviation of the log-data
c		(input)
c
c    ierr	integer flag denoting the status of the run;
c		= 0  if all goes well
c		= +1 if significance level, alpha, is not recognized
c		= -1 if distribution type is unknown
c
c  REVISION HISTORY:
c  1.01	corrected above documentation slightly. (Jun 15/99)
c  1.1	added the lognormal distribution (Feb 28/00)
c  1.2	now perform all math and sums in double precision (Mar 7/00)
c  1.3	corrected exponential and Weibull calculations (Aug 15/00)
c  1.31 return adfit = 0 if distribution type is unknown (May 24/01)
c  1.4	'pt3' should have been 'pt2' below for Weibull dist (Sep 21/01)
c  1.5	added check, error msg, and shift if any exponential/Weibull/lognormal
c	distributed observations are zero (due to roundoff?). Revise Weibull
c	to exp(-(lambda*x)^beta). (Jul 30/06)
c-------------------------------------------------------------------------
      real*4 function adfit(x,z,n,ntype,reject,alpha,p1,p2,ierr)
      implicit real*8 (a-h,o-z)
      real*4 x(*), z(*), alpha, p1, p2
      real*4 adcrit(4,4), tol
      real*4 pt01, pt0175, pt025, pt0375, pt05, pt075, pt1
      real*4 alog, sngl, abs
      logical reject
      data zero/0.d0/, small/0.001d0/, half/0.5d0/
      data one/1.d0/, four/4.d0/, twnty5/25.d0/
      data pt01/0.01/, pt0175/0.0175/, pt025/0.025/
      data pt0375/0.0375/, pt05/0.05/, pt075/0.075/
      data pt1/0.1/, pt2/0.2d0/, pt6/0.6d0/
      data tol/0.001/
c					Anderson-Darling Critical Values
c           Alpha = 0.10  0.05   0.025  0.01

      data adcrit/1.933, 2.492, 3.070, 3.857,	! all par known
     >            0.632, 0.751, 0.870, 1.029,	! normal
     >            1.070, 1.326, 1.587, 1.943,	! exponential
     >            0.637, 0.757, 0.877, 1.038/	! Weibull


   1  format(3a)
   2  format(a,i3,a)
c					assume all will go well
      ierr = 0
c					transfer the data into temp space
      do 10 i = 1, n
         z(i) = x(i)
  10  continue
c					sort the data via Numerical Recipes
      call sort(z,n)
c					compute adjustment

      dn = dble(n)
      f  = one
      if( (ntype .eq. 1) .or. (ntype .eq. 4) ) then
         f  = one + (four/dn) - (twnty5/(dn*dn))
      elseif( ntype .eq. 2 ) then
         f  = one + pt6/dn
      elseif( ntype .eq. 3 ) then
         f  = one + pt2/dsqrt(dn)
      endif
c					convert to sample distribution
      ad  = zero
      dp1 = dble(p1)
      dp2 = dble(p2)
      nz  = 0
      if( (iabs(ntype) .eq. 1) ) then			! normal dist
         do 20 i = 1, n
            k = n - i + 1
            t1 = (dble(z(i)) - dp1)/dp2
            t2 = (dble(z(k)) - dp1)/dp2
            r1 = dphi(t1)
            r2 = dphi(t2)
            z1 = dlog(r1)
            z2 = dlog(one - r2)
            ad = ad + dble(2*i-1)*(z1 + z2)
  20     continue
      elseif( iabs(ntype) .eq. 2 ) then			! exponential dist
         do 30 i = 1, n
            k = n - i + 1
            t1 = one - dexp(-dp1*dble(z(i)))
            if( t1 .eq. zero ) then
               if( nz .eq. 0 ) then			! only give 1 warning
                 write(6,1)'Warning: zeroes found in data set to which'
                 write(6,1)'         the exponential distribution is'
                 write(6,1)'         being fitted (Anderson-Darling'
                 write(6,1)'         goodness-of-fit test (adfit)).'
                 write(6,1)'         Zeroes will be shifted to the'
                 write(6,1)'         mean/1000. This may bias the test.'
                 write(6,1)'         In addition, this dataset may not'
                 write(6,1)'         be exponentially distributed. Are'
                 write(6,1)'         the zeroes due to rounding?'
               endif
               nz = nz + 1
               t1 = small/dp1
            endif
            z1 = dlog(t1)
            z2 = -dp1*dble(z(k))
            ad = ad + dble(2*i-1)*(z1 + z2)
  30     continue
      elseif( iabs(ntype) .eq. 3 ) then			! Weibull dist
         do 40 i = 1, n
            k = n - i + 1
            t1 = one - dexp( -((dp1*dble(z(i)))**dp2) )
            if( t1 .eq. zero ) then
               if( nz .eq. 0 ) then			! only give 1 warning
                 write(6,1)'Warning: zeroes found in data set to which'
                 write(6,1)'         the Weibull distribution is'
                 write(6,1)'         being fitted (Anderson-Darling'
                 write(6,1)'         goodness-of-fit test (adfit)).'
                 write(6,1)'         Zeroes will be shifted to'
                 write(6,1)'         0.001/lambda. This may bias the'
                 write(6,1)'         test. In addition, this dataset'
                 write(6,1)'         may not be Weibully distributed.'
                 write(6,1)'         Are the zeroes due to rounding?'
               endif
               nz = nz + 1
               t1 = small/dp1
            endif
            z1 = dlog(t1)
            z2 = -((dp1*dble(z(k)))**dp2)
            ad = ad + dble(2*i-1)*(z1 + z2)
  40     continue
      elseif( iabs(ntype) .eq. 4 ) then			! lognormal dist
         do 50 i = 1, n
            k = n - i + 1
            if( (z(i) .eq. zero) .or. (z(k) .eq. zero) ) then
               if( nz .eq. 0 ) then			! only give 1 warning
                 write(6,1)'Warning: zeroes found in data set to which'
                 write(6,1)'         the lognormal distribution is'
                 write(6,1)'         being fitted (Anderson-Darling'
                 write(6,1)'         goodness-of-fit test (adfit)).'
                 write(6,1)'         Zeroes will be shifted to the'
                 write(6,1)'         mean/1000. This may bias the test.'
                 write(6,1)'         In addition, this dataset may not'
                 write(6,1)'         be lognormally distributed. Are'
                 write(6,1)'         the zeroes due to rounding?'
               endif
               nz = nz + 1
               u  = dexp(dp1 + half*dp2*dp2)
               t1 = small*u
            endif
            zl1 = dlog( dble(z(i)) )
            zl2 = dlog( dble(z(k)) )
            t1 = (zl1 - dp1)/dp2
            t2 = (zl2 - dp1)/dp2
            r1 = dphi(t1)
            r2 = dphi(t2)
            z1 = dlog(r1)
            z2 = dlog(one - r2)
            ad = ad + dble(2*i-1)*(z1 + z2)
  50     continue
      else
         write(6,2)'Error: unknown distribution type, ',ntype,
     >             ', in adfit.'
         ierr  = -1
         adfit = 0.0
         return
      endif
c					compute adjusted A-D statistic
      adfit = sngl( -f*(dn + ad/dn) )
c					decode significance level
      if( abs(alpha-pt1) .lt. tol ) then
         i = 1
      elseif( abs(alpha-pt05)  .lt. tol ) then
         i = 2
      elseif( abs(alpha-pt025) .lt. tol ) then
         i = 3
      elseif( abs(alpha-pt01)  .lt. tol ) then
         i = 4
      else
         write(6,1)'Warning: unrecognized significance level in adfit.'
         write(6,1)'         Using next closest.'
         ierr = 1
         if( alpha .gt. pt075 ) then
            i = 1
         elseif( alpha .gt. pt0375 ) then
            i = 2
         elseif( alpha .gt. pt0175 ) then
            i = 3
         else
            i = 4
         endif
      endif
c					decode distribution type
      if( ntype .lt. 0 ) then
         j = 1
      elseif( ntype .eq. 4 ) then
         j = 2
      else
         j = 1 + ntype
      endif
c					perform rejection test
      reject = ( adfit .gt. adcrit(i,j) )
c					all done, return
      return
      end
