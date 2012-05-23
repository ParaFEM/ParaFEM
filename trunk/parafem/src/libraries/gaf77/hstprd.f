c  *********************************************************************
c  *                                                                   *
c  *                         subroutine hstprd                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.8
c  Written by Gordon A. Fenton, TUNS, Aug 12, 1997
c  Latest Update: Aug 2, 2006
c
c  PURPOSE  creates the predicted pdf corresponding to a histogram produced by
c           histp and hstbkt
c
c  DESCRIPTION
c  This routine takes as input a series of locations along the x-axis and the
c  parameters of an assumed distribution and produces the probability density
c  function (or mass) values at the prescribed x locations. It also produces
c  a title for the function.
c
c  ARGUMENTS
c
c    preq	real vector of length at least nbkt containing the true pdf
c		or pmf values at each point specified by xf. (output)
c
c      xp	real vector of length at least nbkt which on output will
c		contain the x locations corresponding to each preq value.
c		This is simply a copy of the vector xf (simply provide xf
c		twice in the calling argument list avoids having to allocate
c		additional memory if the copy is not required). (output)
c
c      xf	real vector of length at least nbkt containing the locations
c		at which the pdf values are desired. (input)
c
c    nbkt	the number of points at which the pdf values are desired.
c		(input)
c
c    tkey	character string of length at least 128 which on output
c		will contain a key string corresponding to the provided
c		pdf (containing it's pertinent parameters and p_crit value
c		if lcrit is true). (output)
c
c     sub	character string containing the character to be used to
c		subscript the parameters (normal, lognormal, and exponential
c		distributions only) provided in the key string.
c		NOTE: if the character is a greek letter, it must be prepended
c		by `%b' which is translated into a leading `\' which, in
c		turn, is then treated by libPSLIB routines as an escape into
c		the greek alphabet. For example, the string `%bq' will
c		eventually (after processing by plotps or display) turn
c		into a greek theta. (input)
c
c   lcrit	logical flag which is true if the p_crit value corresponding
c		to the Chi-Square goodness-of-fit test (see chifit) is to
c		be also placed in the key string. If lcrit is false, pvalue
c		is ignored. (input)
c
c   idist	integer type code specifying the type of distribution to
c		use in creating preq. The value of idist can be one of the
c		following;
c		   = 0	for a normal distribution
c		   = 1	for a lognormal distribution
c		   = 2	for an exponential distribution
c		   = 3	for a Beta distribution
c		   = 4  for a Gamma distribution
c		   = 5  for a uniform distribution
c		   = 6  for a binomial (n,p) distribution (discrete)
c		   = 7  for a geometric (p) distribution (discrete)
c		   = 8  for a negative binomial (m,p) distribution (discrete)
c		   = 9  for a Poisson (lambda*t) distribution (discrete)
c		(input)
c
c      ue	the raw sample average of the data. (input)
c
c      se	the raw sample standard deviation of the data. (input)
c
c      p1	the first parameter of the assumed distribution, p1 takes on
c		the following meanings depending on idist; (input)
c		  if idist = 0, p1 is the mean of the data
c		  if idist = 1, p1 is the mean of the log-data
c		  if idist = 2, p1 is the lambda parameter of the exp. dist.
c		  if idist = 3, p1 is the alpha par. of the Beta [0,1] dist.
c		  if idist = 4, p1 is the alpha par. of the Gamma dist.
c		  if idist = 5, p1 is the lower bound of the distribution
c		  if idist = 6, p1 is the prob of success, p, of the binomial
c		  if idist = 7, p1 is the prob of success, p, of the geometric
c		  if idist = 8, p1 is the prob of success, p, of the neg binom
c		  if idist = 9, p1 is the mean lambda*t of the Poisson dist
c		  if idist = 99, p1 is undefined
c		Notes:
c		  1) in the case of the Beta distribution, idist = 3, the
c		     distribution is assumed to be bounded on the interval
c		     [p3,p4], while the parameters alpha and beta which are
c		     provided in p1 and p2 are those for the transformed
c		     interval [0,1].
c		
c      p2	the second parameter of the assumed distribution, p2 takes on
c		the following meanings depending on idist; (input)
c		  if idist = 0, p2 is the standard deviation (s.d.) of the data
c		  if idist = 1, p2 is the s.d. of the log-data
c		  if idist = 2, p2 is the s.d. of the data
c		  if idist = 3, p2 is the beta parameter of the Beta [0,1] dist
c		  if idist = 4, p2 is the beta parameter of the Gamma dist
c		  if idist = 5, p2 is the upper bound of the distribution
c		  if idist = 6, p2 is the s.d. of the data
c		  if idist = 7, p2 is the s.d. of the data
c		  if idist = 8, p2 is the s.d. of the data
c		  if idist = 9, p2 is the s.d. of the data
c		  if idist = 99, p2 is undefined
c
c      p3	the third parameter of the assumed distribution is defined
c		only for idist = +/- 3, 6 or 8, as follows;
c		  3: if idist = 3 (Beta distribution0 then p3 is estimated
c		     from the data to be the `expected' lower bound,
c		     p3 = xmin - (xmax - xmin)/(n-2). Alternatively, if
c		     idist = -3, then p3 is assumed to contain the prescribed
c		     distribution lower bound on input and it is not changed.
c		  6: if idist = 6 (binomial distribution), then p3 is assumed
c		     to contain the number of trials, n, in the binomial test
c		     on input. Note that since p3 is real, n is determined by
c		     rounding p3 to the nearest integer.
c		  8: if idist = 8 (negative binomial distribution), then p3 is
c		     assumed to contain the number of successes, m, of the
c		     negative binomial test. Note that since p3 is real, n is
c		     determined by rounding p3 to the nearest integer.
c		(input)
c
c      p4	the fourth parameter of the assumed distribution, at this time
c		used only for idist = 3, the Beta distribution: p4 is the
c		upper bound of the distribution. (input)
c
c  pvalue	estimated Chi-squared goodness-of-fit critical significance
c		level or p-value. Ignored unless lcrit is true. (input)
c
c  REVISION HISTORY:
c  1.1	allowed multiple character subscripts (Aug 24/97)
c  1.2	replaced \\'s with two char(92)'s for portability (Dec 1/99)
c  1.3	changed pcrit to pvalue to conform to standard usage (Aug 27/00)
c  1.4	changed '\' escapes to %b for latest version of printv (Dec 31/01)
c  1.41	modified documentation for sub above (Sep 6/02)
c  1.5	corrected computation of the Beta distribution (Feb 16/05)
c  1.6	added discrete distributions binomial, geometric, neg binom, and
c	Poisson (Jun 30/05)
c  1.7	correctly obtain parameter n for idist 6 from p3 (Nov 3/05)
c  1.8	eliminated unused arguments xm and dm (Aug 2/06)
c-------------------------------------------------------------------------
      subroutine hstprd(preq,xp,xf,nbkt,tkey,sub,lcrit,idist,
     >                  ue,se,p1,p2,p3,p4,pvalue)
      dimension preq(*), xp(*), xf(*)
      character*(*) tkey, sub
      character key(2,8)*128, fkey*128
      logical lcrit
      data half/0.5/, one/1.0/
      data rt2pi/2.506628274631000502/
      data key/'%bm_ = %f, %bs_ = %f, p-val = %f',
     >         '%bm_ = %f, %bs_ = %f',
     >         '%bm_{ln } = %f, %bs_{ln } = %f, p-val = %f',
     >         '%bm_{ln } = %f, %bs_{ln } = %f',
     >         '%bl_ = %f, p-val = %f',
     >         '%bl_ = %f',
     >         '%ba = %f, %bb = %f, p-val = %f',
     >         '%ba = %f, %bb = %f',
     >         '%ba = %f, %bb = %f, p-val = %f',
     >         '%ba = %f, %bb = %f',
     >         'L = %f, U = %f, p-val = %f',
     >         'L = %f, U = %f',
     >         'p = %f, p-val = %f',
     >         'p = %f',
     >         '%blt = %f, p-val = %f',
     >         '%blt = %f'/

c					get length of subscript string
      ls = lnblnk(sub)
      ic = 2
      if( lcrit ) ic = 1

      if( idist .eq. 0 ) then		! normal distribution
         if( ls .gt. 0 ) then
            fkey(1:4)             = key(ic,1)(1:4)
            fkey(5:5+ls-1)        = sub(1:ls)
            fkey(5+ls:15+ls)      = key(ic,1)(5:15)
            fkey(15+ls+1:15+2*ls) = sub(1:ls)
            fkey(15+2*ls+1:)      = key(ic,1)(16:)
         else
            fkey(1:4)   = key(ic,1)(1:4)
            fkey(5:5)   = 'X'
            fkey(6:16)  = key(ic,1)(4:15)
            fkey(17:17) = 'X'
            fkey(18:)   = key(ic,1)(16:)
         endif
         if( lcrit ) then
            call prins3(tkey,fkey,p1,p2,pvalue)
         else
            call prins2(tkey,fkey,p1,p2)
         endif
         dv = one/(p2*rt2pi)
         do 10 k = 1, nbkt
            xx      = (xf(k) - ue)/se
            xp(k)   = xf(k)
            preq(k) = dv*exp(-half*xx*xx)
  10     continue
      elseif( idist .eq. 1 ) then	! lognormal distribution
         if( ls .gt. 0 ) then
            fkey(1:8)             = key(ic,2)(1:8)
            fkey(9:9+ls-1)        = sub(1:ls)
            fkey(9+ls:24+ls)      = key(ic,2)(9:24)
            fkey(24+ls+1:24+2*ls) = sub(1:ls)
            fkey(24+2*ls+1:)      = key(ic,2)(25:)
         else
            fkey(1:8)   = key(ic,2)(1:8)
            fkey(9:9)   = 'X'
            fkey(10:25) = key(ic,2)(9:24)
            fkey(26:26) = 'X'
            fkey(27:)   = key(ic,2)(25:)
         endif
         if( lcrit ) then
            call prins3(tkey,fkey,p1,p2,pvalue)
         else
            call prins2(tkey,fkey,p1,p2)
         endif
         dv = one/(p2*rt2pi)
         do 20 k = 1, nbkt
            xx      = (alog(xf(k)) - p1)/p2
            xp(k)   = xf(k)
            preq(k) = dv*exp(-half*xx*xx)/xp(k)
  20     continue

      elseif( idist .eq. 2 ) then	! exponential distribution
         if( ls .gt. 0 ) then
            fkey(1:4) = key(ic,3)(1:4)
            fkey(5:5+ls-1) = sub(1:ls)
            fkey(5+ls:) = key(ic,3)(5:)
         else
            fkey(1:4) = key(ic,3)(1:4)
            fkey(5:5) = 'X'
            fkey(6:)  = key(ic,3)(5:)
         endif
         if( lcrit ) then
            call prins2(tkey,fkey,p1,pvalue)
         else
            call prins1(tkey,fkey,p1)
         endif
         do 30 k = 1, nbkt
            xp(k)   = xf(k)
            preq(k) = p1*exp(-p1*xp(k))
  30     continue
      elseif( idist .eq. 3 ) then	! Beta distribution
         if( lcrit ) then
            call prins3(tkey,key(1,4),p1,p2,pvalue)
         else
            call prins2(tkey,key(2,4),p1,p2)
         endif
         dv = gamln(p1+p2) - gamln(p1) - gamln(p2)
         dv = exp( dv - (p1+p2-one)*alog(p4-p3) )
         a1 = p1 - one
         b1 = p2 - one
         do 40 k = 1, nbkt
            xp(k)   = xf(k)
            preq(k) = dv*((xp(k) - p3)**a1)*((p4 - xp(k))**b1)
  40     continue
      elseif( idist .eq. 4 ) then	! Gamma distribution
         if( lcrit ) then
            call prins3(tkey,key(1,5),p1,p2,pvalue)
         else
            call prins2(tkey,key(2,5),p1,p2)
         endif
         dv = gamln(p1)
         dv = one/(exp( dv )*(p2**p1))
         a1 = p1 - one
         do 50 k = 1, nbkt
            xp(k)   = xf(k)
            preq(k) = dv*(xp(k)**a1)*exp(-xp(k)/p2)
  50     continue
      elseif( idist .eq. 5 ) then	! uniform distribution
         if( lcrit ) then
            call prins3(tkey,key(1,6),p1,p2,pvalue)
         else
            call prins2(tkey,key(2,6),p1,p2)
         endif
         dv = one/(p2 - p1)
         do 60 k = 1, nbkt
            xp(k)   = xf(k)
            preq(k) = dv
  60     continue
      elseif( idist .eq. 6 ) then	! binomial distribution
         if( lcrit ) then
            call prins2(tkey,key(1,7),p1,pvalue)
         else
            call prins1(tkey,key(2,7),p1)
         endif
         n     = int( p3 + half )
         q     = one - p1
         r     = p1/q
         preq(1) = q**n
         xp(1)   = xf(1)
         do 70 k = 2, nbkt
            xp(k)   = xf(k)
            preq(k) = preq(k-1)*r*float(n-k+2)/float(k-1)
  70     continue
      elseif( idist .eq. 7 ) then	! geometric distribution
         if( lcrit ) then
            call prins2(tkey,key(1,7),p1,pvalue)
         else
            call prins1(tkey,key(2,7),p1)
         endif
         q     = one - p1
         do 80 k = 1, nbkt
            xp(k)   = xf(k)
            preq(k) = p1*(q**(k-1))
  80     continue
      elseif( idist .eq. 8 ) then	! negative binomial distribution
         if( lcrit ) then
            call prins2(tkey,key(1,7),p1,pvalue)
         else
            call prins1(tkey,key(2,7),p1)
         endif
         q = one - p1
         m = int( p3 + half )
         j = m
         preq(1) = p1**m
         xp(1)   = xf(1)
         do 90 k = 2, nbkt
            preq(k) = preq(k-1)*q*float(j)/float(k-1)
            xp(k) = xf(k)
            j     = j + 1
  90     continue
      elseif( idist .eq. 9 ) then	! Poisson distribution
         if( lcrit ) then
            call prins2(tkey,key(1,8),p1,pvalue)
         else
            call prins1(tkey,key(2,8),p1)
         endif
         preq(1) = exp(-p3)
         xp(1)   = xf(1)
         do 100 k = 2, nbkt
            preq(k) = preq(k-1)*p3/float(k-1)
            xp(k)   = xf(k)
 100     continue
      endif

      return
      end
