c  *********************************************************************
c  *                                                                   *
c  *                         subroutine hstbkt                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Aug 12, 1997
c  Latest Update: Aug 2, 2006
c
c  PURPOSE  creates histogram bucket centers corresponding to the frequency
c           counts created by histp
c
c  DESCRIPTION
c  This routine takes as input the frequency counts produced by histp and
c  creates the corresponding bucket center locations. The frequency counts
c  are normalized to enclose unit area if nbkt is positive.
c
c  ARGUMENTS
c
c    freq	real vector of length at least nbkt containing the frequency
c		counts produced by histp. On output, these counts are
c		normalized so that the area under the histogram is unity
c		if nbkt is positive. (input/output)
c
c      xf	real vector of length at least nbkt which on output will
c		contain the histogram bucket centers. Normally the buckets
c		are evenly divided between xmin and xmax, however for
c		idist = +1 (lognormal), the buckets are evenly divided in
c		log-space between ln(xmin) and ln(xmax). If idist = 2
c		(exponential distribution) then the buckets widths are
c		computed so that they each have equal area under the
c		assumed distribution. In this case, db is the value
c		of the equi-areas associated with each bucket and the
c		right edge of the i'th bucket, x_{i+1}, in real space can
c		be computed as
c
c			x_{i+1} = -(1/p1)*alog[ ((nbkt-i)*xmn + i*xmx)/nbkt ]
c
c		where xmn = exp(-p1*xmin), xmx = exp(-p1*xmax), and p1 is the
c		lambda parameter of the exponential distribution (one over the
c		mean). If idist corresponds to any of the discrete
c		distributions, then the bucket centers are at integer values
c		as follows;
c		   0, 1, 2, ..., nbkt-1 for idist = 6
c		   1, 2, 3, ..., nbkt   for idist = 7
c		   m, m+1, m+2, ..., m+nbkt-1 for idist = 8
c		   0, 1, 2, ..., nkbt-1 for idist = 9
c		(output)
c
c    nbkt	the number of buckets into which the frequency counts are
c		placed. If nbkt is negative, its absolute value is used,
c		but the histogram is then not normalized. (input)
c
c      nd	integer containing the original total number of observations
c		contained in the frequency count. (input)
c
c   idist	integer type code specifying the type of distribution to
c		use in creating histogram buckets. The value of idist can
c		be one of the following;
c		   = 0	for a normal distribution
c		   = 1	for a lognormal distribution
c			(if idist = -1, then buckets are defined in real
c			space, as for the normal. Otherwise buckets are
c			defined in log-space)
c		   = 2	for an exponential distribution
c		   = 3	for a Beta distribution
c		   = 4  for a Gamma distribution
c		   = 5  for a uniform distribution
c		   = 6  for a binomial (n,p) distribution (discrete)
c		   = 7  for a geometric (p) distribution (discrete)
c		   = 8  for a negative binomial (m,p) distribution (discrete)
c		   = 9  for a Poisson (lambda*t) distribution (discrete)
c		   = 99 if distribution is unknown
c		(input)
c
c    xmin	the minimum value of the data in data space. (input)
c
c    xmax	the maximum value of the data in data space. (input)
c
c      db	width of the intervals used in the frequency counts (bucket
c		widths). If idist = +1, this width is measured in log-space.
c		If idist = 2, then db is the area enclosed under the assumed
c		distribution for each bucket. If idist = 6, 7, 8, or 9, then
c		db = 1, suitable for discrete distributions. (input)
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
c
c  REVISION HISTORY:
c  1.01	corrected above documentation, allowed for non-normalization of the
c	histogram, see nbkt. (Dec 8/00)
c  1.02	db is input, not output as implied above. (Dec 31/01)
c  1.1	added discrete distributions binomial, geometric, neg binom, and
c	Poisson. (Jun 30/05)
c  1.11	Corrected lognormal normalization. (Jul 14/05)
c  1.2	allow for lognormal buckets to be defined in real space (Aug 2/06)
c-------------------------------------------------------------------------
      subroutine hstbkt(freq,xf,nbkt,nd,idist,xmin,xmax,
     >                  db,p1,p2,p3,p4)
      dimension freq(*), xf(*)
      data half/0.5/, one/1.0/

   1  format(3a)
c					check data
      if( nbkt .eq. 0 ) then
         write(6,1)'hstbkt: can''t produce a histogram with 0 buckets!'
         return
      endif
c					number of buckets to use
      nbk = iabs(nbkt)
c					normalize the histogram?
      if( nbkt .lt. 0 ) then
         div = one			! no
      else
         div = one/(db*float(nd))	! yes
      endif
c					compute bucket centers
c					and (optionally) normalize frequencies
      if( idist .eq. 1 ) then		! for the lognormal distribution
         if( nbkt .lt. 0 ) then		! don't normalize the histogram
            do 10 k = 1, nbk
               xf(k)   = xmin*exp( db*(float(k) - half) )
  10        continue
         else				! normalize the histogram
            a = xmin*float(nd)
            do 15 k = 1, nbk
               xf(k)   = xmin*exp( db*(float(k) - half) )
               r       = exp(db*float(k)) - exp(db*float(k-1))
               freq(k) = freq(k)/(a*r)
  15        continue
         endif
      elseif( idist .eq. 2 ) then	! for the exponential distribution
         xmn  = exp(-p1*xmin)
         xmx  = exp(-p1*xmax)
         rbkt = float(nbk)
         xl   = xmin
         do 20 k = 1, nbk
            f = float(k)
            xx = ((rbkt-f)*xmn + f*xmx)/rbkt
            xr = -alog(xx)/p1
            xf(k) = half*(xl + xr)
            if( nbkt .ge. 0 ) freq(k) = freq(k)/(float(nd)*(xr-xl))
            xl = xr
  20     continue
      elseif( (idist .eq. 6) .or. (idist .eq. 9) ) then
         do 30 k = 1, nbk		! for the binomial or Poisson dists
            xf(k)   = float(k-1)
            freq(k) = div*freq(k)
  30     continue
      elseif( idist .eq. 7 ) then	! for the geometric distribution
         do 40 k = 1, nbk
            xf(k)   = float(k)
            freq(k) = div*freq(k)
  40     continue
      elseif( idist .eq. 8 ) then	! for the negative binomial dist
         m = int(p3 + half)
         do 50 k = 1, nbk
            xf(k)   = float(k+m-1)
            freq(k) = div*freq(k)
  50     continue
      else				! for idist = 0, 3, 4, 5, or -1 dists
         do 60 k = 1, nbk		! eg normal, Beta, Gamma, bounded,
            xf(k)   = xmin + db*(float(k)-half)! or lognormal with real buckets
            freq(k) = div*freq(k)
  60     continue
      endif

      return
      end
