c  *********************************************************************
c  *                                                                   *
c  *                       subroutine histp                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.8
c  Written by Gordon A. Fenton, TUNS, May 10, 1994
c  Latest Update: Aug 2, 2005
c
c  PURPOSE  to compute a histogram for a given sequence of observations
c           and set up associated ideal cumulative distribution.
c
c  This routine accepts a given sequence "x" of "npt" observations and
c  determines the frequency count in each of "nbkt" buckets uniformly
c  spaced from the minimum to the maximum of the data. The frequency
c  counts are thus determined in the nbkt intervals
c
c    [xmin,xmin+db], (xmin+db,xmin+2*db], ... [xmax-db,xmax]
c
c  where db is the bucket width, with the following exceptions;
c  1) If the data comes from a lognormal distribution and idist = +1
c     (see idist below), then the buckets are uniformly spaced between
c     ln(xmin) and ln(xmax) in log space.
c  2) If the data comes from an exponential distribution, then the
c     buckets widths are computed so that they each have equal area
c     under the assumed distribution. In this case, db is the value
c     of the equi-areas associated with each bucket and the right edge
c     of the i'th bucket, x_{i+1}, in real space can be computed as
c
c	x_{i+1} = -(1/p1)*alog[ ((nbkt-i)*xmn + i*xmx)/nbkt ]
c
c     where xmn = exp(-p1*xmin), xmx = exp(-p1*xmax), and p1 is the
c     lambda parameter of the exponential distribution (one over the
c     mean).
c  3) If the data comes from one of the discrete distributions, idist = 6,
c     7, 8, 9, or 98, then
c	a) if idist = 6 (binomial, with parameters n and p), then nbkt would
c	   normally be set equal to n+1 by the calling routine (see also p3)
c	   Buckets are assumed to correspond to i = 0, 1, ..., nbkt-1.
c	b) if idist = 7 (geometric, with parameter p), then nbkt would
c	   normally be set equal xmax. Buckets are assumed to correspond
c	   to i = 1, 2, ..., nbkt.
c	c) if idist = 8 (negative binomial, with parameters m and p), then
c	   nbkt would normally be set equal to [xmax - m + 1] by the calling
c	   routine. Buckets are assumed to correspond to i = m, m+1, ...,
c	   [nbkt+m-1].
c	d) if idist = 9 (Poisson, with parameter p3=lambda*t), then nbkt would
c	   normally be set equal to xmax+1. Buckets are assumed to correspond
c	   to i = 0, 1, ..., nbkt-1.
c	e) if idist = 98 (unknown discrete dist), then nbkt would normally be
c	   set equal to [xmax - xmin + 1]. Buckets are assumed to correspond
c	   to i = xmin, xmin+1, ..., xmax.
c
c  Assuming the data comes from a known distribution (see idist below), the
c  ideal cumulative probabilities associated with the left edge of the first
c  bucket followed by the right edge of the first, second, ..., nbkt'th bucket
c  are returned in the vector `Fx' (thus Fx is a vector of length (nbkt+1)).
c  Arguments to this routine are as follows;
c
c	x	real vector of length at least npt containing the raw data.
c		(input)
c
c    freq	real vector of length at least nbkt which on output will
c		contain the frequency counts associated with the intervals
c		[xmin,xmin+db], (xmin+db,xmin+2*db], ..., (xmax-db,xmax],
c		for a total of nbkt intervals, where
c		   db = (xmax - xmin)/nbkt
c		Note that if idist = 1, then the buckets and their widths
c		are computed in log-space. In this case, the intervals are
c		  [alog(xmin), alog(xmin)+db], ... (alog(xmax)-db,alog(xmax)]
c		with db = (alog(xmax) - alog(xmin))/nbkt. If idist = 2, then
c		see note 2 above. If idist = 6,7,8, 9, or 98, then see note 3
c		above. (output)
c
c    xmin	the minimum value of the data in data space. (input)
c
c    xmax	the maximum value of the data in data space. (input)
c
c     npt	the number of data values in x. (input)
c
c    nbkt	the number of buckets into which the frequency counts are
c		placed. For continuous distributions, db is adjusted so that
c		the range between xmin and xmax is suitably divided into
c		nbkt buckets. However, for discrete distributions, db = 1,
c		so that nbkt must be properly set by the calling routine to
c		reflect the number of discrete values of interest. See notes
c		above. (input)
c
c      db	width of the intervals used in the frequency counts (bucket
c		widths). If idist = 1, this width is measured in log-space.
c		If idist = 2, then db is the area enclosed under the assumed
c		distribution for each bucket. If idist = 6, 7, 8, 9, or 98,
c		then db = 1, suitable for discrete distributions. (output)
c
c   lprob	logical flag which is true if the cumulative probabilities
c		associated with the left edge of the first bucket and the
c		right edge of each bucket are to be computed (assuming the
c		data comes from a `known' distribution). If lprob is false,
c		the remaining arguments to this routine are ignored. (input)
c
c   idist	integer type code specifying the type of distribution to
c		use in creating Fx. The value of idist can be one of the
c		following;
c		   = 0	for a normal distribution
c		   = 1	for a lognormal distribution
c			(if idist = -1, then buckets are computed in real
c			space, as for the normal. Otherwise buckets are
c			defined in log-space, see comments above)
c		   = 2	for an exponential distribution
c		   = 3	for a Beta distribution
c		   = 4  for a Gamma distribution
c		   = 5  for a uniform distribution
c		   = 6  for a binomial (n,p) distribution (discrete)
c		   = 7  for a geometric (p) distribution (discrete)
c		   = 8  for a negative binomial (m,p) distribution (discrete)
c		   = 9  for a Poisson (lambda*t) distribution (discrete)
c		   = 98 if distribution is discrete but unknown
c			(Fx is not produced)
c		   = 99 if distribution is continuous but unknown
c			(Fx is not produced)
c		(input)
c
c      Fx	real vector of length at least (nbkt + 1) which on output
c		will contain the idealized cumulative probabilities associated
c		with the left edge of the first bucket and the right edge of
c		each bucket (that is, Fx(1) contains the cumulative
c		probability at the left edge of bucket 1, Fx(2) contains the
c		cumulative probability at the right edge of bucket 1, Fx(3)
c		contains the cumulative probability at the right edge of
c		bucket 2, and so on to Fx(nbkt+1) containing the cumulative
c		probability at the right edge of bucket nbkt. This vector is
c		only produced if the distribution type is known (see idist).
c		(output)
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
c		  if idist = 98, p1 is undefined
c		  if idist = 99, p1 is undefined
c		Notes:
c		  1) in the case of the Beta distribution (idist = 3) the
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
c		  if idist = 98, p2 is undefined
c		  if idist = 99, p2 is undefined
c
c      p3	the third parameter of the assumed distribution is defined
c		only for idist = 3, 6 or 8, as follows;
c		  3: if idist = 3 (Beta distribution) then p3 is the
c		     distribution lower bound.
c		  6: if idist = 6 (binomial distribution), then p3 is assumed
c		     to contain the number of trials, n, in the binomial test
c		     on input. Note that since p3 is real, n is determined by
c		     rounding p3 to the nearest integer.
c		  8: if idist = 8 (negative binomial distribution), then p3 is
c		     assumed to contain the number of successes, m, of the
c		     negative binomial test. Note that since p3 is real, m is
c		     determined by rounding p3 to the nearest integer.
c		(input)
c
c      p4	the fourth parameter of the assumed distribution, at this time
c		defined only for idist = 3 (Beta distribution): p4 is the
c		distribution upper bound. (input)
c
c      xm	takes on the following meanings depending on idist (output);
c		  if idist = 0, xm is the number of standard deviations xmin
c				is from the mean u
c		  if idist = 1, xm is the number of standard deviations
c				ln(xmin) is from the mean in log-space
c		  if idist =-1, xm is xmin
c		  if idist = 2, xm is the mean failure rate (= p1)
c		  if idist = 3, xm is the mean derived from alpha and beta
c				and scaled to the data range
c		  if idist = 4, xm is the mean derived from alpha and beta
c		  if idist = 5, xm is the width of the dist. (= p2 - p1)
c
c      dm	takes on the following meanings depending on idist (output);
c		  if idist = 0, dm is the number of standard deviations the
c				bucket width is
c		  if idist = 1, dm is the number of standard deviations the
c				bucket width is in log-space
c		  if idist =-1, dm is the bucket width
c		  if idist = 2, dm is not defined
c		  if idist = 3, dm is the standard deviation derived from
c				alpha and beta and scaled to the data range
c		  if idist = 4, dm is the standard deviation derived from
c				alpha and beta
c		  if idist = 5, dm is not defined
c
c  REVISION HISTORY:
c  1.1	improved the handling of interval boundaries.
c  1.2	can now directly specify one of a set of `known' distributions
c  1.3	separated mean and sd from distribution parameters.
c  1.31	corrected Gamma dist calculation of Fx.
c  1.32 corrected histogram for exponential.
c  1.4	implemented `equal area' buckets for the exponential (Oct 7/95)
c  1.5	eliminated unused arguments `u' and `s' (Dec 5/96)
c  1.6	added discrete distributions binomial, geometric, neg binom, and
c	Poisson (Jun 26/05)
c  1.7	correctly obtain parameter n for idist 6 from p3 (Nov 3/05)
c  1.8	allow lognormal buckets to be set in real space (Aug 2/06)
c----------------------------------------------------------------------------
      subroutine histp( x, freq, xmin, xmax, npt, nbkt, db, lprob,
     >                  idist, Fx, p1, p2, p3, p4, xm, dm )
      real x(*), freq(*), Fx(*)
      real*8 dbeta, dgmdst, du, ds, dt, dx, dexp, dble
      logical lprob
      data zero/0.0/, half/0.5/, one/1.0/, eps/1.e-06/

   1  format(a)
c						initialize frequency count
      do 10 i = 1, nbkt
         freq(i) = zero
  10  continue
      de  = float(nbkt)*eps
c						build histogram (and prob's)
      if( idist .eq. 0 ) then				! normal dist
         db  = (xmax - xmin)/float(nbkt)
         dbi = one/db
         xm  = (xmin - p1)/p2
         dm  = db/p2
         do 20 i = 1, npt
            m       = 1 + int( dbi*(x(i)-xmin) - de )
            freq(m) = freq(m) + one
  20     continue
         if( lprob ) then
            Fx(1) = phi(xm)
            do 30 k = 1, nbkt
               Fx(k+1) = phi(xm + dm*float(k))
  30        continue
         endif
      elseif( iabs(idist) .eq. 1 ) then			! lognormal dist
         if( idist .eq. -1 ) then		! use real-space buckets
            db  = (xmax - xmin)/float(nbkt)
            dbi = one/db
            xm  = xmin
            dm  = db
            do 40 i = 1, npt
               m       = 1 + int( dbi*(x(i)-xmin) - de )
               freq(m) = freq(m) + one
  40        continue
            if( lprob ) then
               if( xmin .le. zero ) then
                  Fx(1) = zero
               else
                  Fx(1) = phi((alog(xmin)-p1)/p2)
               endif
               do 50 k = 1, nbkt
                  xm = xmin + db*float(k)
                  Fx(k+1) = phi((alog(xm)-p1)/p2)
  50           continue
            endif
         else					! use log-space buckets
            if( xmax .le. zero ) then
             write(6,1)'*** histp: non-positive xmax in lognormal dist.'
             write(6,1)'           histogram not produced.'
             return
            endif
            if( xmin .le. zero ) then
             write(6,1)'*** histp: non-positive xmin in lognormal dist.'
             write(6,1)'           histogram not produced.'
             return
            endif
            xmx = alog(xmax)
            xmn = alog(xmin)
            db  = (xmx - xmn)/float(nbkt)
            dbi = one/db
            xm  = (xmn - p1)/p2
            dm  = db/p2
            do 60 i = 1, npt
               m       = 1 + int( dbi*(alog(x(i))-xmn) - de )
               freq(m) = freq(m) + one
  60        continue
            if( lprob ) then
               Fx(1) = phi(xm)
               do 70 k = 1, nbkt
                  Fx(k+1) = phi(xm + dm*float(k))
  70           continue
            endif
         endif
      elseif( idist .eq. 2 ) then			! exponential dist.
         if( xmax .lt. zero ) then
            write(6,1)'*** histp: negative xmax in exponential dist.'
            write(6,1)'           histogram not produced.'
            return
         endif
         if( xmin .lt. zero ) then
            write(6,1)'*** histp: negative xmin in exponential dist.'
            write(6,1)'           histogram not produced.'
            return
         endif
         xmn = exp(-p1*xmin)
         db  = (xmn - exp(-p1*xmax))/float(nbkt)
         dbi = one/db
         xm  = p1
         do 80 i = 1, npt
            m = 1 + int( dbi*(xmn - exp(-p1*x(i))) - de)
            freq(m) = freq(m) + one
  80     continue
         if( lprob ) then
            dx = dble(p1)*dble(xmin)
            du = dexp( -dx )
            dx = dble(p1)*dble(xmax)
            ds = dexp( -dx )
            dx = dble(nbkt)
            Fx(1) = sngl( 1.d0 - du )
            do 90 k = 1, nbkt
               Fx(k+1) = sngl( (dx - dble(nbkt-k)*du - dble(k)*ds)/dx )
  90        continue
         endif
      elseif( idist .eq. 3 ) then			! Beta distribution
         r   = p4 - p3					! [p3,p4] --> [0,1]
         db  = (xmax - xmin)/float(nbkt)
         dbi = one/db
         us  = p1 + p2
         xm  = p3 + r*p1/us
         dm  = r*sqrt( (p1*p2)/(us*us*(us + one)) )
         do 100 i = 1, npt
            m       = 1 + int( dbi*(x(i)-xmin) - de )
            freq(m) = freq(m) + one
 100     continue
         if( lprob ) then
            dt    = 1.d0/(dble(p4)-dble(p3))
            du    = dble(p1)
            ds    = dble(p2)
            do 110 k = 1, nbkt + 1
               dx = dt*(dble(xmin) + dble(k-1)*dble(db) - dble(p3))
               Fx(k) = dbeta(dx,du,ds)
 110        continue
         endif
      elseif( idist .eq. 4 ) then			! Gamma distribution
         if( xmax .lt. zero ) then
            write(6,1)'*** histp: negative xmax in Gamma dist.'
            write(6,1)'           histogram not produced.'
            return
         endif
         if( xmin .lt. zero ) then
            write(6,1)'*** histp: negative xmin in Gamma dist.'
            write(6,1)'           histogram not produced.'
            return
         endif
         db  = (xmax - xmin)/float(nbkt)
         dbi = one/db
         xm  = p1*p2
         dm  = p2*sqrt(p1)
         do 120 i = 1, npt
            m       = 1 + int( dbi*(x(i)-xmin) - de )
            freq(m) = freq(m) + one
 120     continue
         if( lprob ) then
            du    = dble(p1)
            ds    = dble(p2)
            dt    = dble(xmin)
            Fx(1) = dgmdst(dt,du,ds)
            do 130 k = 1, nbkt
               dt = dble(xmin) + dble(k)*dble(db)
               Fx(k+1) = dgmdst(dt,du,ds)
 130        continue
         endif
      elseif( idist .eq. 5 ) then			! uniform distribution
         xm  = p2 - p1
         db  = (xmax - xmin)/float(nbkt)
         dbi = one/db
         do 140 i = 1, npt
            m       = 1 + int( dbi*(x(i)-xmin) - de )
            freq(m) = freq(m) + one
 140     continue
         if( lprob ) then
            Fx(1) = (xmin - p1)/xm
            do 150 k = 1, nbkt
               Fx(k+1) = (xmin + db*float(k) - p1)/xm
 150        continue
         endif
c GOTOHERE: need to verify the following carefully...
c	    add error messages if p = 0.0 or p = 1.0?
      elseif( idist .eq. 6 ) then			! binomial distribution
         db = one
         do 160 i = 1, npt
            m = 1 + int(x(i) + half)
            if( m .le. nbkt ) freq(m) = freq(m) + one
 160     continue
         if( lprob ) then
            n     = int( p3 + half )
            q     = one - p1
            r     = p1/q
            prb   = q**n
            Fx(1) = zero
            Fx(2) = prb
            do 170 k = 3, nbkt+1
               prb = prb*r*float(n-k+3)/float(k-2)
               Fx(k) = Fx(k-1) + prb
 170        continue
         endif
      elseif( idist .eq. 7 ) then			!geometric distribution
         db   = one
         do 180 i = 1, npt
            m = int(x(i) + half)
            if( m .ge. 1 .and. m .le. nbkt ) freq(m) = freq(m) + one
 180     continue
         if( lprob ) then
            q     = one - p1
            Fx(1) = zero
            Fx(2) = p1
            do 190 k = 3, nbkt+1
               Fx(k) = Fx(k-1) + p1*(q**(k-2))
 190        continue
         endif
      elseif( idist .eq. 8 ) then			!neg binom distribution
         m  = int( p3 + half )
         db = one
         do 200 i = 1, npt
            j = 1 + int( x(i) + half ) - m
            if( j .le. nbkt ) freq(j) = freq(j) + one
 200     continue
         if( lprob ) then
            q     = one - p1
            prb   = p1**m
            Fx(1) = zero
            Fx(2) = prb
            j = m
            do 210 k = 3, nbkt+1
               prb   = prb*q*float(j)/float(k-2)
               Fx(k) = Fx(k-1) + prb
               j     = j + 1
 210        continue
         endif
      elseif( idist .eq. 9 ) then			! Poisson distribution
         db = one
         do 220 i = 1, npt
            m = 1 + int(x(i) + half)
            if( m .le. nbkt ) freq(m) = freq(m) + one
 220     continue
         if( lprob ) then
            prb   = exp(-p3)
            Fx(1) = zero
            Fx(2) = prb
            do 230 k = 3, nbkt+1
               prb = prb*p3/float(k-2)
               Fx(k) = Fx(k-1) + prb
 230        continue
         endif
      elseif( idist .eq. 98 ) then		! unknown discrete dist
         m  = int(xmin + half)
         db = one
         do 240 i = 1, npt
            j = 1 + int( x(i) + half ) - m
            if( j .le. nbkt ) freq(j) = freq(j) + one
 240     continue
      else					! unknown continuous dist
         xm  = xmax - xmin
         db  = xm/float(nbkt)
         dbi = one/db
         do 250 i = 1, npt
            m       = 1 + int( dbi*(x(i)-xmin) - de )
            freq(m) = freq(m) + one
 250     continue
      endif

      return
      end
