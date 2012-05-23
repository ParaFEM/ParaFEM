c  *********************************************************************
c  *                                                                   *
c  *                       subroutine statms                           *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.7
c  Written by Gordon A. Fenton, TUNS, 1992
c  Latest Update: Jul 30, 2006
c
c  PURPOSE  estimates distribution parameters for a set of data by method
c           of moments and finds its range
c
c  Arguments to this routine are as follows;
c
c	x	real vector of length at least npt containing the raw data.
c		(input)
c
c     npt	integer giving the number of elements in x. (input)
c
c   idist	integer type code specifying the type of distribution
c		x is assumed to follow;
c		   = 0	for a normal distribution
c		   = 1	for a lognormal distribution
c		   = 2	for an exponential distribution
c		   = 3	for a Beta distribution (-3 if bounds p3, p4 are
c			prescribed on input)
c		   = 4  for a Gamma distribution
c		   = 5  for a uniform distribution
c		   = 6  for a binomial (n,p) distribution (discrete)
c		   = 7  for a geometric (p) distribution (discrete)
c		   = 8  for a negative binomial (m,p) distribution (discrete)
c		   = 9  for a Poisson (lambda*t) distribution (discrete)
c		   = 98 if distribution is discrete but unknown
c		   = 99 if distribution is continuous but unknown
c		(input)
c
c	u	estimated mean of the raw data. (output)
c
c	s	estimated standard deviation of the raw data. (output)
c
c      p1	the first parameter of the assumed distribution, p1 takes on
c		the following meanings depending on idist; (output)
c		  if idist = 0,  p1 is the mean of the data
c		  if idist = 1,  p1 is the mean of the log-data
c		  if idist = 2,  p1 is the lambda parameter of the exp. dist.
c		  if idist = 3,  p1 is the alpha of the Beta [0,1] dist.
c		  if idist = 4,  p1 is the alpha of the Gamma dist.
c		  if idist = 5,  p1 is the expected lower bound of the data
c		  if idist = 6,  p1 is the prob of success, p, of the binomial
c		  if idist = 7,  p1 is the prob of success, p, of the geometric
c		  if idist = 8,  p1 is the prob of success, p, of the neg binom
c		  if idist = 9,  p1 is the mean lambda*t of the Poisson dist
c		  if idist = 98, p1 is undefined
c		  if idist = 99, p1 is undefined
c		Notes:
c		  1) for idist = 1, the lognormal distribution, p1 is
c		     calculated as the average of the log-data (not by
c		     transforming u and s). In general, unless the data
c		     exactly matches a lognormal distribution, the two
c		     sets of parameters (p1,p2) and (u,s) will not precisely
c		     transform to each other under the lognormal
c		     transformations.
c		  2) in the case of the Beta distribution, idist = 3, the
c		     distribution is assumed to be bounded on the interval
c		     [p3,p4], while the parameters alpha and beta which are
c		     actually returned are those for the transformed interval
c		     [0,1]. If the bounds [p3,p4] are unknown, then the
c		     uniform distribution bounds (see next) are used.
c		  3) for the uniform distribution, idist = 5, the bounds
c		     are the `expected' bounds, that is
c			p1 =  xmin - (xmax-xmin)/(n-2)
c			p2 =  xmax + (xmax-xmin)/(n-2)
c		
c      p2	the second parameter of the assumed distribution, p2 takes on
c		the following meanings depending on idist; (output)
c		  if idist = 0,  p2 is the s.d. of the data
c		  if idist = 1,  p2 is the s.d. of the log-data
c		  if idist = 2,  p2 is the s.d. of the data
c		  if idist = 3,  p2 is the beta par of the Beta [0,1] dist
c		  if idist = 4,  p2 is the beta parameter of the Gamma dist
c		  if idist = 5,  p2 is the expected upper bound of the data
c		  if idist = 6,  p2 is the s.d. of the data
c		  if idist = 7,  p2 is the s.d. of the data
c		  if idist = 8,  p2 is the s.d. of the data
c		  if idist = 9,  p2 is the s.d. of the data
c		  if idist = 98, p2 is undefined
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
c		(input/output)
c
c      p4	the fourth parameter of the assumed distribution, at this time
c		defined only for idist = 3 (Beta distribution): p4 is
c		estimated from the data to be the `expected' upper bound,
c		p4 = xmax + (xmax - xmin)/(n-2). Alternatively, if
c		idist = -3, then p4 is assumed to contain the prescribed
c		distribution upper bound on input and it is not changed.
c		(input/output)
c
c    xmin	minimum of the data. (output)
c
c    xmax	maximum of the data. (output)
c
c    ierr	integer flag which is zero if all goes well, -i if the
c		i'th element does not fit the assumed distribution, in which
c		case the distribution parameters p1 through p4 are not
c		defined here. (output)
c
c  REVISION HISTORY:
c  1.2	added ability to specify a range of distributions.
c  1.3	separated estimation of mean and standard deviation of the raw data
c	from the estimation of the distribution parameters.
c  1.31	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c  1.4	now compute average separately from sample variance to avoid
c	catastrophic cancellation (Mar 7/00)
c  1.5	set unused parameters to zero (May 26/01)
c  1.6	added discrete distributions binomial, geometric, neg binom, and
c	Poisson (Jun 20/05)
c  1.7	add error test to ensure non-zero values for lognormal, geometric,
c	and negative binomial. (Jul 30/06)
c----------------------------------------------------------------------------
      subroutine statms( x,npt,idist,u,s,p1,p2,p3,p4,xmin,xmax,ierr )
      dimension x(*)
      real*8 d, dx, dr, dq, sm, sv, dble, dsqrt, dlog
      logical lbzero, lgzero
      data zero/0.0/, half/0.5/, one/1.0/

c					set defaults
      p1 = zero
      p2 = zero
c					check sample size
      if( npt .le. 0 ) then
         ierr = 2
         return
      elseif( npt .eq. 1 ) then
         u = x(1)
         s = zero
         ierr = 1
         return
      else
         ierr = 0
      endif
c					these dists are bounded below by 0

      lgzero = (idist .eq. 1) .or. (idist .eq. 7) .or. (idist .eq. 8)
      lbzero = (idist .eq. 1) .or. (idist .eq. 2) .or. (idist .eq. 4)
     >    .or. ((idist .ge. 6) .and. (idist .le. 9))

c					compute min, max, u, s
      if( lbzero .and. x(1) .lt. zero ) ierr = -1
      if( lgzero .and. x(1) .eq. zero ) ierr = -1
      xmin = x(1)
      xmax = x(1)
      sm   = dble(x(1))				! estimate the mean
      do 10 i = 2, npt
         if( ierr .eq. 0 ) then			! only record first error
            if( lbzero .and. x(i) .lt. zero ) ierr = -i
            if( lgzero .and. x(i) .eq. zero ) ierr = -i
         endif
         d  = dble( x(i) )
         sm = sm + d
         if( x(i) .lt. xmin ) xmin = x(i)
         if( x(i) .gt. xmax ) xmax = x(i)
  10  continue
      d  = sm/dble(npt)			! the sample average
      u  = sngl(d)
      dr = dble(x(1)) - d			! estimate the variance
      sv = dr*dr
      do 20 i = 2, npt
         dr = dble(x(i)) - d
         sv = sv + dr*dr
  20  continue
      sv  = sv/dble(npt-1)		! the sample variance
      s   = sngl( dsqrt(sv) )
      if( ierr .lt. 0 ) return
c					compute distribution parameters
      if( idist .eq. 0 ) then		! normal distribution
         p1 = u
         p2 = s
      elseif( idist .eq. 1 ) then	! lognormal distribution
         d  = dble(x(1))			! estimate the mean
         sm = dlog(d)
         do 30 i = 2, npt
            d  = dble( x(i) )
            sm = sm + dlog(d)
  30     continue
         d  = sm/dble(npt)		! the sample average of the log-data
         p1 = sngl(d)
         dq = dble(x(1))			! estimate the variance
         dr = dlog(dq) - d
         sv = dr*dr
         do 40 i = 2, npt
            dq = dble(x(i))
            dr = dlog(dq) - d
            sv = sv + dr*dr
  40     continue
         sv = sv/dble(npt-1)
         p2 = sngl( dsqrt(sv) )
      elseif( idist .eq. 2 ) then	! exponential distribution
         p1 = one/u
      elseif( iabs(idist) .eq. 3 ) then	! Beta distribution
         if( idist .eq. 3 ) then
            if( npt .gt. 2 ) then
               p3 = xmin - (xmax - xmin)/float(npt-2)
               p4 = xmax + (xmax - xmin)/float(npt-2)
            else
               p3 = xmin
               p4 = xmax
            endif
         endif
         dx = dble(p4) - dble(p3)		! need to transform to [0,1]
         sm = (d - dble(p3))/dx
         sv = sv/(dx*dx)
         dr = (1.d0 - sm)/sm
         dq = 1.d0 + dr
         sm = ((dr/(sv*dq*dq)) - 1.d0)/dq
         sv = sm*dr
         p1 = sngl(sm)
         p2 = sngl(sv)
      elseif( idist .eq. 4 ) then		! Gamma distribution
         p1 = sngl(d*d/sv)
         p2 = sngl(sv/d)
      elseif( idist .eq. 5 ) then		! uniform distribution
         if( npt .gt. 2 ) then
            p1 = xmin - (xmax - xmin)/float(npt-2)
            p2 = xmax + (xmax - xmin)/float(npt-2)
         else
            p1 = xmin
            p2 = xmax
         endif
      elseif( idist .eq. 6 ) then		! binomial (n,p)
         p1 = u/p3				! p3 contains n
      elseif( idist .eq. 7 ) then
         p1 = one/u
      elseif( idist .eq. 8 ) then
         p1 = p3/u
      elseif( idist .eq. 9 ) then
         p1 = u
      endif

      return
      end
