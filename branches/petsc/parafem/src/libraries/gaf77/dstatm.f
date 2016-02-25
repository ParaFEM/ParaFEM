c  *********************************************************************
c  *                                                                   *
c  *                       subroutine statms                           *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.4
c  Written by Gordon A. Fenton, TUNS, 1992
c  Latest Update: Mar 7, 2000
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
c		   = 99 if distribution is unknown
c		(input)
c
c	u	estimated mean of the raw data. (output)
c
c	s	estimated standard deviation of the raw data. (output)
c
c      p1	the first parameter of the assumed distribution, p1 takes on
c		the following meanings depending on idist; (output)
c		  if idist = 0, p1 is the mean of the data
c		  if idist = 1, p1 is the mean of the log-data
c		  if idist = 2, p1 is the lambda parameter of the exp. dist.
c		  if idist = 3, p1 is the alpha of the Beta [0,1] dist.
c		  if idist = 4, p1 is the alpha of the Gamma dist.
c		  if idist = 5, p1 is the expected lower bound of the data
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
c		     [0,1].
c		  3) for the uniform distribution, idist = 5, the bounds
c		     are the `expected' bounds, that is
c			p1 =  xmin - (xmax-xmin)/(n-2)
c			p2 =  xmax + (xmax-xmin)/(n-2)
c		
c      p2	the second parameter of the assumed distribution, p2 takes on
c		the following meanings depending on idist; (output)
c		  if idist = 0, p2 is the standard deviation (s.d.) of the data
c		  if idist = 1, p2 is the s.d. of the log-data
c		  if idist = 2, p2 is the s.d. of the data
c		  if idist = 3, p2 is the beta parameter of the Beta [0,1] dist
c		  if idist = 4, p2 is the beta parameter of the Gamma dist
c		  if idist = 5, p2 is the expected upper bound of the data
c		  if idist = 99, p2 is undefined
c
c      p3	the third parameter of the assumed distribution, at this time
c		defined only for idist = 3, the Beta distribution: p3 is the
c		`expected' lower bound, p3 = xmin - (xmax - xmin)/(n-2). If
c		idist = -3, then p3 is the prescribed distribution lower
c		bound on input and it is not changed. Otherwise p3 is estimated
c		from the data. (input/output)
c
c      p4	the fourth parameter of the assumed distribution, at this time
c		defined only for idist = 3, the Beta distribution: p4 is the
c		`expected' upper bound, p4 = xmax + (xmax - xmin)/(n-2). If
c		idist = -3, then p4 is the prescribed distribution upper
c		bound on input and it is not changed. Otherwise p4 is estimated
c		from the data. (input/output)
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
c----------------------------------------------------------------------------
      subroutine dstatm( x,npt,idist,u,s,p1,p2,p3,p4,xmin,xmax,ierr )
      implicit real*8 (a-h,o-z)
      dimension x(*)
      logical lbzero
      data zero/0.d0/, one/1.d0/

      if( npt .le. 0 ) then
         ierr = -(npt + 1)
         return
      elseif( npt .eq. 1 ) then
         u = x(1)
         ierr = 1
         return
      else
         ierr = 0
      endif
c					these dists are bounded below by 0

      lbzero = (idist .eq. 1) .or. (idist .eq. 2) .or. (idist .eq. 4)

c					compute min, max, u, s
      if( lbzero .and. x(1) .lt. zero ) ierr = -1
      xmin = x(1)
      xmax = x(1)
      sm   = x(1)				! estimate the mean
      do 10 i = 2, npt
         if( lbzero .and. x(i) .lt. zero ) ierr = -i
         sm = sm + x(i)
         if( x(i) .lt. xmin ) xmin = x(i)
         if( x(i) .gt. xmax ) xmax = x(i)
  10  continue
      u  = sm/dble(npt)			! the sample average
      dr = x(1) - u				! estimate the variance
      sv = dr*dr
      do 20 i = 2, npt
         dr = x(i) - u
         sv = sv + dr*dr
  20  continue
      sv = sv/dble(npt-1)		! the sample variance
      s  = dsqrt(sv)
      if( ierr .lt. 0 ) return
c					compute distribution parameters
      if( idist .eq. 0 ) then		! normal distribution
         p1 = u
         p2 = s
      elseif( idist .eq. 1 ) then	! lognormal distribution
         sm = dlog(x(1))			! estimate the mean
         do 30 i = 2, npt
            sm = sm + dlog(x(i))
  30     continue
         p1 = sm/dble(npt)		! the sample average of the log-data
         dr = dlog(x(1)) - p1
         sv = dr*dr
         do 40 i = 2, npt
            dr = dlog(x(i)) - p1
            sv = sv + dr*dr
  40     continue
         sv = sv/dble(npt-1)
         p2 = dsqrt(sv)
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
         dx = p4 - p3			! need to transform to [0,1]
         sm = (u - p3)/dx
         sv = sv/(dx*dx)
         dr = (one - sm)/sm
         dq = one + dr
         sm = ((dr/(sv*dq*dq)) - one)/dq
         sv = sm*dr
         p1 = sm
         p2 = sv
      elseif( idist .eq. 4 ) then		! Gamma distribution
         p1 = u*u/sv
         p2 = sv/u
      elseif( idist .eq. 5 ) then		! uniform distribution
         if( npt .gt. 2 ) then
            p1 = xmin - (xmax - xmin)/float(npt-2)
            p2 = xmax + (xmax - xmin)/float(npt-2)
         else
            p1 = xmin
            p2 = xmax
         endif
      endif

      return
      end
