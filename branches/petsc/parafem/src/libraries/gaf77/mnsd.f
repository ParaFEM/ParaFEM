c  *********************************************************************
c  *                                                                   *
c  *                          subroutine mnsd                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, Dalhousie University, Apr 26, 2006
c  Latest Update: Apr 26, 2006
c
c  PURPOSE  computes the sample mean and standard deviation of a data set
c           along with the standard deviations of these estimates
c
c  DESCRIPTION
c  This routine computes the sample mean,
c
c            1   n
c	a = --- sum x(i)
c            n  i=1
c
c  the sample standard deviation
c
c                    1    n              2
c       s = sqrt(  ----- sum [ x(i) - a ]  )
c                   n-1  i=1
c
c  the estimated standard deviation of the sample mean
c
c               s
c       sa = -------
c            sqrt(n)
c
c  and the estimated standard deviation of the sample variance (s^2)
c
c                    2      2
c       sv = sqrt( ----- ) s
c                   n-1
c
c  ARGUMENTS
c	x	real vector of length at least n containing the sample
c		of the population whose mean and standard deviation are
c		to be estimated. (input)
c
c	n	integer containing the size of the sample contained in
c		the vector x. (input)
c
c	a	the sample average. (output)
c
c	s	the sample standard deviation. If n < 2, this is set to
c		zero. (output)
c
c      sa	the estimated standard deviation of the sample average.
c		If n < 2, this is set to zero. (output)
c
c      sv	the estimated standard deviation of the sample variance.
c		(Sorry, at the moment I don't know the estimated standard
c		deviation of the sample standard deviation) If n < 2, this
c		is set to zero. (output)
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      subroutine mnsd(x,n,a,s,sa,sv)
      real*4 x(*)
c				sample mean (average)
      a = 0.0
      do 10 i = 1, n
         a = a + x(i)
  10  continue
      a = a/float(n)
c				sample standard deviation
      s  = 0.0			! by default this is zero if n < 2
      sa = 0.0			! by default this is zero if n < 2
      sv = 0.0			! by default this is zero if n < 2
      if( n .gt. 1 ) then
         v = 0.0
         do 20 i = 1, n
            d = x(i) - a
            v = v + d*d
  20     continue
         v = v/float(n-1)
         s = sqrt(v)
c				standard devation of the average
         sa = s/sqrt(float(n))
c				standard deviation of v
         sv = v*sqrt(2.0/float(n-1))
      endif
c				all done, close up shop
      return
      end
