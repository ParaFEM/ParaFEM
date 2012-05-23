c  *********************************************************************
c  *                                                                   *
c  *                        subroutine dbest1                          *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, Dalhousie University, June 16, 2006
c  Latest Update: Jun 16, 2006
c
c  PURPOSE  estimates intercept and slope by least squares of the straight
c           line y = a_1 + a_2*x given a set of data
c
c  This routine takes as input a set of n observations of (x,y) and fits
c  a straight line to the data via least squares regression. The intercept
c  and slope values are returned. Arguments to this routine are as follows;
c
c       x	real*8 vector of length at least n containing the observed
c		x values. (input)
c
c	y	real*8 vector of length at least n containing the observed
c		y values. (input)
c
c     tmp	temporary real*8 vector of length at least 3*n used as
c		workspace.
c
c	n	the number of observations. (input)
c
c      a1	the best fit intercept. (output)
c
c      a2	the best fit slope. (output)
c
c      r2	r-squared value denoting the proportion of the total error
c		accounted for by the regression model (ie the coefficient of
c		determination). (output)
c
c    ierr	integer error flag. (output)
c		  =  0 if all goes well
c		  = -1 if n < 2 so that a line is not defined
c		  = -2 if the line cannot be found (as in a vertical line)
c-------------------------------------------------------------------------
      subroutine dbest1(x,y,tmp,n,a1,a2,r2,ierr)
      implicit real*8 (a-h,o-z)
      dimension x(*), y(*), tmp(*)
      data one/1.d0/
c				transfer x and y into tmp
      do 10 i = 1, n
         tmp(i) = one
         tmp(n+i) = x(i)
         tmp(2*n+i) = y(i)
  10  continue
c				call qrhsq to solve the least squares problem

      call dqrgrs( tmp, n, tmp(2*n+1), n, 2, se, r2, ierr )

      a1 = tmp(2*n+1)
      a2 = tmp(2*n+2)

      return
      end

