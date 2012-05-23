c  *********************************************************************
c  *                                                                   *
c  *                           function chinv                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Dec. 23, 1992
c  Latest Update: Dec 5, 1996
c
c  PURPOSE  returns the inverse Chi-square distribution function for a
c           given probability
c
c  This function employs the inverse gamma distribution dgminv, specialized
c  to produce the Chi-Square distribution results, to compute the inverse
c  Chi-Square distribution function, that is the point x such that
c  p = P[ X <= x], where X follows a Chi-square distribution with n
c  degrees of freedom. Arguments to the function are as follows;
c
c     p    the desired probability level such that p = P[ X <= x] (x is
c          returned). (input)
c
c     n    the number of degrees of freedom associated with the Chi-square
c          distribution. (input)
c
c  ierr    integer error flag which is 0 if all goes well, -1 if convergence
c          not achieved.
c
c  REVISION HISTORY:
c  1.1	eliminated unused variable `one' (Dec 5/96)
c---------------------------------------------------------------------------
      real function chinv( p, n, ierr )
      data two/2.0/

      pa = float(n)/two
      chinv = gminv(p,pa,two,ierr)

      return
      end
