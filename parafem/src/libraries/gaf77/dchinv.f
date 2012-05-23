c  *********************************************************************
c  *                                                                   *
c  *                          function dchinv                          *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Dec. 23, 1992
c  Latest Update: Dec 4, 1998
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
c  1.1	eliminated unused local variable `one' (Dec 5/96)
c  1.2	explicitly call dgminv to get ierr returned (Dec 4/98)
c---------------------------------------------------------------------------
      real*8 function dchinv( p, n, ierr )
      implicit real*8 (a-h,o-z)
      data two/2.d0/
c	NOTE: Linux's f77 compiler doesn't pass args in the following by
c	      address... thus `i' can't be returned in `ierr'...
c	      we have to specifically call dbetnv below.
c     gminv(z1,z2,z3,i) = dgminv(z1,z2,z3,i)
      float(i) = dble(i)

      pa = float(n)/two
      dchinv = dgminv(p,pa,two,ierr)

      return
      end
