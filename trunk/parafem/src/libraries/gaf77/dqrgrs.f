c  *********************************************************************
c  *                                                                   *
c  *                         subroutine dqrgrs                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, Apr. 9, 1993
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  solve the general least squares problem via QR decomposition
c
c  This program takes a matrix [X] of the m basis functions evaluated at
c  each of the n data points (so that X is n x m in size), and a vector
c  of observed values of y (of length n), and computes the regression
c  coefficients which minimize the sum of squared errors between the
c  observed and predicted values of y. It does this using the QR decomposition
c  algorithm which is generally more accurate than the direct solution of
c  the normal equations via Gaussian Elimination.
c  The prediction equation is assumed to have the form
c
c     y = a_1*X_1(x) + a_2*X_2(x) + ... + a_m*X_m(x)
c
c  where X_i(x) are the basis functions - these are fixed functions of x
c  (which may be a vector itself) and contain no unknown parameters. The
c  matrix [X] is composed of elements
c
c             X_1(x_1)   X_2(x_1) ....  X_m(x_1)
c             X_1(x_2)   X_2(x_2) ....  X_m(x_2)
c             X_1(x_3)   X_2(x_3) ....  X_m(x_3)
c                .          .     .        .
c                .          .      .       .
c                .          .       .      .
c             X_1(x_n)   X_2(x_n) ....  X_m(x_n)
c
c  where x_1, x_2, ... are the 1st, 2nd,... data point locations. Note
c  that m must be less than or equal to n.
c  Arguments to this routine are as follows;
c
c     X   real array of size at least n x m which contains the elements of
c         the basis functions evaluated at each data location (see above).
c         Note that the contents of X are destroyed by this routine.
c         (input/destroyed)
c
c    ix   leading dimension of the array X as specified in the calling routine.
c         (input)
c
c     y   real vector of length at least n containing the observations. On
c         output, y will contain the elements of {a} in its first m locations.
c         (input/output)
c
c     n   number of data points. (input)
c
c     m   number of basis functions. (input)
c
c    se   residual sum of squared errors. (output)
c
c    r2   r-squared value denoting the proportion of the total error
c         accounted for by the regression model (ie the coefficient of
c         determination). (output)
c
c  ierr   integer flag which is assigned the following possible values (output)
c          =  0 if all goes well
c          = -1 if m > n so that no solution can be found
c          = -2 if the basis functions are not linearly independent on the
c               set of data (the system of equations is singular).
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      subroutine dqrgrs( X, ix, y, n, m, se, r2, ierr )
      implicit real*8 (a-h,o-z)
      dimension X(ix,*), y(*)
      data zero/0.d0/
      float(i) = dble(i)
c					check m versus n
      if( m .gt. n ) then
         ierr = -1
         return
      endif
c					compute SSTO
      sy  = zero
      syy = zero
      do 10 i = 1, n
         sy  = sy + y(i)
         syy = syy + y(i)*y(i)
  10  continue
      ssto = syy - sy*sy/float(n)
c					compute the QR decomposition of X
      call dqrhde( X, ix, n, m )
c					compute the least square regression
      se = dqrhsq( X, ix, y, n, m )
      if( se .lt. zero ) then
         ierr = -2
         return
      endif
c					compute r-squared
      r2 = (ssto - se)/ssto
c					all done
      ierr = 0
      return
      end
