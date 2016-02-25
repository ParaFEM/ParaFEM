c  *********************************************************************
c  *                                                                   *
c  *                         subroutine dcvit2                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 2.11
c  Written by Gordon A. Fenton, TUNS, July 18, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  computes the initial k1*k2 x k1*k2 covariance matrix (and a 9x9
c           submatrix) between local averages. Used by LAS2G.
c          
c
c  This routine forms the covariance matrix Q between local averages arranged
c  on a k1 x k2 grid using a user provided variance function "dvarfn" that
c  returns the variance of a local average (note that this is the point
c  variance sigma^2 times the traditionally "variance" function, as defined
c  by Vanmarcke in "Random Fields", MIT Press, 1984). The grid and element
c  numbering schemes appear as follows (for k1 = k2 = 5)
c
c          | Dx |
c     ---  --------------------------
c      Dy  | 21 | 22 | 23 | 24 | 25 |
c     ---  --------------------------
c          | 16 | 17 | 18 | 19 | 20 |
c          --------------------------           ----------------
c          | 11 | 12 | 13 | 14 | 15 |           |  7 |  8 |  9 |
c          --------------------------           ----------------
c          |  6 |  7 |  8 |  9 | 10 |           |  4 |  5 |  6 |
c          --------------------------           ----------------
c          |  1 |  2 |  3 |  4 |  5 |           |  1 |  2 |  3 |
c          --------------------------           ----------------
c
c                  Q Array                          R Array
c
c  with the k1 x k2 array numbering on the left and a basic 3 x 3 subarray
c  shown on the right. If we call Z_i the local average of a random process
c  over the i'th cell, then for E[Z] = 0, the elements of Q are defined by
c
c           Q(i,j) = E[Z_i*Z_j]
c
c  Note that the elements of R are simply extracted directly from the
c  appropriate elements of Q (for example, R(1,5) = Q(1,7)) if k1 and k2 are
c  greater than 2. Note also that since both Q and R are symmetric, ONLY THE
c  UPPER TRIANGULAR VALUES are stored herein. Finally note that the random
c  process is assumed to be quadrant symmetric (that is the covariance
c  function is even in each individual lag component) to reduce the total
c  number of covariance calculations.
c
c  Arguments to this routine are as follows;
c
c    dvarfn     external double precision function which returns the
c               variance of a local average of size X x Y. This function
c               is referenced using the call "V = dvarfn(X,Y)" (other
c               parameters to the function must be passed via common blocks).
c               On each invocation, this routine calls DVARFN approximately
c               225 times (via DCVAA2).
c
c    Q          double precision array of size at least k1*k2 x k1*k2 which
c               on output will contain the covariance matrix between the local
c               average elements shown above on the left. Note that since Q is
c               symmetric, only the upper triangular values are actually
c               stored (compatible with DCHOL2). (output)
c
c    iq         leading dimension of the array Q as specified in the calling
c               routine. (input)
c
c    R          double precision array of size at least 9 x 9 which on output
c               will contain the covariance matrix between the local average
c               elements shown above on the right. Note that since R is
c               symmetric, only the upper triangular values are actually
c               stored (compatible with DSIFA). (output)
c
c    ir         leading dimension of the array R as specified in the calling
c               routine. (input)
c
c    k1         number of local average cells in the x direction
c               (horizontally). (input)
c
c    k2         number of local average cells in the y direction
c               (vertically). (input)
c
c    Dx         x-dimension of the local average cells. (input)
c
c    Dy         y-dimension of the local average cells. (input)
c
c  Required:
c   1) from libGAFsim: DCVAA2
c
c  REVISION HISTORY:
c  2.1	eliminated unused local variables (Dec 5/96)
c  2.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine dcvit2( dvarfn, Q, iq, R, ir, k1, k2, Dx, Dy )
      implicit real*8 (a-h,o-z)
      dimension Q(iq,*), R(ir,*)
      integer kx(9), ky(9)
      external dvarfn
      data kx/0,1,2,0,1,2,0,1,2/, ky/0,0,0,1,1,1,2,2,2/

c					first form the essential elements of Q
      ii = 0
      do 20 j = 1, k2
         ty = dble(j-1)
         do 10 i = 1, k1
            ii = ii + 1
            tx = dble(i-1)
            Q(1,ii) = dcvaa2( dvarfn, Dx, Dy, tx, ty )
  10     continue
  20  continue
c					now distribute into the upper triangle
      kk = k1*k2
      do 40 j = 2, kk
         mxj = mod(j-1,k1)
         myj = (j-1)/k1
         do 30 i = 2, j
            mxi    = mod(i-1,k1)
            myi    = (i-1)/k1
            m      = 1 + iabs(mxj - mxi) + k1*iabs(myj - myi)
            Q(i,j) = Q(1,m)
  30     continue
  40  continue
c					form R (9 x 9)
      if( k1 .lt. 3 .or. k2 .lt. 3 ) then
c						find essential elements of R
         ii = 0
         do 60 j = 1, 3
            ty = dble(j-1)
            do 50 i = 1, 3
               ii = ii + 1
               tx = dble(i-1)
               R(1,ii) = dcvaa2( dvarfn, Dx, Dy, tx, ty )
  50        continue
  60     continue
c						and distribute into upper tri.
         do 70 j = 2, 9
         do 70 i = 2, j
            m  = 1 + iabs(kx(j)-kx(i)) + 3*iabs(ky(j)-ky(i))
            R(i,j) = R(1,m)
  70     continue
      else
c						or extract elements of R from Q
         do 80 j = 1, 9
         do 80 i = 1, j
            m      = 1 + iabs(kx(j)-kx(i)) + k1*iabs(ky(j)-ky(i))
            R(i,j) = Q(1,m)
  80     continue
      endif

      return
      end
