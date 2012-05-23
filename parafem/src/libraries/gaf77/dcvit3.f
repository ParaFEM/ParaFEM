c  *********************************************************************
c  *                                                                   *
c  *                         subroutine dcvit3                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 2.0
c  Written by Gordon A. Fenton, TUNS, Aug. 2, 1993
c
c  PURPOSE  computes the initial k1*k2*k3 x k1*k2*k3 covariance matrix (and a
c           27x27 submatrix) between local averages. Used by LAS3G.
c          
c
c  This routine forms the covariance matrix R0 between local averages arranged
c  on a k1 x k2 x k3 grid using a user provided variance function "dvarfn" that
c  returns the variance of a local average (note that this is the point
c  variance sigma^2 times the traditionally "variance" function, as defined
c  by Vanmarcke in "Random Fields", MIT Press, 1984). The grid and element
c  numbering schemes appear as follows (for k1 = k2 = 5, k3 = 1 plane shown)
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
c                  R0 Array                          R Array
c
c  with the k1 x k2 array numbering on the left and a basic 3 x 3 subarray
c  shown on the right. In the depth direction (not shown), the numbering
c  continues; ie the element immediately behind element 1 on the left is
c  element number 26, etc. If we call Z_i the local average of a random process
c  over the i'th cell, then for E[Z] = 0, the elements of R0 are defined by
c
c           R0(i,j) = E[Z_i*Z_j]
c
c  Note that the elements of R are simply extracted directly from the
c  appropriate elements of R0 (for example, R(1,5) = R0(1,7)) if k1, k2 and k3
c  are greater than 2. Note also that since both R0 and R are symmetric, ONLY
c  THE UPPER TRIANGULAR VALUES are stored herein. Finally note that
c  the random process is assumed to be quadrant symmetric (that is the
c  covariance function is even in each individual lag component) to reduce
c  the total number of covariance calculations.
c
c  Arguments to this routine are as follows;
c
c    dvarfn     external double precision function which returns the
c               variance of a local average of size X x Y x Z. This function
c               is referenced using the call "V = dvarfn(X,Y,Z)" (other
c               parameters to the function must be passed via common blocks).
c
c    R0          double precision array of size at least k1*k2*k3 x k1*k2*k3
c               which on output will contain the covariance matrix between
c               the local average elements shown above on the left. Note
c               that since R0 is symmetric, only the upper triangular values
c               are actually stored (compatible with DCHOL2). (output)
c
c    iq         leading dimension of the array R0 as specified in the calling
c               routine. (input)
c
c    R          double precision array of size at least 27 x 27 which on output
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
c    k3         number of local average cells in the z direction
c               (depth). (input)
c
c    Dx         x-dimension of the local average cells. (input)
c
c    Dy         y-dimension of the local average cells. (input)
c
c    Dz         z-dimension of the local average cells. (input)
c
c  Required:
c   1) from libGAFsim:	DCVAA3
c----------------------------------------------------------------------------
      subroutine dcvit3( dvarfn, R0, iq, R, ir, k1, k2, k3, Dx, Dy, Dz )
      implicit real*8 (a-h,o-z)
      dimension R0(iq,*), R(ir,*)
      external dvarfn
c					first form the essential elements of R0
      ii = 0
      do 30 k = 1, k3
         tz = dble(k-1)
         do 20 j = 1, k2
            ty = dble(j-1)
            do 10 i = 1, k1
               ii = ii + 1
               tx = dble(i-1)
               R0(1,ii) = dcvaa3( dvarfn, Dx, Dy, Dz, tx, ty, tz )
  10        continue
  20     continue
  30  continue
c					now distribute into the upper triangle
      kk  = k1*k2*k3
      k12 = k1*k2
      do 50 n = 2, kk
         n3 = (n-1)/k12
         n2 = (n - k12*n3 - 1)/k1
         do 40 m = 2, n
            m3 = (m-1)/k12
            m2 = (m - k12*m3 - 1)/k1
            i3 = m3 - n3
            i2 = m2 - n2
            i1 = (m-n) - k1*i2 - k12*i3
            k  = 1 + iabs(i1) + k1*iabs(i2) + k12*iabs(i3)
            R0(m,n) = R0(1,k)
  40     continue
  50  continue
c					form R (27 x 27)

      if( k1 .lt. 3 .or. k2 .lt. 3 .or. k3 .lt. 3 ) then
c						find essential elements of R
         ii = 0
         do 80 k = 1, 3
            tz = dble(k-1)
            do 70 j = 1, 3
               ty = dble(j-1)
               do 60 i = 1, 3
                  tx = dble(i-1)
                  ii = ii + 1
                  R(1,ii) = dcvaa3( dvarfn, Dx, Dy, Dz, tx, ty, tz )
  60           continue
  70        continue
  80     continue
c						and distribute into upper tri.
         do 100 n = 2, 27
            n3 = (n-1)/9
            n2 = (n - 9*n3 - 1)/3
            do 90 m = 2, n
               m3 = (m-1)/9
               m2 = (m - 9*m3 - 1)/3
               i3 = m3 - n3
               i2 = m2 - n2
               i1 = (m-n) - 3*i2 - 9*i3
               k  = 1 + iabs(i1) + 3*iabs(i2) + 9*iabs(i3)
               R(m,n) = R(1,k)
  90        continue
 100     continue
      else
c						or extract R from R0
         do 120 n = 1, 27
            n3 = (n-1)/9
            n2 = (n - 9*n3 - 1)/3
            do 110 m = 1, n
               m3 = (m-1)/9
               m2 = (m - 9*m3 - 1)/3
               i3 = m3 - n3
               i2 = m2 - n2
               i1 = (m-n) - 3*i2 - 9*i3
               k  = 1 + iabs(i1) + k1*iabs(i2) + k12*iabs(i3)
               R(m,n) = R0(1,k)
 110        continue
 120     continue
      endif

      return
      end

