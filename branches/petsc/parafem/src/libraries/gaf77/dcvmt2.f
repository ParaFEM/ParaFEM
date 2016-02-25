c  *********************************************************************
c  *                                                                   *
c  *                       subroutine dcvmt2                           *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.11
c  Written by Gordon A. Fenton, TUNS, July 17, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   computes the covariance matrix between local averages needed
c            by LAS2G.
c
c  This routine computes the three covariance matrices, R (9x9), B (4x4),
c  and S (9x4), required by LAS2G. The matrices R and B represent the
c  covariances between a 3 x 3 array of equal sized local average elements
c  and a 2 x 2 subarray respectively, numbered as follows;
c
c              |-- Dx--|   R Array                         B Array
c        ---   -------------------------           -------------------------
c         |    |       |       |       |           |       |       |       |
c         Dy   |   7   |   8   |   9   |           |       |       |       |
c         |    |       |       |       |           |       |       |       |
c        ---   -------------------------           -------------------------
c              |       |       |       |           |       |       |       |
c              |   4   |   5   |   6   |           |   3   |   4   |       |
c              |       |       |       |           |       |       |       |
c              -------------------------           -------------------------
c              |       |       |       |           |       |       |       |
c              |   1   |   2   |   3   |           |   1   |   2   |       |
c              |       |       |       |           |       |       |       |
c              -------------------------           -------------------------
c
c                                        Figure 1
c
c  Notice that the array B is just equal to selected elements of the array R
c  (ie. B(1,1) = R(1,1), B(1,4) = R(1,5), etc.). S represents the 9x4
c  covariance matrix between a doubly large array and an interior subdivided
c  2 x 2 array as shown, where the numbering of the larger array is the same
c  as on the left in Figure 1.
c
c              |--2Dx--|
c        ---   -------------------------
c         |    |       |       |       |
c        2Dy   |       |       |       |
c         |    |       |       |       |
c        ---   -------------------------
c              |       | 3 | 4 |       |
c              |       |-------|       |                     Figure 2
c              |       | 1 | 2 |       |
c              -------------------------
c              |       |       |       |
c              |       |       |       |
c              |       |       |       |
c              -------------------------
c
c  Finally note that the random process is assumed to be quadrant symmetric
c  (that is the covariance function is even in each individual lag component)
c  to reduce the total number of covariance calculations.
c
c  Arguments to this routine are as follows;
c
c    dvarfn     external double precision function which returns the
c               variance of a local average of size X x Y. This function
c               is referenced using the call "V = dvarfn(X,Y)" (other
c               parameters to the function must be passed via common blocks).
c               On each invocation, this routine calls DVARFN approximately
c               225 times (via DCVAA2 and DCVAB2).
c
c    R          double precision array of size at least 9 x 9 which on output
c               will contain the covariance matrix between the local average
c               elements shown above in Figure 1, left. Note that since R is
c               symmetric, only the upper triangular values are actually
c               stored (this is to be compatible with DSIFA). (output)
c
c    ir         leading dimension of the array R as specified in the calling
c               routine. (input)
c
c    B          double precision array of size at least 4 x 4 which on output
c               will contain the covariance matrix between the local average
c               elements shown above in Figure 1, right. Note that since B is
c               symmetric, only the upper triangular values are actually
c               stored (to be compatible with DCHOL2). (output)
c
c    ib         leading dimension of the array B as specified in the calling
c               routine. (input)
c
c    S          double precision array of size at least 9 x 4 which on output
c               will contain the covariances between the set of subdivided
c               central elements and the parent 3 x 3 set of doubly large
c               elements (see Figure 2). (output)
c
c    is         leading dimension of the array S as specified in the calling
c               routine. (input)
c
c    Dx         x-dimension of the subdivided elements. The dimension of the
c               parent elements is assumed to be 2*Dx. (input)
c
c    Dy         y-dimension of the subdivided elements. The dimension of the
c               parent elements is assumed to be 2*Dy. (input)
c
c    lformR     logical flag which is true if the matrix R is to be formed
c               (in some cases, only the smaller matrix B is required). (input)
c
c  REVISION HISTORY:
c  1.1	eliminated unused local variable `m' (Dec 5/96)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine dcvmt2( dvarfn, R, ir, B, ib, S, is, Dx, Dy, lformR )
      implicit real*8 (a-h,o-z)
      dimension R(ir,*), B(ib,*), S(is,*)
      logical lformR
      external dvarfn
      data zero/0.d0/, half/0.5d0/, one/1.d0/, onept5/1.5d0/
      data two/2.d0/, twopt5/2.5d0/
c					first find B (upper triangle only)
      B(1,1) = dcvaa2(dvarfn,Dx,Dy,zero,zero)
      B(1,2) = dcvaa2(dvarfn,Dx,Dy,one,zero)
      B(1,3) = dcvaa2(dvarfn,Dx,Dy,zero,one)
      B(1,4) = dcvaa2(dvarfn,Dx,Dy,one,one)
      B(2,2) = B(1,1)
      B(2,3) = B(1,4)
      B(2,4) = B(1,3)
      B(3,3) = B(1,1)
      B(3,4) = B(1,2)
      B(4,4) = B(1,1)
c					now find R optionally (upper triangle)
      if( lformR ) then
         R(1,1) = B(1,1)
         R(2,2) = B(1,1)
         R(3,3) = B(1,1)
         R(4,4) = B(1,1)
         R(5,5) = B(1,1)
         R(6,6) = B(1,1)
         R(7,7) = B(1,1)
         R(8,8) = B(1,1)
         R(9,9) = B(1,1)

         R(1,2) = B(1,2)
         R(2,3) = B(1,2)
         R(4,5) = B(1,2)
         R(5,6) = B(1,2)
         R(7,8) = B(1,2)
         R(8,9) = B(1,2)

         R(1,3) = dcvaa2(dvarfn,Dx,Dy,two,zero)
         R(4,6) = R(1,3)
         R(7,9) = R(1,3)

         R(1,4) = B(1,3)
         R(2,5) = B(1,3)
         R(3,6) = B(1,3)
         R(4,7) = B(1,3)
         R(5,8) = B(1,3)
         R(6,9) = B(1,3)

         R(1,5) = B(1,4)
         R(2,6) = B(1,4)
         R(4,8) = B(1,4)
         R(5,9) = B(1,4)
         R(2,4) = B(1,4)
         R(3,5) = B(1,4)
         R(5,7) = B(1,4)
         R(6,8) = B(1,4)

         R(1,6) = dcvaa2(dvarfn,Dx,Dy,two,one)
         R(4,9) = R(1,6)
         R(3,4) = R(1,6)
         R(6,7) = R(1,6)
         R(1,7) = dcvaa2(dvarfn,Dx,Dy,zero,two)
         R(2,8) = R(1,7)
         R(3,9) = R(1,7)
         R(1,8) = dcvaa2(dvarfn,Dx,Dy,one,two)
         R(2,9) = R(1,8)
         R(2,7) = R(1,8)
         R(3,8) = R(1,8)
         R(1,9) = dcvaa2(dvarfn,Dx,Dy,two,two)
         R(3,7) = R(1,9)
      endif
c				now find S
      S(1,1) = dcvab2(dvarfn,Dx,Dy, onept5, onept5 )
      S(2,1) = dcvab2(dvarfn,Dx,Dy, half,   onept5 )
      S(3,1) = dcvab2(dvarfn,Dx,Dy, twopt5, onept5 )
      S(4,1) = dcvab2(dvarfn,Dx,Dy, onept5, half   )
      S(5,1) = dcvab2(dvarfn,Dx,Dy, half,   half   )
      S(6,1) = dcvab2(dvarfn,Dx,Dy, twopt5, half   )
      S(7,1) = dcvab2(dvarfn,Dx,Dy, onept5, twopt5 )
      S(8,1) = dcvab2(dvarfn,Dx,Dy, half,   twopt5 )
      S(9,1) = dcvab2(dvarfn,Dx,Dy, twopt5, twopt5 )
      S(1,2) = S(3,1)
      S(2,2) = S(2,1)
      S(3,2) = S(1,1)
      S(4,2) = S(6,1)
      S(5,2) = S(5,1)
      S(6,2) = S(4,1)
      S(7,2) = S(9,1)
      S(8,2) = S(8,1)
      S(9,2) = S(7,1)
      S(1,3) = S(7,1)
      S(2,3) = S(8,1)
      S(3,3) = S(9,1)
      S(4,3) = S(4,1)
      S(5,3) = S(5,1)
      S(6,3) = S(6,1)
      S(7,3) = S(1,1)
      S(8,3) = S(2,1)
      S(9,3) = S(3,1)
      S(1,4) = S(9,1)
      S(2,4) = S(8,1)
      S(3,4) = S(7,1)
      S(4,4) = S(6,1)
      S(5,4) = S(5,1)
      S(6,4) = S(4,1)
      S(7,4) = S(3,1)
      S(8,4) = S(2,1)
      S(9,4) = S(1,1)

      return
      end
