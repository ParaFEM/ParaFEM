c  ************************************************************************
c  *                                                                      *
c  *                        Subroutine dorthH                             *
c  *                                                                      *
c  ************************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, Dec. 1990
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  creates an orthogonal transformation matrix from a given Nth row
c
c  This routine creates an orthogonal transformation matrix of size N x N
c  whose Nth row is the provided unit vector D. The algorithm proceeds by
c  successive rotations of the identity matrix; first in the (x1,x2) plane,
c  then the (x2,x3) plane, etc. until the nth basis vector is aligned with
c  the vector D. In each rotation, the (i+1)th basis vector is aligned with
c  the projection of D onto the (x_i,x_i+1) plane. Arguments to the routine
c  are as follows;
c
c    H   real array of size at least N x N which on output will contain
c        the orthogonal transformation matrix with D as its Nth row. (output)
c
c   ih   integer containing the leading dimension of H exactly as specified in
c        the calling routine. (input)
c
c    D   real vector of length at least N which on input contains the
c        unit vector to which the N'th basis vector will be aligned.
c        On output, D should be composed of the N elements {0,0,...,1}.
c        The degree to which numerically produced D differs from the
c        above is an indication of the round-off errors incurred in the
c        routine. (input/output)
c
c    N   integer giving the length of the vector D and the size of H. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine dorthH( H, ih, D, N )
      implicit real*8 (a-h,o-z)
      dimension H(ih,*), D(*)
      data zero/0.d0/, one/1.d0/
      sqrt(x) = dsqrt(x)

c				Initialize H(1,1)
      H(1,1) = one
c				perform (N-1) rotations
      do 50 i = 1, N-1
         k = i+1
c				check if we need to rotate in this plane
         if( D(i) .eq. zero ) then
            do 10 j = k, N
               H(i,j) = zero
  10        continue
            do 20 j = 1, i
               H(k,j) = zero
  20        continue
            H(k,k) = one
         else
c				compute rotation matrix (only need T11 and T21)
            x   = D(i)*D(i) + D(k)*D(k)
            x   = one/sqrt(x)
            T11 = x*D(k)
            T21 = x*D(i)
c				update transformation matrix H
            do 30 j = 1, i
               x      = H(i,j)
               H(i,j) = T11*x
               H(k,j) = T21*x
  30        continue
            H(i,k) = -T21
            H(k,k) =  T11
            do 40 j = i+2, N
               H(i,j) = zero
  40        continue
c				update unit vector D (relative to new basis)
            x    = D(i)
            D(i) = T11*x - T21*D(k)
            D(k) = T21*x + T11*D(k)
         endif
  50  continue

      return
      end
