c  *******************************************************************
c  *                                                                 *
c  *                        subroutine axpy                          *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 03/11/78 1.01
c  Written by Jack Dongarra, Linpack
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   constant times a vector plus a vector (LINPACK)
c
c  This routine computes the output vector
c
c         dy(1)        = dy(1)        + da*dx(1)
c         dy(1+incy)   = dy(1+incy)   + da*dx(1+incx)
c         dy(1+2*incy) = dy(1+2*incy) + da*dx(1+2*incx)
c          .              .               .
c          .              .               .
c
c  Arguments are
c
c   n      the number of elements in the output vector
c
c   da     the constant
c
c   dx     the vector being multiplied by da
c
c   incx   the stride through the vector dx
c
c   dy     on input, dy contains the additive vector. On output, dy contains
c          the result.
c
c   incy   the stride through the vector dy
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine axpy(n,da,dx,incx,dy,incy)
      real dx(*), dy(*), da, zero
      integer i, incx, incy, ix, iy, n
      data zero/0.0/
c
      if( n .le. 0 )return
      if( da .eq. zero ) return
      if( incx .eq. 1 .and. incy .eq. 1 ) then
         do 10 i = 1, n
            dy(i) = dy(i) + da*dx(i)
  10     continue
      else
         ix = 1
         iy = 1
         if( incx .lt. 0 )ix = (-n+1)*incx + 1
         if( incy .lt. 0 )iy = (-n+1)*incy + 1
         do 20 i = 1, n
            dy(iy) = dy(iy) + da*dx(ix)
            ix = ix + incx
            iy = iy + incy
  20     continue
      endif

      return
      end
