c  *******************************************************************
c  *                                                                 *
c  *                       subroutine adot                           *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 03/11/78 1.01
c  Written by Jack Dongarra, Linpack
c  Latest Update: Jun 9, 1999
c
c  PURPOSE    forms the dot product of two vectors (LINPACK)
c
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-------------------------------------------------------------------------
      real function adot(n,dx,incx,dy,incy)
      real dx(*), dy(*), zero
      integer i, incx, incy, ix, iy, n
      data zero/0.0/
c
      adot = zero
      if( n .le. 0 ) return

      if( incx .eq. 1 .and. incy .eq. 1 ) then
         do 10 i = 1,n
            adot = adot + dx(i)*dy(i)
  10     continue
      else
         ix = 1
         iy = 1
         if( incx .lt. 0 ) ix = (-n+1)*incx + 1
         if( incy .lt. 0 ) iy = (-n+1)*incy + 1
         do 20 i = 1,n
            adot = adot + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
  20     continue
      endif

      return
      end
