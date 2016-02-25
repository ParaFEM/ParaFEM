c  *******************************************************************
c  *                                                                 *
c  *                       subroutine addot                          *
c  *                                                                 *
c  *******************************************************************
c  Single/Mixed Precision Version 03/11/78 1.01
c  Written by Jack Dongarra, Linpack
c  Latest Update: Jun 9, 1999
c
c  PURPOSE    forms the dot product of two vectors (LINPACK)
c             (high precision version)
c
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-------------------------------------------------------------------------
      real function addot(n,dx,incx,dy,incy)
      real dx(*), dy(*), zero
      real*8 dble, s
      integer i, incx, incy, ix, iy, n
      data zero/0.0/
c
      addot = zero
      if( n .le. 0 ) return

      s = 0.d0
      if( incx .eq. 1 .and. incy .eq. 1 ) then
         do 10 i = 1,n
            s = s + dble(dx(i))*dble(dy(i))
  10     continue
      else
         ix = 1
         iy = 1
         if( incx .lt. 0 ) ix = (-n+1)*incx + 1
         if( incy .lt. 0 ) iy = (-n+1)*incy + 1
         do 20 i = 1,n
            s  = s + dble(dx(ix))*dble(dy(iy))
            ix = ix + incx
            iy = iy + incy
  20     continue
      endif
c					convert back to single precision
      addot = s

      return
      end
