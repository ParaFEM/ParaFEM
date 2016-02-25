c  *******************************************************************
c  *                                                                 *
c  *                       subroutine aswap                          *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 06/17/77 1.01
c  Written by Jack Dongarra, Linpack
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   interchanges two vectors (LINPACK)
c
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-------------------------------------------------------------------------
      subroutine  aswap (n,dx,incx,dy,incy)
      real dx(*), dy(*), dtemp
      integer i, incx, incy, ix, iy, m, mp1, n

      if( n .le. 0 ) return
      if( incx .eq. 1 .and. incy .eq. 1 ) then
         m = mod(n,3)
         if( m .ne. 0 ) then
            do 10 i = 1, m
              dtemp = dx(i)
              dx(i) = dy(i)
              dy(i) = dtemp
  10        continue
            if( n .lt. 3 ) return
         endif
         mp1 = m + 1
         do 20 i = mp1, n, 3
            dtemp   = dx(i)
            dx(i)   = dy(i)
            dy(i)   = dtemp
            dtemp   = dx(i + 1)
            dx(i+1) = dy(i + 1)
            dy(i+1) = dtemp
            dtemp   = dx(i + 2)
            dx(i+2) = dy(i + 2)
            dy(i+2) = dtemp
  20     continue
      else
         ix = 1
         iy = 1
         if( incx .lt. 0 ) ix = (-n+1)*incx + 1
         if( incy .lt. 0 ) iy = (-n+1)*incy + 1
         do 30 i = 1, n
            dtemp  = dx(ix)
            dx(ix) = dy(iy)
            dy(iy) = dtemp
            ix     = ix + incx
            iy     = iy + incy
  30     continue
      endif

      return
      end
