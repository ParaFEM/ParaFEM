c  *******************************************************************
c  *                                                                 *
c  *                       subroutine iamax                          *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 06/17/77 1.01
c  Written by Jack Dongarra (?), Linpack
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   find the index of the maximum absolute element in a vector of
c            length n (LINPACK)
c
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c------------------------------------------------------------------------
      integer function iamax(n,dx,incx)
      real dx(*), dmax
      integer i, incx, ix, n

      if( n .le. 0 ) then
         iamax = 0
         return
      endif
      iamax = 1
      if( n .eq. 1 )return
      if( incx .eq. 1 ) then
         dmax = abs( dx(1) )
         do 10 i = 2, n
            if (abs( dx(i) ) .gt. dmax ) then
               iamax = i
               dmax   = abs( dx(i) )
            endif
  10     continue
      else
         ix = 1
         dmax = abs(dx(1))
         ix = ix + incx
         do 20 i = 2, n
            if (abs( dx(ix) ) .gt. dmax ) then
               iamax = i
               dmax   = abs(dx(ix))
            endif
            ix = ix + incx
  20     continue
      endif

      return
      end
