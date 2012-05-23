c  *******************************************************************
c  *                                                                 *
c  *                       subroutine dswap                          *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 06/17/77 1.01
c  Written by Jack Dongarra, Linpack
c  Modified by Gordon A. Fenton, Aug. 24. 1993
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   interchanges two vectors.
c
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-------------------------------------------------------------------------
      subroutine  dswap (n,dx,dy)
      real*8 dx(*), dy(*), dtemp
      integer i, m, mp1, n

      if( n .le. 0 ) return

      m = mod(n,3)
      if( m .ne. 0 ) then
         do 10 i = 1, m
           dtemp = dx(i)
           dx(i) = dy(i)
           dy(i) = dtemp
  10     continue
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
  20  continue

      return
      end
