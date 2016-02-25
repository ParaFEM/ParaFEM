c  *********************************************************************
c  *                                                                   *
c  *                        Subroutine minmx                           *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, Jan. 21, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  finds the minimum and maximum values in a vector.
c
c  Find the minimum and maximum of n values stored in z
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      subroutine minmx( z, n, zmax, zmin )
      integer n
      real z(*), zmax, zmin

      zmax = z(1)
      zmin = z(1)

      do 10 i = 2, n
         if( z(i) .gt. zmax ) zmax = z(i)
         if( z(i) .lt. zmin ) zmin = z(i)
  10  continue

      return
      end
