c  *********************************************************************
c  *                                                                   *
c  *                       Subroutine izero                            *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.01
c  Written by Gordon A. Fenton, Princeton, Feb. 18, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  zeroes an integer vector of length n
c
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine izero( iv, n )
      integer n, iv(*)
c
      do 10 i = 1, n
         iv(i) = 0
  10  continue
c
      return
      end
