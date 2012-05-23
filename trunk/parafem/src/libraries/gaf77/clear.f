c  *********************************************************************
c  *                                                                   *
c  *                       Subroutine clear                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, Feb. 18, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  sets a real vector of length n to zero
c
c  Arguments to the routine are as follows;
c
c    x   real vector of length at least n which is set to contain 0's. (output)
c
c    n   number of elements of x to zero. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine clear( x, n )
      dimension x(*)
      data zero/0.0/
c
      do 10 i = 1, n
         x(i) = zero
  10  continue
c
      return
      end
