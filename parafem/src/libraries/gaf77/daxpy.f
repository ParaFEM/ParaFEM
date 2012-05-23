c  *******************************************************************
c  *                                                                 *
c  *                       subroutine daxpy                          *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 03/11/78 1.01
c  Written by Jack Dongarra, Linpack.
c  Modified by Gordon A. Fenton, Aug. 24, 1993
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   constant times a vector plus a vector.
c
c  This routine computes the output vector
c
c         dy(1)        = dy(1)        + da*dx(1)
c         dy(2)        = dy(2)        + da*dx(2)
c         dy(3)        = dy(3)        + da*dx(3)
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
c   dy     on input, dy contains the additive vector. On output, dy contains
c          the result.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine daxpy(n,da,dx,dy)
      real*8 dx(*), dy(*), da, zero
      integer i, n
      data zero/0.d0/

      if( n .le. 0 .or. da .eq. zero ) return

      do 10 i = 1, n
         dy(i) = dy(i) + da*dx(i)
  10  continue

      return
      end
