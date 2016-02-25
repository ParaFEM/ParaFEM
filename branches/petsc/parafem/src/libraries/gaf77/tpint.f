c  ************************************************************************
c  *                                                                      *
c  *                          Function tpint                              *
c  *                                                                      *
c  ************************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, Princeton, December 1989
c
c  PURPOSE  returns the integral of a user supplied function by trapezoidal
c           peicewise approximation.
c
c  Computes the integral of a user supplied function `fnc' between the bounds
c  `a' and `b' by discretizing over the range and summing the area under the
c  piece-wise linear approximation. The number of sampling points is given
c  by `n'.
c-----------------------------------------------------------------------------
      real function tpint( fnc, a, b, n )
      external fnc
      data zero/0.0/, half/0.5/

      dx = (b - a)/float(n-1)
      f1 = fnc(a)
      s  = zero
      do 10 i = 1, n
         x  = a + float(i)*dx
         f2 = fnc(x)
         s  = s + (f1 + f2)
         f1 = f2
  10  continue

      tpint = half*dx*s

      return
      end
