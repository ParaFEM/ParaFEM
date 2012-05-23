c  *******************************************************************
c  *                                                                 *
c  *                      Subroutine FFT2D                           *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, Mar. 23, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  compute the fast Fourier transform of a 2-D array using row-column
c           1-D transforms.
c
c  Calculates the discrete fourier transform of a two-dimensional
c  array A(N1,N2) + iB(N1,N2). The array A is overwritten with the A(k1,k2)
c  fourier coefficients and B(N1,N2) is used to store the B(k1,k2)
c  coefficients. The processing is carried out by subroutine fft1d.
c  Both N1 and N2 must be powers of 2 and so only the parameters m1, m2
c  where N1 = 2**m1, N2 = 2**m2, are required in this routine.
c
c  The input and output variables are:
c
c   A    real array of size (N1 x N2). On input A contains the observed
c        field values and on output A contains the first fourier coefficient.
c        (input/output)
c
c   ia   column dimension of array A as specified in the calling routine.
c        (input)
c
c   B    real array of size (N1 x N2). On output B contains the second
c        fourier coefficient. (If the inverse fourier transform is being run,
c        B contains the second fourier coefficients on input). (input/output)
c
c   ib   column dimension of array B as specified in the calling routine.
c        (input)
c
c   m1   integer such that N1 = 2**m1. (input)
c
c   m2   integer such that N2 = 2**m2. (input)
c
c ta,tb  temporary single precision vectors of length at least N2
c
c   INV  logical input flag. INV = .false. if normal fourier transform is being
c        run. INV = .true. if inverse fourier transform is being run ( in which
c        case, A and B contain fourier coefficients on input and A contains
c        'field' values on output). (input)
c
c   Note: this program is adapted from D.E. Newland.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      Subroutine fft2d( A, ia, B, ib, m1, m2, ta, tb, INV )
      logical INV
      real*4 A(ia,*), B(ib,*), ta(*), tb(*)

c                          first transform in the N2 direction
      N1 = 2**m1
      N2 = 2**m2
      do 20 i = 1, N1
c
         do 10 j = 1, N2
            ta(j) = A(i,j)
            tb(j) = B(i,j)
  10     continue
c
         call fft1d ( ta, tb, m2, INV )
c
         do 20 j = 1, N2
            A(i,j) = ta(j)
            B(i,j) = tb(j)
  20  continue
c
c  Now we can pass A and B by columns directly (this routine is most efficient
c  if N2 < N1.
c
      do 30 j = 1,N2
         call fft1d ( A(1,j), B(1,j), m1, INV )
  30  continue
c
      return
      end
