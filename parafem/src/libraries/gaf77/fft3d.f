c  *******************************************************************
c  *                                                                 *
c  *                      Subroutine FFT3D                           *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, Mar. 23, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  computes the fast Fourier transform of a 3-D array using
c           row-column-leaf 1-D transforms.
c
c  Calculates the discrete fourier transform of a three-dimensional
c  array A(N1,N2,N3) + iB(N1,N2,N3). The array A is overwritten with the
c  A(k1,k2,k3) fourier coefficients and B(N1,N2,N3) is used to store the
c  B(k1,k2,k3) coefficients. The processing is carried out by subroutine
c  fft1d. N1, N2, and N3 must be powers of 2 and so only the parameters m1,
c  m2, and m3, where N1 = 2**m1, N2 = 2**m2, and N3 = 2**m3, are required
c  in this routine as input.
c
c  The input and output variables are:
c
c   A    real array of size (N1 x N2 x N3). On input A contains the observed
c        field values and on output A contains the first fourier coefficient
c        (cosine coefficients). (input/output)
c
c ia,ja  first and second index dimension of array A as specified in the
c        calling routine. (input)
c
c   B    real array of size (N1 x N2 x N3). On output B contains the second
c        fourier coefficient (sine coefficients). (If the inverse fourier
c        transform is being run, B contains the second fourier coefficients
c        on input). (input/output)
c
c ib,jb  first and second index dimension of array B as specified in the
c        calling routine. (input)
c
c   m1   integer such that N1 = 2**m1. (input)
c
c   m2   integer such that N2 = 2**m2. (input)
c
c   m3   integer such that N3 = 2**m3. (input)
c
c ta,tb  two temporary real vectors of length max(N2,N3)
c
c   INV  logical input flag. INV = .false. if normal fourier transform is being
c        run. INV = .true. if inverse fourier transform is being run ( in which
c        case, A and B contain fourier coefficients on input and A contains
c        'field' values on output).
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      Subroutine fft3d( A, ia, ja, B, ib, jb, m1, m2, m3, ta, tb, INV )
      dimension A(ia,ja,*), B(ib,jb,*), ta(*), tb(*)
      logical INV
c				first transform in the N3 direction
      N1 = 2**m1
      N2 = 2**m2
      N3 = 2**m3
      do 20 i = 1, N1
      do 20 j = 1, N2
         do 10 k = 1, N3
            ta(k) = A(i,j,k)
            tb(k) = B(i,j,k)
  10     continue

         call fft1d ( ta, tb, m3, INV )

         do 20 k = 1, N3
            A(i,j,k) = ta(k)
            B(i,j,k) = tb(k)
  20  continue
c				now transform in the N2 direction
      do 40 i = 1, N1
      do 40 k = 1, N3
         do 30 j = 1, N2
            ta(j) = A(i,j,k)
            tb(j) = B(i,j,k)
  30     continue

         call fft1d ( ta, tb, m2, INV )

         do 40 j = 1, N2
            A(i,j,k) = ta(j)
            B(i,j,k) = tb(j)
  40  continue
c				finally transform in the N1 direction
      do 50 j = 1, N2
      do 50 k = 1, N3
         call fft1d ( A(1,j,k), B(1,j,k), m1, INV )
  50  continue

      return
      end
