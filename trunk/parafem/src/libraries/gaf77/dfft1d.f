c  ******************************************************************
c  *                                                                *
c  *                   Subroutine dfft1d                            *
c  *                                                                *
c  ******************************************************************
c  Double Precision Version 1.12
c  Written by Gordon A. Fenton, Princeton, Jan. 6, 1988.
c  Latest Update: Apr 26, 2006
c
c  PURPOSE  compute the fast Fourier transform (or inverse) of a 1-D vector
c
c  Computes the fast Fourier transform of an observed series x(t)
c  which can be expressed as
c
c             x(t_i) = a(t_i) + i*b(t_i)
c
c  where a(t_i) and b(t_i) are the real and imaginary parts of x(t) at the
c  time t_i. The vectors A and B are overwritten with the calculated fourier
c  coefficients A(k) and B(k) corresponding to the frequencies 2*pi*k/(N*s)
c  where s is the time increment. This routine can also perform the inverse
c  fourier transform by setting the flag INV. The forward transform is
c  defined as
c                 1    N-1
c         X(k) = --- * sum x(t)*exp{-i*2*pi*k*t/N}
c                 N    t=0
c
c  where X(k) = A(k) + i*B(k) are the Fourier coefficients, and the inverse
c  transform is
c
c                 N-1
c         x(t) =  sum X(k)*exp{i*2*pi*k*t/N}	(complex case)
c                 k=0
c
c                           N-1
c              = A(0) + 2 * sum [A(k)*cos(w(k)*t) - B(k)*sin(w(k)*t)]
c                           k=1
c
c  if x(t) is real. x(t) is computed if INV = .true.
c
c  Input and output variables are;
c
c      A   real vector of length N. On input A contains the real part of the
c          observed series. On output A contains the fourier coefficients A(k)
c          corresponding to the frequencies 2*pi*k/(N*s) where s is the time
c          increment.
c
c      B   real vector of length N. On input B contains the imaginary
c          components of the observed series. On ouput, B contains the fourier
c          coefficients B(k).
c
c      m   the value of m is such that N = 2**m
c
c    INV   logical flag which is false if the fourier coefficients are to be
c          calculated, true if the inverse transform is required.
c
c   Note: this routine is modified from a routine presented by D.E. Newland
c         (originally written by J.W. Cooley et al.)
c
c  REVISION HISTORY:
c  1.1	eliminated argument 'N' to agree with real*4 version. (Jul 1/97)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c  1.12	corrected documentation above to read A(k)*cos - B(k)*sin (Apr 26/06)
c------------------------------------------------------------------------------
      Subroutine dfft1d( A, B, m, INV )
      implicit real*8 (a-h,o-z)
      dimension A(*), B(*)
      logical INV
      real*8 ninv
      data pi/3.141592653589793d0/, zero/0.d0/, one/1.d0/
      float(i) = dble(i)
c
      sygn = one
      N    = 2**m
      if ( .not. INV ) then
         sygn = -one
         ninv = one/float(N)
         do 10 j = 1,N
            A(j) = A(j)*ninv
            B(j) = B(j)*ninv
  10     continue
      endif
c                                      reorder vector
      nby2 = n/2
      j = 1
      do 30 l = 1, n-1
         if( l .lt. j ) then
             t = a(j)
             a(j) = a(l)
             a(l) = t
             t = b(j)
             b(j) = b(l)
             b(l) = t
         endif
         k = nby2
  20     if( k .lt. j ) then
             j = j - k
             k = k/2
             go to 20
         endif
         j = j + k
  30  continue
c                                 Calculate fourier coefficients
      me = 1
      do 50 mm = 1,m
         k = me
         pibyk = pi/dble(k)
         me = 2*me
         wr = dcos(pibyk)
         wi = sygn*dsin(pibyk)
c
         qr    = one
         qi    = zero
         do 50 j = 1, k
            do 40 l = j, N, me
               sr = A(l+k)*qr - B(l+k)*qi
               si = B(l+k)*qr + A(l+k)*qi
               A(l+k) = A(l) - sr
               B(l+k) = B(l) - si
               A(l) = A(l) + sr
               B(l) = B(l) + si
  40        continue
            tr = qr*wr - qi*wi
            ti = qr*wi + qi*wr
            qr = tr
            qi = ti
  50  continue
c
      return
      end
