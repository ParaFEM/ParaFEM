c  ***************************************************************
c  *                                                             *
c  *                  Subroutine pspct                           *
c  *                                                             *
c  ***************************************************************
c  Single Precision Version 1.22
c  Written by Gordon A. Fenton, Princeton, May 14, 1989.
c  Latest Update: Apr 26, 2006
c
c  PURPOSE  estimates power spectral density using Fourier transforms
c
c  Obtains the power spectral estimates from the discrete fourier
c  transform according to the procedure described by D.E. Newland in
c  "An Introduction to Random Vibrations and Spectral Analysis",
c  pg 142-145.  The input-ouput variables are;
c
c    X   real vector of length N. On input X contains the observed series
c        and on output X contains the power spectral estimates (one-sided
c        or two-sided, see ONESD) X(k) , k = 1,2, ..., N/2+1. Note that
c        X(1) becomes the power spectral estimate corresponding to zero
c        frequency and X(N/2+1) corresponds to the upper Nyquist frequency
c        1/2dt Hz or pi/dt rad/sec.  To construct the two-sided spectra from
c        this vector (for ONESD = .FALSE. ) simply set X(-k) = X(k).
c
c    A   real vector of length N. On output A contains the cosine fourier
c        coefficients A(k) as in
c
c                               N-1
c             X(t) = A(0) + 2 * sum { A(k)cos(w(k)*t) - B(k)sin(w(k)*t) }
c                               k=1
c
c        (output)
c
c    B   real vector of length N. On output B contains the sine fourier
c        coefficients B(k).
c
c    M   input integer such that N = 2**M
c
c    N   the length of vector X. N must be a power of 2.
c
c    L   the number of zeroes padding the end of vector X.
c
c   DT   the real input sampling interval between discrete values of X.
c
c  NSMTH integer input to indicate the number of adjacent spectral estimates to
c        average in the smoothing process. If NSMTH is zero, then no smoothing
c        is performed. If NSMTH is i, then the i previous and the i subsequent
c        spectral estimates are averaged together (2i+1 values in total) and
c        the result assigned to the central value. Averaging is performed
c        with equal weighting.
c
c  LZERO logical input flag. Set LZERO to .TRUE. if the series is to be zero
c        meaned and false if it is already a zero mean process prior to
c        entering this routine.
c
c  ONESD logical input flag. Set ONESD to .TRUE. if one-sided power spectral
c        estimates are desired (frequency range becomes 0 to pi/dt where dt is
c        the sampling interval). If ONESD is .FALSE., then the two sided
c        spectral density is returned with frequency range -pi/dt to pi/dt.
c        In the latter case, the spectral estimates are still only stored in
c        the first N/2+1 positions of X but is understood that X(-i) = X(i).
c        The one-sided spectra, G(w), is obtained from the (theoretical)
c        two-sided spectra, S(w), from the relationship
c
c                G(w) = 2*S(w), w >= 0.0
c
c        Note that the factor of 2 is applied also when w = 0, even though,
c        strictly speaking, the power associated with w = 0 is G(w)*dw/2 since
c        the frequency increment at w = 0 is half that used elsewhere. This
c        should be accounted for explicitly in any simulation procedure.
c
c  REVISION HISTORY:
c  1.2	eliminated unused variable `half' (Dec 5/96)
c  1.21	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c  1.22	corrected documentation above to read A(k)*cos - B(k)*sin (Apr 26/06)
c-----------------------------------------------------------------------------
      subroutine pspct( X, A, B, M, N, L, DT, NSMTH, LZERO, ONESD )
      dimension X(*), A(*), B(*)
      integer N, L, NSMTH, M
      logical LZERO, ONESD
      data pi/3.1415926535897932384/
      data zero/0.0/, one/1.0/, two/2.0/
c                                     zero-mean the series if desired
      if ( LZERO ) then
          u = averg(X,N)
          do 10 i = 1,N
             A(i) = X(i) - u
  10      continue
      else
          do 20 i = 1, N
             A(i) = X(i)
  20      continue
      endif
c                                     calculate the fourier coefficients

      call fft1d( A, B, M, .false. )

c                                     calculate the spectral coefficients
c                                    (X(1) corresponds to zero frequency)
      na = N/2
      nb = na + 1
      t = DT*float(N*N)/(two*pi*float(N-L))
      if (ONESD .and. NSMTH .eq. 0) t = two*t
      do 30 i = 1, nb
         X(i) = t*(A(i)*A(i) + B(i)*B(i))
  30  continue

c                                     smooth the continuous spectrum
      if ( NSMTH .eq. 0 ) return
      as = one/float(2*NSMTH + 1)
      if ( ONESD ) as = two*as
c                                     for k = 1;
      sum = X(1)
      do 40 i = 2, NSMTH+1
         sum = sum + two*X(i)
  40  continue
      X0 = as*sum
c                                     for k = 2 to NSMTH
      do 70 k = 2,NSMTH
         sum = X(1)
         do 50 i = 2, 2-k+NSMTH
            sum = sum + X(i)
  50     continue
         do 60 i = 2, k+NSMTH
            sum = sum + X(i)
  60     continue
         X(k+na) = as*sum
  70  continue
c                                     for k = NSMTH+1  to  N/2 - NSMTH + 1
      do 90 k = NSMTH+1, nb-NSMTH
         sum = zero
         do 80 i = k-NSMTH, k+NSMTH
            sum = sum + X(i)
  80     continue
         X(k+na) = as*sum
  90  continue
c                                     for k = N/2-NSMTH+2  to N/2
      do 110 k = nb-NSMTH+1, na
         sum = zero
         do 100 i = k-NSMTH, nb
            sum = sum + X(i)
 100     continue
         X(k+na) = as*sum
 110  continue
c                                     for k = N/2+1
      sum = zero
      do 120  i = nb-NSMTH, nb
         sum = sum + X(i)
 120  continue
      X(nb) = as*sum
c                     put the smoothed values back into X(i), i=1,2,...,N/2+1
c                     corresponding to frequencies 0 to pi/DT rad/sec.
      X(1) = X0
      do 130 k = 2, na
         X(k) = X(k+na)
 130  continue
c
      return
      end
