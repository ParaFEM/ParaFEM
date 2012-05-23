c   ***************************************************************
c   *                                                             *
c   *                  Subroutine dxspct                          *
c   *                                                             *
c   ***************************************************************
c   Double Precision Version 1.1
c   Written by Gordon A. Fenton, Princeton, May 4, 1989.
c   Latest Update: Dec. 1, 1998
c
c   PURPOSE  compute the power spectral and cross-spectral estimates of a pair
c            of time histories via Fourier transforms.
c
c   Obtains the power spectral and cross-spectral estimates from the discrete
c   fourier transform according to the procedure described by D.E. Newland in
c   "An Introduction to Random Vibrations and Spectral Analysis",
c   pg 142-145.  The input-ouput variables are;
c
c    X   real vector of length N. On input X contains the first observed series
c        and on output X contains its power spectral estimates (one-sided
c        or two-sided, see ONESD) X(k) , k = 1,2, ..., N/2+1. Note that
c        X(1) becomes the power spectral estimate corresponding to zero
c        frequency and X(N/2+1) corresponds to the upper Nyquist frequency
c        1/2dt Hz or pi/dt rad/sec.  To construct the two-sided spectra from
c        this vector (for ONESD = .FALSE. ) simply set X(-k) = X(k).
c
c    Y   real vector of length N. On input Y contains the second observed
c        series and on output Y contains its power spectral estimates (one-
c        sided or two-sided, see ONESD) Y(k) , k = 1,2, ..., N/2+1. Note that
c        Y(1) becomes the power spectral estimate corresponding to zero
c        frequency and Y(N/2+1) corresponds to the upper Nyquist frequency
c        1/2dt Hz or pi/dt rad/sec.  To construct the two-sided spectra from
c        this vector (for ONESD = .FALSE. ) simply set Y(-k) = Y(k).
c
c    XY  real vector of length N+2. On output, XY contains the cross-spectral
c        estimates of series X with Y. Since these are typically complex
c        values, every second value is the imaginary component. That is XY(1)
c        is the real component of the cross-spectra corresponding to zero
c        frequency and XY(2) is its respective imaginary component.
c
c    AX, AY   real vectors of length N. On output AX or AY contain the cosine
c             fourier coefficients A(k) defined by
c                                 X(t) = SUM { A(k)cosWkt + B(k)sinWkt }
c             for the series X and Y respectively.
c
c    BX, BY   real vectors of length N. On output BX or BY contain the sine
c             fourier coefficients B(k) for series X and Y respectively.
c
c    M   input integer such that N = 2**M
c
c    N   the length of vector X and Y. N must be a power of 2.
c
c    L   the number of zeroes padding the end of vectors X and Y (assumed to be
c        the same).
c
c   DT   the real input sampling interval between discrete values of X or Y.
c
c  NSMTH integer input to indicate the number of adjacent spectral estimates to
c        average in the smoothing process. If NSMTH is zero, then no smoothing
c        is performed. If NSMTH is i, then the i previous and the i subsequent
c        spectral estimates are averaged together (2i+1 values in total) and
c        the result assigned to the central value. Averaging is performed
c        with equal weighting.
c
c LZERO  logical input flag. Set LZERO to .TRUE. if the series are to be mean
c        zeroed and .FALSE. if they are already mean zeroed prior to entering
c        this routine.
c
c  ONESD logical input flag. Set ONESD to .TRUE. if one-sided power spectral
c        estimates are desired (frequency range becomes 0 to pi/dt where dt is
c        the sampling interval). If ONESD is .FALSE., then the two sided
c        spectral density is returned with frequency range -pi/dt to pi/dt.
c        In the latter case, the spectral estimates are still only stored in
c        the first N/2+1 positions of X but is understood that X(-i) = X(i).
c
c  REVISION HISTORY:
c  1.1	corrected argument list in call to dfft1d (Dec 1/98)
c---------------------------------------------------------------------------
      Subroutine dxspct ( X, Y, XY, AX, BX, AY, BY, M, N,
     *                   L, DT, NSMTH, LZERO, ONESD )
      implicit real*8 (a-h,o-z)
      dimension X(*), Y(*), XY(*), AX(*), BX(*), AY(*), BY(*)
      logical LZERO, ONESD
      data pi/3.1415926535897932384d0/
      data zero/0.d0/, half/0.5d0/, one/1.d0/, two/2.d0/
      float(i) = dble(i)
c                                             zero-mean the series
      if ( LZERO ) then
          u = daverg(X,N-L)
          v = daverg(Y,N-L)
          do 10 i = 1,N-L
          AY(i) = Y(i) - v
  10      AX(i) = X(i) - u
      else
          do 20 i = 1,N
          AY(i) = Y(i)
  20      AX(i) = X(i)
      endif
c                                       calculate the fourier coefficients

      call dfft1d( AX, BX, M, .false. )
      call dfft1d( AY, BY, M, .false. )

c  calculate the spectral coefficients (X(1) corresponds to zero frequency)

      na = N/2
      nb = na + 1
      t = DT*float(N*N)/(two*pi*float(N-L))
      if (ONESD .and. NSMTH .eq. 0) t = two*t
      do 30 i = 1, nb
         j = 2 * i
         Y(i)    = t*(AY(i)*AY(i) + BY(i)*BY(i))
         XY(j-1) = t*(AX(i)*AY(i) + BX(i)*BY(i))
         XY(j)   = t*(AY(I)*BX(i) - AX(i)*BY(i))
         X(i)    = t*(AX(i)*AX(i) + BX(i)*BX(i))
  30  continue
      if (ONESD .and. NSMTH .eq. 0) then
          X(1)  = half*X(1)
          Y(1)  = half*Y(1)
          XY(1) = half*XY(1)
          XY(2) = half*XY(2)
      endif
c                                            smooth the continuous spectrum
      if ( NSMTH .eq. 0 ) return
      as = one/float(2*NSMTH + 1)
c                                            for k = 1;
      sumx = X(1)
      sumy = Y(1)
      do 40 i = 2, NSMTH+1
         sumy = sumy + two*Y(i)
         sumx = sumx + two*X(i)
  40  continue
      X0 = as*sumx
      Y0 = as*sumy

      if ( ONESD ) as = two*as
c                                            for k = 2 to NSMTH
      do 70 k = 2,NSMTH
         sumx = X(1)
         sumy = Y(1)
         do 50 i = 2, 2-k+NSMTH
            sumy = sumy + Y(i)
            sumx = sumx + X(i)
  50     continue
         do 60 i = 2, k+NSMTH
            sumy = sumy + Y(i)
            sumx = sumx + X(i)
  60     continue
         Y(k+na) = as*sumy
         X(k+na) = as*sumx
  70  continue
c                                        for k = NSMTH+1  to  N/2 - NSMTH + 1
      do 90 k = NSMTH+1, nb-NSMTH
         sumx = zero
         sumy = zero
         do 80 i = k-NSMTH, k+NSMTH
            sumy = sumy + Y(i)
            sumx = sumx + X(i)
  80     continue
         Y(k+na) = as*sumy
         X(k+na) = as*sumx
  90  continue
c                                               for k = N/2-NSMTH+2  to N/2
      do 110 k = nb-NSMTH+1, na
         sumx = zero
         sumy = zero
         do 100 i = k-NSMTH, nb
            sumy = sumy + Y(i)
            sumx = sumx + X(i)
 100     continue
         Y(k+na) = as*sumy
         X(k+na) = as*sumx
 110  continue
c                                               for k = N/2+1
      sumx = zero
      sumy = zero
      do 120  i = nb-NSMTH, nb
         sumy = sumy + Y(i)
         sumx = sumx + X(i)
 120  continue
      X(nb) = as*sumx
      Y(nb) = as*sumy

c  put the smoothed values back into X(i), i=1,2,...,N/2+1
c  corresponding to frequencies 0 to pi/DT rad/sec.

      X(1) = X0
      Y(1) = Y0
      do 130 k = 2, na
         Y(k) = Y(k+na)
         X(k) = X(k+na)
 130  continue
c
c   repeat the above except now for the cross-spectra
c
c                                   for k = 1;
      sumr = XY(1)
      sumi = XY(2)
      do 140 i = 2, NSMTH+1
         j = 2 * i
         sumi = sumi + two*XY(j)
         sumr = sumr + two*XY(j-1)
 140  continue
      X0r = as*sumr
      X0i = as*sumi
c
      if ( ONESD ) as = two*as
c                                   for k = 2 to NSMTH
      do 170 k = 2,NSMTH
         sumr = XY(1)
         sumi = XY(2)
         do 150 i = 2, 2-k+NSMTH
            j = 2 * i
            sumi = sumi + XY(j)
            sumr = sumr + XY(j-1)
 150     continue
         do 160 i = 2, k+NSMTH
            j = 2 * i
            sumi = sumi + XY(j)
            sumr = sumr + XY(j-1)
 160     continue
c                           store in unused portions of X and Y
         Y(k+na) = as*sumi
         X(k+na) = as*sumr
 170  continue
c                                   for k = NSMTH+1  to  N/2 - NSMTH + 1
      do 190 k = NSMTH+1, nb-NSMTH
         sumr = zero
         sumi = zero
         do 180 i = k-NSMTH, k+NSMTH
            j = 2 * i
            sumi = sumi + XY(j)
            sumr = sumr + XY(j-1)
 180     continue
         Y(k+na) = as*sumi
         X(k+na) = as*sumr
 190  continue
c                                   for k = N/2-NSMTH+2  to N/2
      do 210 k = nb-NSMTH+1, na
         sumr = zero
         sumi = zero
         do 200 i = k-NSMTH, nb
            j = 2 * i
            sumi = sumi + XY(j)
            sumr = sumr + XY(j-1)
 200     continue
         Y(k+na) = as*sumi
         X(k+na) = as*sumr
 210  continue
c                                   for k = N/2+1
      sumr = zero
      sumi = zero
      do 220  i = nb-NSMTH, nb
         j = 2 * i
         sumi = sumi + XY(j)
         sumr = sumr + XY(j-1)
 220  continue
      XY(N+1) = as*sumr
      XY(N+2) = as*sumi
c
c  put the smoothed values back into XY(i), i=1,2,...,N+2
c  corresponding to frequencies 0 to pi/DT rad/sec with odd values real,
c  even values imaginary.
c
      XY(1) = X0r
      XY(2) = X0i
      do 230 k = 2, na
         j = 2 * k
         XY(j) = Y(k+na)
         XY(j-1) = X(k+na)
 230  continue
c
      return
      end
