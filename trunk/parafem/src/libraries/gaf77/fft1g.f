c  *********************************************************************
c  *                                                                   *
c  *                       Subroutine fft1g                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 2.11
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  generates a 1-D Gaussian random process using the Fast Fourier
c           Transform
c
c  This program generates a 1-D random process with given spectral density
c  function via the FFT algorithm. Note that the output process will have
c  a symmetric covariance structure (about the field mid- point). To improve
c  the representation of the total variance of the process, frequencies
c  symmetrically above the Nyquist frequency (pi/Dx, where Dx is the
c  sampling interval) are folded into the calculation of the variance
c  associated with each frequency. Thus the spectral content between 0 and
c  2*pi/Dx is included in the process. The spectral content between 2*pi/Dx
c  and 4*pi/Dx is checked (if init = -1) to ensure that it is less than 10%
c  of that between 0 and 2*pi/Dx (if not, ierr is set to -10). This is not
c  to be considered a full check, but may indicate a possible source of
c  error. The user should ensure that the full spectral content is
c  adequately represented by the portion between 0 and 2*pi/Dx. Only a
c  simple trapezoidal integration scheme is used to check spectral content.
c  Arguments are as follows,
c
c       Z     real vector of length at least N1 which on output will
c             contain the desired real random sequence. (output)
c
c      N1     number of points in the desired process (must be a power of 2).
c             This number should not change from call to call unless INIT is
c             set to 1 at the time of the change. (input)
c
c      D1     physical length of the process. (input)
c
c    PSFN     external function which computes the (one-sided) power spectral
c             value at a given frequency/wavenumber. The call to this function
c             appears as
c
c                         V = psfn( x )
c
c             where x is the frequency/wave number.
c
c   KSEED     on the first to call to this routine (or when INIT = 1), the
c             pseudo-random number generator is initialized using this
c             integer seed. If KSEED = 0, then a random seed is determined
c             using the clock time that the routine is entered. In any
c             case, the actual seed used is returned in the value of
c             KSEED (see the GAF77 library function ISEED). (input/output)
c
c    INIT     integer whose absolute value is 1 if the parameters of the
c             process are to be initialized on this call. Subsequent calls
c             for other realizations of the same process should used |INIT|
c             not equal to 1. Note that if INIT = -1, only the initialization
c             is performed (the realization is not computed). In this case,
c             the  power spectral content below (2*pi/Dx) is compared to that
c             above (2*pi/Dx) and an error code returned if the latter exceeds
c             10% of the former. (input)
c
c      ZB     temporary real vector of length at least N1 used as workspace.
c
c      SD     real vector of length at least 1+N1/2 which will contain the
c             standard deviations of the Fourier coefficients so that
c
c                     Var[A_1] = SD(1)^2
c                     Var[A_k] = Var[B_k] = (SD(k))^2,  k = 2, 3, ..., N1/2
c                     Var[A_(1+N1/2)] = SD(1+N1/2)
c
c             where we use the fact that A(N1-k+2) = A(k), k = 2, 3,..,N1/2
c             and B(N1-k+2) = -B(k), k = 2, 3, ..,N1/2. Also B(1) = 0 and
c             B(N1/2+1) = 0 if the process is to be real (where A and B are
c             the vectors of Fourier coefficients X = A + i*B). If subsequent
c             realizations of the same process are desired, then the
c             values stored in SD should not be changed between calls
c             to this routine. (input/output)
c
c    IERR     integer flag which on return will have value 0 if all goes well.
c              = -1   if n1 is not a power of two - execution continues, but
c                     uses the next lower power of two for the field size.
c              = -10  if spectral content beyond (2*pi/Dx) exceeds 10%
c                     (see pstol).
c              = -11  if both of the above errors occur.
c
c  Notes:
c   1) Timing information is broadcast through the block common /FFTTYM/: TI
c      is the time taken to initialize this generator and TS is the cumulative
c      generation time in seconds.
c
c  Requires:
c   1) from libGAFsim:	VGAUS, ISEED, RANDF, NTWOM, SECOND, FFT1D, TPINT
c   2) user provided external spectral density function.
c
c  REVISION HISTORY:
c  2.1	eliminated unused local variables `one' and `pi' (Dec 5/96)
c  2.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c ===========================================================================
      subroutine fft1g( z, n1, d1, psfn, kseed, init, zb, sd, ierr )
      dimension z(*), zb(*), sd(*)
      external psfn
      common/ffttym/ ti, ts
      save n1b2, n1h, m1, iperr
      data pstol/0.1/
      data zero/0.0/, quart/0.25/, half/0.5/, two/2.0/
      data twopi/6.2831853071795864769/
      data ifirst/1/
c ------------------------------------- initialization -----------------
      if( ifirst .eq. 1 .or. iabs(init) .eq. 1 ) then
         ti     = second()
         ifirst = 0
         kseed  = iseed(kseed)
         m1     = ntwom(n1)
         iperr  = 0
         if( (2**m1) .ne. n1 ) iperr = -1
         n1b2   = n1/2
         n1h    = 1 + n1b2
c							for known PSD ...
         dw      = twopi*float(n1-1)/(d1*float(n1))
         qdw     = quart*dw
         sd(1)   = sqrt(half*dw*psfn(zero))
         sd(n1h) = sqrt(dw*psfn(float(n1b2)*dw))
         do 40 i = 2, n1b2
            w1    = dw*float(i-1)
            w2    = dw*float(n1-i+1)
            sd(i) = sqrt(qdw*(psfn(w1) + psfn(w2)))
  40     continue

         ts = 0.
         if( init .eq. -1 ) then
c							check spectral content
            wm1 = dw*float(n1-1)
            wm2 = two*wm1
            psint = tpint( psfn, zero, wm1, n1 )
	    psext = tpint( psfn, wm1, wm2, n1 )
            if( abs(psext/psint) .gt. pstol ) iperr = iperr - 10
            ierr = iperr
            ti = second() - ti
            return
         endif
         ti = second() - ti
      endif
c ------------------------------------------------------------------------
c                      generate the process
      tt   = second()
      ierr = iperr
c					generate random fourier coefficients
      call vgaus( Z,  sd, n1h )
      call vgaus( ZB, sd, n1b2 )
c					introduce symmetry
      ZB(1)   = zero
      ZB(n1h) = zero
      do 50 i = 2, n1b2
         Z (n1-i+2) =  Z(i)
         ZB(n1-i+2) = -ZB(i)
  50  continue
c					generate the field
      call fft1d( Z, ZB, m1, .true. )

      ts = ts + (second() - tt)

      return
      end
