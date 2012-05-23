c  *********************************************************************
c  *                                                                   *
c  *                       subroutine fft2g                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 2.11
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  generates a 2-D Gaussian quadrant symmetric random process
c           using the Fast Fourier Transform
c
c  This program generates a 2-d random field using the FFT method. The
c  field will be Gaussian, homogeneous, and quadrant symmetric with
c  structure defined by an external user defined function which returns
c  the spectral density of the process. Note that quadrant symmetry implies
c  that S(w1,w2) = S(-w1,w2) = S(w1,-w2) = etc.
c  The process mean will be zero (over the ensemble).
c  Note that the final field will have a symmetric covariance
c  structure due to the periodicity of the FFT (about the midpoint of the
c  field). To improve representation of the total variance of the process,
c  variances associated with frequencies above the Nyquist limits (symmetric
c  about pi/Dx, where Dx is the sampling interval) are added to the variances
c  associated with frequencies below the Nyquist limit. Thus the spectral
c  content between 0 and 2*pi/Dx is included in the process. This routine
c  does not check to ensure that the spectral content beyond 2*pi/Dx is
c  insignificant - this must be ensured by the user.
c  Arguments to the function are described as follows,
c
c       Z   real array of size at least N1 x N2 which on output will contain
c           the desired 2-D realization. (output)
c
c      N1   number of field points desired in the x direction. (input)
c
c      N2   number of field points desired in the y direction. (input)
c
c      D1   physical length of the field in the x direction. (input)
c
c      D2   physical length of the field in the y direction. (input)
c
c    PSFN   external function which returns the power spectral value at a
c           given frequency/wavenumber. The call to this function appears
c           as follows
c
c                  R = psfn( X, Y )
c
c           where R is the power spectra at the frequency/wavenumber given
c           by the components (X,Y).
c          
c   KSEED   on the first to call to this routine (or when INIT = 1), the
c           pseudo-random number generator is initialized using this
c           integer seed. If KSEED = 0, then a random seed is determined
c           using the clock time that the routine is entered. In any
c           case, the actual seed used is returned in the value of
c           KSEED (see the GAF77 library function ISEED). (input/output)
c
c    INIT   integer whose absolute value is 1 if the parameters of the
c           process are to be initialized on this call. Subsequent calls
c           for other realizations of the same process should used |INIT|
c           not equal to 1. Note that if INIT = -1, only the initialization
c           is performed (the realization is not computed). (input)
c
c      ZB   temporary real array of size at least N1 x N2.
c
c      iz   column dimension of the arrays Z and ZB as specified in
c           the calling routine. (input)
c
c      SD   real array of size at least (1+N1/2) x (1+N2/2) which will contain
c           the standard deviations of the fourier coefficients. If subsequent
c           realizations of the same process are desired, SD should not be
c           modified between calls to this routine. (input/output)
c
c      is   column dimension of the array SD as specified in the calling
c           routine. (input)
c
c   t1,t2   two temporary real vectors each of length at least max(N1,N2).
c
c    ierr   integer flag which will have value 0 on return if all goes well.
c            = -1 if n1 is not a power of two
c            = -2 if n2 is not a power of two
c            = -3 if n1 and n2 are not powers of two
c           (execution continues, but field size is set to the next lower
c           power of two)
c
c  Notes:
c   1) Timing information is broadcast through the block common /FFTTYM/: TI
c      is the time taken to initialize this generator and TS is the cumulative
c      generation time in seconds.
c
c  Requires:
c   1) from libGAFsim:	VGAUS, ISEED, RANDF, NTWOM, SECOND, FFT2D
c   2) user provided external spectral density function.
c
c  REVISION HISTORY:
c  2.1	eliminated unused local variables `one', `two', and `pi' (Dec 5/96)
c  2.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c ======================================================================
      subroutine fft2g( z, n1, n2, d1, d2, psfn, kseed, init,
     >                  zb, iz, sd, is, t1, t2, ierr )
      dimension z(iz,*), zb(iz,*), sd(is,*), t1(*), t2(*)
      external psfn
      common/ffttym/ ti, ts
      save n1b2, n2b2, n1h, n2h, m1, m2, iperr
      data zero/0.0/, eighth/0.125/, half/0.5/
      data twopi/6.2831853071795864769/
      data rttwo/1.4142135623730950488/
      data ifirst/1/
c --------------------------------------------------------------------------
c					initialization

      if( ifirst .eq. 1 .or. init .eq. 1 ) then
         ti     = second()
         ifirst = 0
         kseed  = iseed( kseed )
         iperr  = 0
         m1     = ntwom(n1)
         m2     = ntwom(n2)
         if( (2**m1) .ne. n1 ) iperr = -1
         if( (2**m2) .ne. n2 ) iperr = iperr - 2
         n1b2   = n1/2
         n2b2   = n2/2
         n1h    = 1 + n1b2
         n2h    = 1 + n2b2
c						compute coefficient std. dev's
         dw1 = twopi*float(n1-1)/(d1*float(n1))
         dw2 = twopi*float(n2-1)/(d2*float(n2))
         f0  = eighth*dw1*dw2
         f1  = half*f0
         w2a = zero
         w2b = zero
         do 20 j = 1, n2h
            f2  = half*f1
            w1a = zero
            w1b = zero
            do 10 i = 1, n1h
               sd(i,j) = sqrt(f2*(psfn(w1a,w2a) + psfn(w1a,w2b)
     >                          + psfn(w1b,w2a) + psfn(w1b,w2b)))
               f2  = f1
               w1a = float(i)*dw1
               w1b = float(n1-i)*dw1
  10        continue
            do 15 i = 2, n1b2
               sd(n1-i+2,j) = sd(i,j)
  15        continue
            f1  = f0
            w2a = float(j)*dw2
            w2b = float(n2-j)*dw2
  20     continue
c						corner corrections...
         sd(1,1)     = rttwo*sd(1,1)
         sd(n1h,1)   = rttwo*sd(n1h,1)
         sd(1,n2h)   = rttwo*sd(1,n2h)
         sd(n1h,n2h) = rttwo*sd(n1h,n2h)

         ts = 0.
         ti = second() - ti
      endif
c --------------------------------------------------------------------------
c				generate the field
      tt   = second()
      ierr = iperr
c				generate random fourier coefficients by columns
      do 30 j = 1, n2h
         call vgaus( Z(1,j),  sd(1,j), n1 )
         call vgaus( zb(1,j), sd(1,j), n1 )
  30  continue
c				introduce symmetry for real Z
      zb(1,1)     = zero
      zb(n1h,1)   = zero
      zb(1,n2h)   = zero
      zb(n1h,n2h) = zero
      do 50 i = 2, n1h
         i2 = n1 - i + 2
         Z(i2,1)  = Z(i,1)
         zb(i2,1) = -zb(i,1)
  50  continue
      do 60 j = 2, n2h
         j2 = n2 - j + 2
         Z(1,j2)  = Z(1,j)
         zb(1,j2) = -zb(1,j)
         do 60 i = 2, n1h
            i2        = n1 - i + 2
            Z(i2,j2)  =  Z(i,j)
            zb(i2,j2) = -zb(i,j)
            Z(i,j2)   =  Z(i2,j)
            zb(i,j2)  = -zb(i2,j)
  60  continue
c				construct the field by inverse FFT
      call fft2d( Z, iz, ZB, iz, m1, m2, t1, t2, .true. )

      ts = ts + second() - tt
      return
      end
