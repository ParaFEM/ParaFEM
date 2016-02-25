c  *********************************************************************
c  *                                                                   *
c  *                     Subroutine fft3g                              *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 2.11
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  generates a 3-D Gaussian quadrant symmetric random process
c           using the Fast Fourier Transform
c
c  This program generates a 3-D random field using the FFT method. The
c  field will be Gaussian, homogeneous, and quadrant symmetric with
c  structure defined by an external user defined function which returns
c  the spectral density of the process. Note that quadrant symmetry implies
c  that S(w1,w2,w3) = S(-w1,w2,w3) = S(w1,-w2,w3) = etc.
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
c       Z   real array of size at least N1 x N2 x N3 which on output will
c           contain the desired 3-D realization. (output)
c
c      N1   number of field points desired in the x direction. (input)
c
c      N2   number of field points desired in the y direction. (input)
c
c      N2   number of field points desired in the z direction. (input)
c
c      D1   physical length of the field in the x direction. (input)
c
c      D2   physical length of the field in the y direction. (input)
c
c      D3   physical length of the field in the z direction. (input)
c
c    PSFN   external function which returns the power spectral value at a
c           given frequency/wavenumber. The call to this function appears
c           as follows
c
c                  R = psfn( X, Y, Z )
c
c           where R is the power spectra at the frequency/wavenumber given
c           by the components (X,Y,Z).
c          
c   KSEED   on the first to call to this routine (or when INIT = 1), the
c           pseudo-random number generator is initialized using this
c           integer seed. If KSEED = 0, then a random seed is determined
c           using the clock time that the routine is entered. In any
c           case, the actual seed used is returned in the value of
c           KSEED (see the GAF77 library function ISEED). (input/output)
c
c    INIT   integer whose value is 1 if the parameters of the process
c           are to be initialized on this call. Subsequent calls for
c           other realizations of the same process should used INIT not
c           equal to 1. (input)
c
c      ZB   temporary real array of size at least N1 x N2 x N3.
c
c      iz   leading dimension of the arrays Z and ZB as specified in
c           the calling routine. (input)
c
c      jz   second dimension of the arrays Z and ZB as specified in the
c           calling routine. (input)
c
c      SD   real array of size at least N1 x N2 x (1+N3/2) which
c           will contain the standard deviations of the fourier coefficients.
c           If subsequent realizations of the same process are desired, SD
c           should not be modified between calls to this routine.
c           (input/output)
c
c      is   leading dimension of the array SD as specified in the calling
c           routine. (input)
c
c      js   second dimension of the array SD as specified in the calling
c           routine. (input)
c
c   t1,t2   two temporary real vectors each of length at least max(N2,N3).
c
c    IERR   integer flag which on return will have value 0 if all goes well.
c            = -1 if n1 is not a power of two
c            = -2 if n2 is not a power of two
c            = -3 if n1 and n2 are not powers of two
c            = -4 if n3 is not a power of two
c            = -5 if n1 and n3 are not powers of two
c            = -6 if n2 and n3 are not powers of two
c            = -7 if n1, n2, and n3 are not powers of two
c           (execution continues, but field size is set to the next lower
c           power of two)
c
c  Notes:
c   1) Timing information is broadcast through the block common /FFTTYM/: TI
c      is the time taken to initialize this generator and TS is the cumulative
c      generation time in seconds.
c   2) THIS ROUTINE HAS NOT BEEN THOROUGHLY TESTED. IF YOU HAVE PROBLEMS WITH
c      IT, LET ME KNOW. Feb. 4, 1993.
c      (gordon@random.field.tuns.ca)
c
c  Requires:
c   1) from libGAFsim:	VGAUS, ISEED, RANDF, NTWOM, SECOND, FFT3D
c   2) user provided external spectral density function.
c
c  REVISION HISTORY:
c  2.1	eliminated unused local variables `one', `two', and `pi' (Dec 5/96)
c  2.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c ======================================================================
      subroutine fft3g ( z, n1, n2, n3, d1, d2, d3, psfn, kseed, init,
     >                   zb, iz, jz, sd, is, js, t1, t2, ierr )
      dimension z(iz,jz,*), zb(iz,jz,*), sd(is,js,*), t1(*), t2(*)
      integer n1, n2, n3, init, kseed
      external psfn
      save n1b2, n2b2, n3b2, n1h, n2h, n3h, m1, m2, m3, iperr
      common/FFTTYM/ ti, ts
      data zero/0.0/, sixtnt/0.0625/, half/0.5/
      data twopi/6.2831853071795964769/
      data rttwo/1.4142135623730950488/
      data ifirst/1/
c --------------------------------------------- initialization ----------

      if( ifirst .eq. 1 .or. init .eq. 1 ) then
         ti     = second()
         ifirst = 0
         kseed  =  iseed( kseed )
c							saved parameters
         m1    = ntwom( n1 )
         m2    = ntwom( n2 )
         m3    = ntwom( n3 )
         n1b2  = n1/2
         n1h   = n1b2 + 1
         n2b2  = n2/2
         n2h   = n2b2 + 1
         n3b2  = n3/2
         n3h   = n3b2 + 1
         iperr = 0
         if( (2**m1) .ne. n1 ) iperr = -1
         if( (2**m2) .ne. n2 ) iperr = iperr - 2
         if( (2**m3) .ne. n3 ) iperr = iperr - 4
c							for known PSD...
         dw1   = twopi*float(n1 - 1)/(d1*float(n1))
         dw2   = twopi*float(n2 - 1)/(d2*float(n2))
         dw3   = twopi*float(n3 - 1)/(d3*float(n3))
         f0    = sixtnt*dw1*dw1*dw3
         f1    = half*f0
         w3a   = zero
         w3b   = zero
         do 30 k = 1, n3h
            f2  = half*f1
            w2a = zero
            w2b = zero
            do 20 j = 1, n2h
               f3  = half*f2
               w1a = zero
               w1b = zero
               do 10 i = 1, n1h
                  tm = f3*(psfn(w1a,w2a,w3a) + psfn(w1a,w2a,w3b)
     >                   + psfn(w1a,w2b,w3a) + psfn(w1a,w2b,w3b)
     >                   + psfn(w1b,w2a,w3a) + psfn(w1b,w2a,w3b)
     >                   + psfn(w1b,w2b,w3a) + psfn(w1b,w2b,w3b))
                  sd(i,j,k) = sqrt(tm)
                  f3  = f2
                  w1a = dw1*float(i)
                  w1b = dw1*float(n1-i)
  10           continue
               f2  = f1
               w2a = dw2*float(j)
               w2b = dw2*float(n2-j)
  20        continue
            do 25 j = 2, n2b2
               j2 = n2 - j + 2
               do 25 i = 2, n1b2
                  i2 = n1 - i + 2
                  sd(i2, j,k) = sd(i,j,k)
                  sd( i,j2,k) = sd(i,j,k)
                  sd(i2,j2,k) = sd(i,j,k)
  25        continue
            f1  = f0
            w3a = dw3*float(k)
            w3b = dw3*float(n3-k)
  30     continue
c							corner corrections...
         sd(1,1,1)       = rttwo*sd(1,1,1)
         sd(1,1,n3h)     = rttwo*sd(1,1,n3h)
         sd(1,n2h,1)     = rttwo*sd(1,n2h,1)
         sd(1,n2h,n3h)   = rttwo*sd(1,n2h,n3h)
         sd(n1h,1,1)     = rttwo*sd(n1h,1,1)
         sd(n1h,1,n3h)   = rttwo*sd(n1h,1,n3h)
         sd(n1h,n2h,1)   = rttwo*sd(n1h,n2h,1)
         sd(n1h,n2h,n3h) = rttwo*sd(n1h,n2h,n3h)

         ts = 0.
         ti = second() - ti
      endif
c ------------------------------------------------------------------------
c                         generate the field
      tt   = second()
      ierr = iperr
c                                              random fourier coefficients
      do 40 k = 1, n3h
      do 40 j = 1, n2
         call vgaus( Z(1,j,k),  sd(1,j,k), n1 )
         call vgaus( zb(1,j,k), sd(1,j,k), n1 )
  40  continue

c                                              introduce symmetry for real Z
      zb(1,1,1)       = zero
      zb(1,1,n3h)     = zero
      zb(1,n2h,1)     = zero
      zb(1,n2h,n3h)   = zero
      zb(n1h,1,1)     = zero
      zb(n1h,1,n3h)   = zero
      zb(n1h,n2h,1)   = zero
      zb(n1h,n2h,n3h) = zero
      do 50 k = 2, n3b2
         k2 = n3 - k + 2
         Z(1,1,k2)      =  Z(1,1,k)
         Z(1,n2h,k2)    =  Z(1,n2h,k)
         Z(n1h,1,k2)    =  Z(n1h,1,k)
         Z(n1h,n2h,k2)  =  Z(n1h,n2h,k)
         ZB(1,1,k2)     = -ZB(1,1,k)
         ZB(1,n2h,k2)   = -ZB(1,n2h,k)
         ZB(n1h,1,k2)   = -ZB(n1h,1,k)
         ZB(n1h,n2h,k2) = -ZB(n1h,n2h,k)
  50  continue
      do 60 j = 2, n2b2
         j2 = n2 - j + 2
         Z(1,j2,1)      =  Z(1,j,1)
         Z(1,j2,n3h)    =  Z(1,j,n3h)
         Z(n1h,j2,1)    =  Z(n1h,j,1)
         Z(n1h,j2,n3h)  =  Z(n1h,j,n3h)
         ZB(1,j2,1)     = -ZB(1,j,1)
         ZB(1,j2,n3h)   = -ZB(1,j,n3h)
         ZB(n1h,j2,1)   = -ZB(n1h,j,1)
         ZB(n1h,j2,n3h) = -ZB(n1h,j,n3h)
  60  continue
      do 70 i = 2, n1b2
         i2 = n1 - i + 2
         Z(i2,1,1)      =  Z(i,1,1)
         Z(i2,1,n3h)    =  Z(i,1,n3h)
         Z(i2,n2h,1)    =  Z(i,n2h,1)
         Z(i2,n2h,n3h)  =  Z(i,n2h,n3h)
         ZB(i2,1,1)     = -ZB(i,1,1)
         ZB(i2,1,n3h)   = -ZB(i,1,n3h)
         ZB(i2,n2h,1)   = -ZB(i,n2h,1)
         ZB(i2,n2h,n3h) = -ZB(i,n2h,n3h)
  70  continue
      do 80 k = 2, n3b2
         k2 = n3 - k + 2
         do 80 j = 2, n2b2
            j2 = n2 - j + 2
            Z(1,j2,k2)    =  Z(1,j,k)
            Z(1,j,k2)     =  Z(1,j2,k)
            Z(n1h,j2,k2)  =  Z(n1h,j,k)
            Z(n1h,j,k2)   =  Z(n1h,j2,k)
            ZB(1,j2,k2)   = -ZB(1,j,k)
            ZB(1,j,k2)    = -ZB(1,j2,k)
            ZB(n1h,j2,k2) = -ZB(n1h,j,k)
            ZB(n1h,j,k2)  = -ZB(n1h,j2,k)
  80  continue
      do 90 k = 2, n3b2
         k2 = n3 - k + 2
         do 90 i = 2, n1b2
            i2 = n1 - i + 2
            Z(i2,1,k2)    =  Z(i,1,k)
            Z(i,1,k2)     =  Z(i2,1,k)
            Z(i2,n2h,k2)  =  Z(i,n2h,k)
            Z(i,n2h,k2)   =  Z(i2,n2h,k)
            ZB(i2,1,k2)   = -ZB(i,1,k)
            ZB(i,1,k2)    = -ZB(i2,1,k)
            ZB(i2,n2h,k2) = -ZB(i,n2h,k)
            ZB(i,n2h,k2)  = -ZB(i2,n2h,k)
  90  continue
      do 100 j = 2, n2b2
         j2 = n2 - j + 2
         do 100 i = 2, n1b2
            i2 = n1 - i + 2
            Z(i2,j2,1)    =  Z(i,j,1)
            Z(i,j2,1)     =  Z(i2,j,1)
            Z(i2,j2,n3h)  =  Z(i,j,n3h)
            Z(i,j2,n3h)   =  Z(i2,j,n3h)
            ZB(i2,j2,1)   = -ZB(i,j,1)
            ZB(i,j2,1)    = -ZB(i2,j,1)
            ZB(i2,j2,n3h) = -ZB(i,j,n3h)
            ZB(i,j2,n3h)  = -ZB(i2,j,n3h)
 100  continue
      do 110 k = 2, n3b2
         k2 = n3 - k + 2
         do 110 j = 2, n2b2
            j2 = n2 - j + 2
            do 110 i = 2, n1b2
               i2 = n1 - i + 2
               Z(i2,j2,k2)  =  Z(i,j,k)
               Z(i2,j,k2)   =  Z(i,j2,k)
               Z(i,j2,k2)   =  Z(i2,j,k)
               Z(i,j,k2)    =  Z(i2,j2,k)
               ZB(i2,j2,k2) = -ZB(i,j,k)
               ZB(i2,j,k2)  = -ZB(i,j2,k)
               ZB(i,j2,k2)  = -ZB(i2,j,k)
               ZB(i,j,k2)   = -ZB(i2,j2,k)
 110  continue
c						now build the random field

      call fft3d( Z, iz, jz, ZB, iz, jz, m1, m2, m3, t1, t2, .true. )

      ts = ts + (second() - tt)
      return
      end
