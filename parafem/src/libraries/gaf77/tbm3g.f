c  *******************************************************************
c  *                                                                 *
c  *                      Subroutine tbm3g                           *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 3.21
c  Written by Gordon A. Fenton, TUNS, Dec. 22, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  generates a three-dimensional random field using the Turning
c           Bands Method.
c
c  A three-dimensional random field is generated here using the turning
c  bands algorithm original developed by Matheron. This implementation
c  employs all eight corners of the field as origins of the 1-d basis
c  line processes - the corner used depends on the orientation
c  of the line.
c  In order to produce the three-dimensional field, a sequence of 1-d
c  processes are generated at random or deterministic orientations.
c  The 1-d process is defined through an integral equation involving
c  the desired 3-d spectral density function. This 1-d function must be
c  known prior to calling this routine and its spectral properties are
c  passed to this routine via the external function referred to here
c  as "teq1d". Realizations of the 1-D process are produced through calls
c  to fft1g, which uses the FFT technique. To avoid periodicity problems
c  in the 1-D covariance structure, the line process produced is twice as
c  long as required - the second half discarded. To extend the frequency
c  range considered, the resolution along the line process is twice that
c  of the destination field (recall the Nyquist frequency is pi/Dx, where
c  Dx is the distance between points along the line process).
c  To improve representation of the total variance of the line process, fft1g
c  adds variances associated with frequencies above the Nyquist limits
c  (symmetric about pi/Dx, where Dx is the sampling interval) to the variances
c  associated with frequencies below the Nyquist limit. Thus the spectral
c  content between 0 and 2*pi/Dx is included in the process. The spectral
c  content between 2*pi/Dx and 4*pi/Dx is checked to ensure that it does
c  not exceed 10% of the variance between 0 and 2*pi/Dx. (see ierr and fft1g).
c  This is not a perfect check, but may give an indication of a possible
c  source of error.
c  Arguments to this routine are as follows;
c
c       Z   real array of size at least N1 x N2 x N3 which on output will
c           contain the desired 3-D realization. (output)
c
c      N1   number of field points desired in the x direction. (input)
c
c      N2   number of field points desired in the y direction. (input)
c
c      N3   number of field points desired in the z direction. (input)
c
c      D1   physical length of the field in the x direction. (input)
c
c      D2   physical length of the field in the y direction. (input)
c
c      D3   physical length of the field in the z direction. (input)
c
c   TEQ1D   external function which returns the power spectral value of the
c           equivalent 1-D process at a given frequency/wavenumber. The call
c           to this function appears as follows
c
c                  R = teq1d( X )
c
c           where R is the power spectra at the frequency/wavenumber given
c           by X.
c          
c   KSEED   on the first to call to this routine (or when INIT = 1), the
c           pseudo-random number generator is initialized using this
c           integer seed. If KSEED = 0, then a random seed is determined
c           using the clock time that the routine is entered. In any
c           case, the actual seed used is returned in the value of
c           KSEED (see the GAF77 library function ISEED). (input/output)
c
c    INIT   integer whose absolute value is 1 if the parameters of the
c           process are to be initialized on this call (which should be
c           done if any of the previous parameters, and/or L, are changed).
c           Subsequent calls for other realizations of the same process
c           should used |INIT| not equal to 1. Note that if INIT = -1,
c           only the initialization is performed (the realization is not
c           computed). (input)
c
c      iz   leading dimension of the array Z as specified in the calling
c           routine. (input)
c
c      jz   second dimension of the array Z as specified in the calling
c           routine. (input)
c
c       L   number of lines (or bands) to use in the generation. (input)
c
c   lrand   logical flag which is true if the lines are to be randomly
c           oriented on the sphere. Otherwise, if L is one of 6, 14, 26, or
c           50, lrand = false results in lines which are approximately
c           evenly distributed on the sphere (for values of L other than
c           those specified, results will be unpredictable). (input)
c
c      sd   real vector of length at least 1 + 2*max(n1,n2,n3) used to contain
c           the standard deviations associated with each frequency
c           representing the 1-D equivalent process. Note that this vector is
c           reused on subsequent calls to produce new realizations, so should
c           not be touched by the calling routine(s). (output/input)
c
c     uni   temporary real vector of length at least 4*max(n1,n2,n3) used to
c           store the 1-D process generated by fft1g.
c
c      t1   temporary real vector of length at least 4*max(n1,n2,n3) used by
c           fft1g.
c
c    IERR   integer flag which on return will have value 0 if all goes well.
c            = -1   if 4*max(n1,n2,n3) is not a power of two - execution
c                   continues, but uses the next lower power of two for the
c                   1-D process length.
c            = -2   if L > 50 (it is reset to 50).
c            = -10  if 1-D spectral content between (2*pi/Dx) and (4*pi/Dx)
c                   exceeds 10% of that between 0 and (2*pi/Dx).
c            = sum  for combinations of the above errors.
c
c  Notes:
c   1) Timing information is broadcast through the block common /TBMTYM/: TI
c      is the time taken to initialize this generator and TS is the cumulative
c      generation time in seconds.
c   2) there is little penalty associated with using np = 4*n rather than
c      just 2*n (about 10% slower). Most of the time is lost in the projection
c      onto the field process.
c
c  Requires:
c   1) from libGAFsim:	ISEED, FFT1G, VGAUS, RANDF, SECOND, NTWOM, FFT1D, TPINT
c   2) user defined external spectral density function (see TEQ1D)
c
c  REVISION HISTORY:
c  3.1	now using new randu function (RAN2 from Numerical Recipes, 2nd Ed)
c	(Oct 14/96)
c  3.2	eliminated unused local variable `half' (Dec 5/96)
c  3.21	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c ======================================================================
      subroutine tbm3g( z, n1, n2, n3, d1, d2, d3, teq1d, kseed, init,
     >                  iz, jz, L, lrand, sd, uni, t1, ierr )
      logical lrand
      real z(iz,jz,*), uni(*), sd(*), t1(*)
      real xstart(2), ystart(2), zstart(2), e1(50), e2(50), e3(50)
      save np, ul, xstart, ystart, zstart, dx, dy, dz, sc
      common/TBMTYM/ ti, ts
      data zero/0.0/, one/1.0/, two/2.0/
      data twopi/6.283185307/
      data ifirst/0/

c  Contains the data for the approximately evenly distributed unit
c  vectors on the unit sphere. To use properly, L should be one of
c  6, 14, 26, or 50 (but not greater than 50).

      data e1/1.00000,-1.000000, 0.000000, 0.000000, 0.000000, 0.000000,
     >       0.577350,-0.577350,-0.577350, 0.577350, 0.577350,-0.577350,
     >      -0.577350, 0.577350, 0.707107,-0.707107,-0.707107, 0.707107,
     >       0.707107,-0.707107,-0.707107,-0.707107, 0.000000, 0.000000,
     >       0.000000, 0.000000, 0.816497, 0.816497, 0.408248, 0.408248,
     >       0.408249, 0.408249,-0.408249,-0.408249,-0.408248,-0.408248,
     >      -0.816497,-0.816497,-0.816497,-0.816497,-0.408248,-0.408248,
     >      -0.408249,-0.408249, 0.408249, 0.408249, 0.408248, 0.408248,
     >       0.816496, 0.816496/
      data e2/0.00000, 0.000000, 1.000000,-1.000000, 0.000000, 0.000000,
     >       0.577350, 0.577350, 0.577350, 0.577350,-0.577350,-0.577350,
     >      -0.577350,-0.577350, 0.707107, 0.707107,-0.707107,-0.707107,
     >       0.000000, 0.000000, 0.000000, 0.000000, 0.707107, 0.707107,
     >      -0.707107,-0.707107, 0.408248,-0.408248, 0.816496,-0.816496,
     >       0.408248,-0.408248, 0.408248,-0.408248, 0.816496,-0.816496,
     >       0.408248,-0.408248, 0.408248,-0.408248, 0.816496,-0.816496,
     >       0.408248,-0.408248, 0.408248,-0.408248, 0.816496,-0.816496,
     >       0.408248,-0.408248/
      data e3/0.00000, 0.000000, 0.000000, 0.000000, 1.000000,-1.000000,
     >       0.577350, 0.577350,-0.577350,-0.577350, 0.577350, 0.577350,
     >      -0.577350,-0.577350, 0.000000, 0.000000, 0.000000, 0.000000,
     >       0.707107, 0.707107, 0.707107,-0.707107, 0.707107,-0.707107,
     >      -0.707107, 0.707107, 0.408249, 0.408249, 0.408248, 0.408248,
     >       0.816496, 0.816496, 0.816497, 0.816497, 0.408248, 0.408248,
     >       0.408249, 0.408249,-0.408249,-0.408249,-0.408248,-0.408248,
     >      -0.816497,-0.816497,-0.816497,-0.816497,-0.408248,-0.408248,
     >      -0.408249,-0.408249/
c
c--------------------------------------------- initialize the generator ---
      if( iabs(init) .eq. 1 .or. ifirst .eq. 0 ) then
         ti        = second()
         ifirst    = 1
         dx        = d1/float(n1 - 1)
         dy        = d2/float(n2 - 1)
         dz        = d3/float(n3 - 1)
         xstart(1) = -dx
         xstart(2) = -d1 - dx
         ystart(1) = -dy
         ystart(2) = -d2 - dy
         zstart(1) = -dz
         zstart(2) = -d3 - dz
         np        = 4*max0( n1, n2, n3 )
         ul        = two*sqrt( d1*d1 + d2*d2 + d3*d3 )
         sc        = float(np-1)/ul
         call fft1g( uni, np, ul, teq1d, kseed, -1, t1, sd, ierr )
         if( .not. lrand .and. L .gt. 50 ) then
            L = 50
            ierr = ierr - 2
         endif
c						normalize sd wrt L
         div = one/sqrt(float(L))
         do 10 i = 1, 1+np/2
            sd(i) = div*sd(i)
  10     continue
         ts = zero
         ti = second() - ti
         if( init .eq. -1 ) return
      endif
c-------------------------------------- generate the field -----------------
c						first clear it
      tt = second()
      do 20 k = 1, n3
      do 20 j = 1, n2
      do 20 i = 1, n1
         z(i,j,k) = zero
  20  continue
c						for each line ...
      do 40 ll = 1, L
c							generate 1-d process
         call fft1g( uni, np, ul, teq1d, kseed, 0, t1, sd, ierr )
c							determine its direction
         if( lrand ) then
             phi   = twopi*randu(0)
             theta = twopi*randu(0)
             u1    = sc*cos(theta)*cos(phi)
             u2    = sc*sin(theta)
             u3    = sc*cos(theta)*sin(phi)
         else
             u1    = sc*e1(ll)
             u2    = sc*e2(ll)
             u3    = sc*e3(ll)
         endif
c							decide which corner
         kx = 1
         ky = 1
         kz = 1
         if( u1 .lt. zero ) kx = 2
         if( u2 .lt. zero ) ky = 2
         if( u3 .lt. zero ) kz = 2
c							loop over field points
         zz = zstart(kz)
         do 30 k = 1, n3
            zz = zz + dz
            y  = ystart(ky)
         do 30 j = 1, n2
            y = y + dy
            x = xstart(kx)
         do 30 i = 1, n1
c						project line onto field point
            x = x + dx
            m = int( u1*x + u2*y +u3*zz ) + 1
            z(i,j,k) = z(i,j,k) + uni(m)
  30     continue
  40  continue
c						all done...
      ts = ts + (second() - tt)
      return
      end
