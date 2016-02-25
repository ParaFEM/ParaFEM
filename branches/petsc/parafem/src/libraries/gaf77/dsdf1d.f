c  *********************************************************************
c  *                                                                   *
c  *                         subroutine dsdf1d                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.5
c  Written by Gordon A. Fenton, TUNS, Jun 25, 1997
c  Latest Update: Jun 19, 2006
c
c  PURPOSE  estimate the spectral density function of a real 1-D process
c
c  DESCRIPTION
c  This routine takes as input a sequence of observations (a record)
c  drawn from a realization of a 1-D random process and estimates its
c  spectral density function. This is done by first computing its
c  periodogram (using a Fast Fourier transform of the zero-meaned
c  process) and then smoothing the periodogram via one of;
c
c	1) averaging over an ensemble (involving multiple calls to
c	   this routine). This is equivalent to applying a Bartlett
c	   window to the entire set of records joined end to end.
c
c	2) applying a centered Daniell smoothing window (simple rectangular).
c	   SDF values near zero are smoothed in a special fashion;
c		s(1) = p(1)	(ie, first ordinate is not smoothed at all)
c		s(2) = [p(1)+p(2)+p(3)]/3
c		.
c		.
c		s(m) = [p(1) + ... + p(2*m-1)]/(2*m-1)
c	   and the remainder are smoothed over nsmuth = 2*m+1 points (m
c	   points on either side). The last m frequency values are smoothed
c	   by double counting values `reflected' around 1+n/2 from beyond
c	   1+n/2.
c
c  ARGUMENTS
c	z	real vector of length at least nn containing the data
c		to be analyzed. If nn is not a power of 2, then nn is
c		increased to the next higher power of 2 and the corresponding
c		number of zeroes are added to the end of z. In this
c		case a warning message is written to standard output.
c		On output, z will have its average removed if lzero is
c		true. (input/output)
c
c      zi	temporary real vector of length at least nn which is
c		used to store the imaginary Fourier coefficients.
c
c      nn	the number of elements in the vector z. In this version,
c		nn must be a power of 2. If it is not, it is increased
c		to the next higher power of 2 and zeroes are added to
c		the end of the vector z. (input/output)
c
c	L	the number of zero pads added to the end of nn to bring
c		it up to a power of 2 by the calling routine. If nn is
c		NOT a power of 2, L is ignored (it is assumed that there
c		are no added zeroes in z) and nn is brought up to being
c		a power of 2 by adding zeroes. The number of zeroes added
c		herein is then returned in L. (input/output)
c
c      dt	real value giving the increment (time or space) between
c		elements of z. (input)
c
c	s	real vector of length at least 1+nn/2 which on output will
c		contain the estimated spectral density function ordinates.
c		If nit > 1, then smoothing is assumed to be taking place
c		by averaging over the ensemble of nit records. Thus, the
c		current SDF estimates are added to the vector s. When
c		nen = nit, the vector s is normalized by dividing by nit
c		to give the average SDF estimates. The SDF is evaluated at the
c		Fourier frequencies w_j = 2*pi*(j-1)/(nn*dt), where dt
c		is the time increment between sequential elements of z, for
c		j = 1, 2, ..., 1+nn/2.
c		If nit = 1, and nsmuth is non-zero, then s will contain the
c		SDF estimates smoothed by averaging over a rectangular
c		window of width nsmuth (see nsmuth and discussion above).
c		(input/output)
c
c    smax	real vector of length at least 1+nn/2 which on output will
c		contain the maximum spectral density value observed over
c		the sequence of nit calls to this routine at each Fourier
c		frequency. This argument is only set if nit > 1. It should
c		not be changed between calls to this routine. (input/output)
c
c
c    smin	real vector of length at least 1+nn/2 which on output will
c		contain the minimum spectral density value observed over
c		the sequence of nit calls to this routine at each Fourier
c		frequency. This argument is only set if nit > 1. It should
c		not be changed between calls to this routine. (input/output)
c
c	p	real vector of length at least 1+nn/2 used to contain the
c		periodogram estimates. These are the raw un-smoothed
c		estimates for this particular input vector z. (output)
c
c   lzero	logical flag which is true if the input process is to be
c		zero-averaged, by subtracting the average from each element,
c		prior to the spectral analysis. This ensures that the power
c		at zero frequency will be zero. If the input process is not
c		zero-averaged, then the significant power present at zero
c		frequency tends to `leak' into adjacent spectral ordinates,
c		particularly after smoothing (nit = 1, nsmuth != 0 case),
c		corrupting their values. However, it may not always be
c		desirable to zero-average the input sequence, particularly
c		for processes which are known to have a pole at the frequency
c		origin (such as fractional Gaussian noise). If the process is
c		NOT zero-averaged, and 'zeroes' are required to pad the
c		process out to a power of 2 in length, then the trailing
c		L elements of z are set to the process average, rather than
c		to zeroes. (input)
c
c   ltape	logical flag which is true if data is to be tapered at
c		both ends using a 'Bingham' cosine taper on the first
c		and last 10% of the record. This is said to slightly
c		improve the shape of the resulting spectral window at
c		the expense of slightly increased effective bandwidth
c		(reduced resolution) and a slight reduction in statistical
c		accuracy. (input)
c
c  nsmuth	integer indicating whether smoothing of the spectral
c		estimates is to take place. In general, this is recommended
c		since raw periodogram estimates are notoriously badly
c		behaved and non-consistent (variance of the estimate does
c		not decrease with increasing nn). Smoothing the periodogram
c		produces a consistent estimate of the SDF as well as a more
c		reasonable appearance overall. The smoothing window used is
c		the Daniell window (simple rectangular). See Priestley,
c		"Spectral Analysis and Time Series", Chapter 6, for details.
c		The parameter nsmuth is interpreted as follows; if nsmuth
c		  = -1	then this routine derives a reasonable value based
c			on an approximation to the spectral bandwidth,
c		  = 0	then no smoothing takes place,
c		  > 0	then smoothing is performed using a centered local
c			average over the specified number of Fourier
c			frequencies. The specified value must be odd, so
c			that the average is performed over the (nsmuth-1)/2
c			values before and after the current value.
c		Note that this parameter is ignored if nit > 1 (in which
c		case `smoothing' is performed by averaging over the ensemble.
c		(input)
c
c     nen	entry number - this is the number of times that this routine
c		has been called. When nen = nit, then the SDF average
c		is properly normalized prior to returning.
c		(input)
c
c     nit	the total number of realizations over which the calculations
c		are to be performed. See nen for details. If nit > 1, then
c		the SDF is assumed to be being smoothed by averaging over
c		an ensemble. Otherwise, the SDF is smoothed via a window
c		if nsmuth is non-zero, see nsmuth. (input)
c
c  REVISION HISTORY:
c  1.1	changed SDF variance in favour of max/min reporting. (Jul 3/97)
c  1.2	corrected G(0) - no longer mult by 2 (for one-sided) at origin
c	(Mar 12/98)
c  1.3	corrected single realization smoothing (Apr 28/03)
c  1.4	since forward FFT used, which no longer divides by N, correct
c	estimate below to no longer multiply by N*N (Apr 5/04)
c  1.5	forward FFT does divide by N, so multiply by N*N below (!) (Jun 19/06)
c-------------------------------------------------------------------------
      subroutine dsdf1d(z,zi,nn,L,dt,s,smax,smin,p,
     >                  lzero,ltape,nsmuth,nen,nit)
      implicit real*8 (a-h,o-z)
      dimension z(*), zi(*), s(*), smax(*), smin(*), p(*)
      logical ltape, lzero
      data zero/0.d0/, half/0.5d0/, one/1.d0/, two/2.d0/
      data    pi/3.1415926535897932384d0/
c     data twopi/6.2831853071795864769d0/
      float(i) = dble(i)
      cos(q)   = dcos(q)
      sqrt(q)  = dsqrt(q)
c     int(q)   = idint(q)

   1  format(5a)
c						deduce original sample length
c						nn = total length,
c						n  = non-padded length
      n = nn
      if( L .gt. 0 ) n = nn - L
         
c						compute sample mean of z
      if( nit .le. 1 ) then
         zz    = z(1)				! when nit = 1, then we will
         zv    = z(1)*z(1)			! be smoothing the SDF.
         z1    = zero				! Here we estimate the size
         do 10 i = 2, n				! of the smoothing window by
            zz    = zz + z(i)			! the rate of fall of the COV
            zv    = zv + z(i)*z(i)		! function... and then use an
            z1    = z1 + z(i)*z(i-1)		! approximation by Priestley
  10     continue				! which assumes an AR(1)
         zm = zz/float(n)			! process.
         zv = zv - zm*zz
         if( zv .le. zero ) then
            write(6,1)'Warning: SDF of a zero-variance process is 0.'
            return
         endif
         z1 = z1 - zm*(zz + zm - z(1) - z(n))
         aa = z1/zv				! estimated AR(1) coefficient
      else
         zz    = z(1)				! when nit > 1, we will smooth
         do 20 i = 2, n				! by averaging over the
            zz    = zz + z(i)			! ensemble.
  20     continue
         zm = zz/float(n)
      endif
c						zero mean z, zero zi
      if( lzero ) then
         zpad = zero
         do 30 i = 1, n
            z(i)  = z(i) - zm
  30     continue
      else
         zpad = zm
      endif
c						check n = 2**m
      m = ntwom(nn)
      if( 2**m .lt. nn ) then			! it isn't; adjust nn and L
         n = nn
         m  = m + 1
         nn = 2**m
         L  = nn - n
         do 40 i = n+1, n+L
            z(i)  = zpad
  40     continue
      endif
c						taper the data?
      f = float(nn*nn)/float(n)
      if( ltape ) then
         ns = n/10
         r  = float(n - 2*(ns + 1))
         do 50 i = 1, ns + 1
            an = pi*float(i-1)/float(ns)
            dc = cos(an)
            d1 = half*(one - dc)
            r = r + two*d1
            z(i) = zpad + (z(i)-zpad)*d1
            z(n-i+1) = zpad + (z(n-i+1)-zpad)*d1
  50     continue
         f = f*float(n)/r
      endif
c						zero imaginary portion
      do 60 i = 1, nn
         zi(i) = zero
  60  continue
c						compute periodogram
      call dfft1d(z,zi,m,.false.)

      f    = f*dt/pi				! one-sided SDF
      nby2 = nn/2
      nb2p = 1 + nby2
      do 70 i = 1, nb2p
         p(i) = f*(z(i)*z(i) + zi(i)*zi(i))
  70  continue
c						correct zero frequency value
      p(1) = p(1)/two				! correct for one-sided mult

c						compute SDF estimate
      if( nsmuth .eq. 0 .and. nit .le. 1 ) then
         do 80 i = 1, nb2p
            s(i)  = p(i)
  80     continue
         return
      endif
c						ensemble averaging?
      if( nit .gt. 1 ) then
         if( nen .eq. 1 ) then
            do 90 i = 1, nb2p
               s(i)    = p(i)
               smax(i) = p(i)
               smin(i) = p(i)
  90        continue
         else
            do 100 i = 1, nb2p
               s(i) = s(i) + p(i)
               if( p(i) .gt. smax(i) ) smax(i) = p(i)
               if( p(i) .lt. smin(i) ) smin(i) = p(i)
 100        continue
         endif
         if( nen .eq. nit ) then
            do 110 i = 1, nb2p
               s(i) = s(i)/float(nit)
 110        continue
         endif
         return
      endif
c						apply smoothing window (nit=1)
      if( nsmuth .eq. -1 ) then			! determine good value
         if( aa .le. zero ) then
            nsmuth = 5
         else
            bh = two*(one - aa)/sqrt(aa)	! see Priestley, pg 517, 537
            nsmuth = 0.357d0*(bh*float(nn))**(0.8d0)
         endif
         nmx = int( half + sqrt( float(nn) ) )	! This is an arbitrary max
         if( nsmuth .gt. nmx ) nsmuth = nmx	! to avoid oversmoothing.
         if( nsmuth .le. 1   ) nsmuth = 3	! Always do some smoothing
      endif					! if we get to this point.
      m =  nsmuth/2
      if( 2*m .eq. nsmuth ) nsmuth = nsmuth + 1
      as = one/float(nsmuth)
c						for frequencies k = 1 to m
      s(1) = p(1)
      a    = zero
      do 112 k = 2, m
         a = a + p(2*k-2) + p(2*k-1)
         s(k) = a/float(2*k-1)
 112  continue
c						for frequency k = 1 + N/2
      a2 = p(nb2p)
      do 120 i = 2, m+1
         a2 = a2 + two*p(nb2p-i+1)
 120  continue
      s(nb2p) = as*a2
c						for k = N/2, N/2-1, ... N/2-m+1
      do 150 k = 2, m
         a2 = p(nb2p)
         do 130 i = 2, 2-k+m
            a2 = a2 + p(nb2p-i+1)
 130     continue
         do 140 i = 2, k+m
            a2 = a2 + p(nb2p-i+1)
 140     continue
         s(nb2p-k+1) = as*a2
 150  continue
c						for k = m+1, ..., N/2+1-m
      a = zero
      do 160 i = 1, 2*m+1
         a = a + p(i)
 160  continue
      s(m+1) = as*a
      do 170 k = m+2, nb2p-m
         a = a + p(k+m) - p(k-m-1)
         s(k) = as*a
 170  continue
c						all done
      return
      end
