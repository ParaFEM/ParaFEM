c  *********************************************************************
c  *                                                                   *
c  *                        subroutine dcor1d                          *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Jun 24, 1997
c  Latest Update: Jul 2, 1997
c
c  PURPOSE  computes the sample mean, variance, and correlation function from
c           a discrete sequence of observations taken from a 1-D
c           random process.
c
c  DESCRIPTION
c  This routine computes the sample mean, variance, and correlation function
c  from a sequence of random observations. Provision has been made to call this
c  routine a number of times for different realizations of the observation
c  sequence and average the mean, variance, and correlation estimates over
c  the ensemble. In the event that the realizations consist of different
c  numbers of observations, the trailing correlations (from the maximum lag
c  of the current sequence to the maximum lag of the longest sequence) are
c  assumed to be zero.
c
c  ARGUMENTS
c  Arguments to this routine are as follows;
c
c	z	real vector of length at least [1 + (n-1)*ns] containing the
c		sequence of values to be examined. (input)
c
c	n	integer giving the number of elements in the sequence to be
c		processed. (input)
c
c      ns	integer giving the 'stride' between elements in z to be
c		processed. That is, the sequence to be examined is z(1),
c		z(1+ns), z(1+2*ns), ... z(1+(n-1)*ns), for a total of n
c		elements. (input)
c
c      zm	real value containing the sample mean. If nit > 1, then
c		the sample mean is added to the current value of zm, so
c		that the average sample average can be calculated on the
c		nit'th call to this routine. (input/output)
c
c    zmax	real value containing the maximum sample mean encountered
c		over the nit > 1 calls to this routine. zmax should not
c		be changed between calls. (input/output)
c
c    zmin	real value containing the minimum sample mean encountered
c		over the nit > 1 calls to this routine. zmin should not
c		be changed between calls. (input/output)
c
c      zv	real value containing the sample variance. If nit > 1, then
c		the sample variance is added to the current value of zv so
c		that the average sample variance can be calculated on the
c		nit'th call to this routine. (input/output)
c
c   zvmax	real value containing the maximum sample variance encountered
c		over the nit > 1 calls to this routine. zvmax should not
c		be changed between calls. (input/output)
c
c   zvmin	real value containing the minimum sample variance encountered
c		over the nit > 1 calls to this routine. zvmin should not
c		be changed between calls. (input/output)
c
c       r	real vector of length at least n, which will contain
c		the estimated correlation function. If nit > 1, then
c		the current estimate is added to the current value of
c		r so that the average estimate can be computed on the
c		nit'th call to this routine. Note that r(1) = 1, corresponding
c		to the correlation at zero lag, r(2) contains the correlation
c		coefficient at lag 1 (index space), etc. (input/output)
c
c    rmax	real vector of length at least n which will contain the
c		maximum correlation encountered over the nit > 1 calls
c		to this routine at each lag. rmax should not be changed
c		between calls to this routine and has no meaning if nit = 1.
c		(input/output)
c
c    rmin	real vector of length at least n which will contain the
c		minimum correlation encountered over the nit > 1 calls
c		to this routine at each lag. rmin should not be changed
c		between calls to this routine and has no meaning if nit = 1.
c		(input/output)
c
c      ni	integer vector of length at least nmax, containing the
c		number of correlation samples summed into the average (r)
c		for each lag. This is used to correctly compute the average
c		correlation at each lag on the nit'th call to this routine.
c		This vector is only used when nit > 1 and only differs from
c		nit at larger lags when samples are not of the same length.
c		(input/output)
c
c   lbias	logical flag which is true if the biased (but reduced
c		variance) estimator should be used. This implies that
c		the sum of deviation products is divided by n rather
c		than by (n - lag - 1) as required by the unbiased
c		estimator. The unbiased estimator is used if lbias
c		is false, in which case, correlations from lag zero up to
c		lag (n-2) are evaluated. Most investigators now use the
c		biased estimator. (input)
c
c     nen	entry number - this is the number of times that this routine
c		has been called. When nen = nit, then the average correlation
c		is properly normalized prior to returning. (input)
c
c     nit	the total number of realizations over which the calculations
c		are to be performed. See nen for details. (input)
c
c    nmax	integer containing the maximum length of all samples, z,
c		processed in the sequence of nit calls to this routine.
c		If nit = 1, this argument is unused. If all samples are
c		of the same length, then set nmax = n. (input)
c
c  REVISION HISTORY:
c  1.1	returns max/min correlations rather than correlation var. (Jul 2/97)
c---------------------------------------------------------------------------
      subroutine dcor1d(z,n,ns,zm,zmax,zmin,zv,zvmax,zvmin,
     >                  r,rmax,rmin,ni,lbias,nen,nit,nmax)
      implicit real*8 (a-h,o-z)
      dimension z(*), r(*), rmax(*), rmin(*)
      integer ni(*)
      logical lbias
      data zero/0.d0/, one/1.d0/, two/2.d0/
      float(i) = dble(i)
c					estimate the mean and variance
      if( n .le. 1 ) then
         if( n .eq. 1 ) zm = z(1)
         return
      endif
      tm = z(1)
      vm = z(1)*z(1)
      do 10 i = 2, n
         tm = tm + z(1+(i-1)*ns)
         vm = vm + z(1+(i-1)*ns)*z(1+(i-1)*ns)
  10  continue
      am = tm/float(n)
      sm = vm - am*tm
      if( sm .lt. zero ) sm = zero
      if( lbias ) then
         vm = sm/float(n)
      else
         vm = sm/float(n-1)
      endif
      if( nen .eq. 1 ) then
         zm    = am
         zmax  = am
         zmin  = am
         zv    = vm
         zvmax = vm
         zvmin = vm
      else
         zm = zm + am
         if( am .gt. zmax ) zmax = am
         if( am .lt. zmin ) zmin = am
         zv = zv + vm
         if( vm .gt. zvmax ) zvmax = vm
         if( vm .lt. zvmin ) zvmin = vm
      endif
c					normalize the mean and variance?
      if( nen .eq. nit ) then
         zm = zm/float(nit)
         zv = zv/float(nit)
      endif
c					initialize for nit > 1 case
      if( (nit .gt. 1) .and. (nen .eq. 1) ) then
         r(1)    = one
         rmax(1) = one
         rmin(1) = one
         ni(1)   = 1
         do 20 i = 2, nmax
            ni(i)   =  0
            r(i)    =  zero
            rmax(i) = -two
            rmin(i) =  two
  20     continue
      elseif( nit .eq. 1 ) then
         r(1) = one
      endif
c					estimate the correlation function
      do 40 i = 2, n			! lag is (i-1)
         v = zero
         do 30 j = i, n
            d1 = z(1+(j-1)*ns) - am
            d2 = z(1+(j-i)*ns) - am
            v  = v + d1*d2
  30     continue
         if( lbias ) then
            v = v/sm
         else
            if( i .lt. n ) then
               v = v/(vm*float(n-i))
            else
               v = zero
            endif
         endif
         if( nit .eq. 1 ) then
            r(i) = v
         else
            ni(i) = ni(i) + 1
            r(i)  = r(i)  + v
            if( v .gt. rmax(i) ) rmax(i) = v
            if( v .lt. rmin(i) ) rmin(i) = v
         endif
  40  continue
c					normalize ensemble averages
      if( (nen .eq. nit) .and. (nit .gt. 1) ) then
         do 50 i = 2, nmax
            if( ni(i) .gt. 0 ) r(i) = r(i)/float(ni(i))
  50     continue
      endif

      return
      end
