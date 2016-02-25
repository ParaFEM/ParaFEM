c  *********************************************************************
c  *                                                                   *
c  *                        subroutine dvfn1d                          *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Jun 25, 1997
c  Latest Update: Mar 6, 1998
c
c  PURPOSE  computes the sample variance function from a discrete sequence of
c           observations taken from a 1-D random process.
c
c  DESCRIPTION
c  This routine computes the sample variance function from a sequence of
c  random observations. The sample variance function gives the variance
c  reduction of the 'local' average of the sample as a function of the 'size'
c  of the local average. The estimate is obtained by estimating the
c  variance of a moving average of size i applied to the original sequence.
c  For example, if Z(I,j) is defined as the local average
c
c                 1    j+I-1
c	Z(I,j) = ---   sum   Z(k)	for I = 1, 2, ..., n   (averaging size)
c                 I    k=j		and j = 1, 2, ..., n-I+1
c
c  using the raw process Z(k), k = 1, 2, ..., n,
c  then the variance of the average of size I is estimated as
c
c                  1    n-I+1
c	 V(I) = ------- sum   (Z(I,j) - m)^2
c                n-I+1  j=1
c
c  where m is the globally estimated mean (m = Z(n,1), which implies that
c  V(n) is always 0). The variance function, VF(i), is obtained by dividing
c  the local average variance, V(i), by the sample variance, V(1), so that
c  VF(1) = 1 always. The variance function estimate is biased (but
c  asymptotically unbiased as the process length goes to infinity).
c
c  Provision has been made to call this routine a number of times for
c  different realizations of the observation sequence and average the
c  variance function estimates over the ensemble.
c
c  In the event that the realizations consist of different numbers of
c  observations, the trailing variance function values (from the end of a
c  given sequence to the end of the longest sequence) are assumed to be zero.
c  In this case, there is some question as to how representative the average
c  variance function will be towards the tail, since VF(n) is guaranteed to
c  equal zero and n may be typically less than nmax. This means that average
c  variance function values towards the tail may be somewhat lower than they
c  might be if all samples were of equal length.
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
c      zm	real value containing the average of the sequence z. (output)
c
c      zv	real value containng the estimated variance of the sequence
c		z (biased or unbiased in the classical sense, see lbias).
c		(output)
c
c      vf	real vector of length at least n, which will contain
c		the estimated variance function. If nit > 1, then
c		the current estimate is added to the current value of
c		vf so that the average estimate can be computed on the
c		nit'th call to this routine.  (input/output)
c
c    vmax	real vector of length at least n, which will contain the
c		maximum variance function value for each sample size
c		observed over the nit calls to this routine. This vector
c		is only useful if nit > 1. The vector vmax should not be
c		modified between calls to this routine. (input/output)
c
c    vmin	real vector of length at least n, which will contain the
c		minimum variance function value for each sample size
c		observed over the nit calls to this routine. This vector
c		is only useful if nit > 1. The vector vmin should not be
c		modified between calls to this routine. (input/output)
c
c      ni	integer vector of length at least nmax, containing the number
c		of var func samples summed into the average (vf) for each
c		sample size. This is used to correctly compute the average
c		var func at each size on the nit'th call to this routine.
c		This vector is only used when nit > 1 and only differs from
c		nit at larger lags when samples are not of the same length.
c		(input/output)
c
c     tmp	temporary real vector of length at least n used as workspace
c		(to store the sequence of local averages of size i = 2, 3,
c		 .. n)
c
c   lbias	logical flag which is true if the process variance estimate
c		is to be biased (division by n) rather than unbiased (division
c		by n-1). (input)
c
c     nen	entry number - this is the number of times that this routine
c		has been called. When nen = nit, then the estimated
c		statistics are properly normalized prior to returning.
c		(input)
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
c  1.1	return max/min rather than variance of estimate. (Jul 2/97)
c  1.2	fixed nmax endpoints (Mar 6/98)
c---------------------------------------------------------------------------
      subroutine dvfn1d(z,n,ns,zm,zv,vf,vmax,vmin,ni,tmp,
     >                  lbias,nen,nit,nmax)
      implicit real*8 (a-h,o-z)
      dimension z(*), vf(*), vmax(*), vmin(*), tmp(*)
      integer ni(*)
      logical lbias
      data zero/0.d0/, one/1.d0/
      float(i) = dble(i)

   1  format(a,$)
   2  format(i5,$)
   3  format()
c					estimate the mean and variance
      if( n .le. 1 ) return
      tm = z(1)
      vm = z(1)*z(1)
      tmp(1) = tm
      do 10 i = 2, n
         tm = tm + z(1+(i-1)*ns)
         vm = vm + z(1+(i-1)*ns) * z(1+(i-1)*ns)
         tmp(i) = tm
  10  continue
      zm = tm/float(n)
      sm = vm - zm*tm
      if( sm .lt. zero ) sm = zero
      if( lbias ) then
         zv = sm/float(n)
      else
         zv = sm/float(n-1)
      endif
c					initialize for nit > 1 case
      if( (nit .gt. 1) .and. (nen .eq. 1) ) then
         ni(1)   = 1
         vf(1)   = one
         vmax(1) = one
         vmin(1) = one
         do 20 i = 2, nmax-1
            ni(i)   = 0
            vf(i)   = zero
            vmax(i) = zero
            vmin(i) = one
  20     continue
         vf(nmax)   = zero
         vmax(nmax) = zero
         vmin(nmax) = zero
      endif
c					estimate variance function values
c					for 1-step averaging (ie. observations)
      if( nit .eq. 1 ) then
         vf(1) = one
         vf(n) = zero
      endif
c					for 2, 3, ... n-1 step averaging
      do 40 i = 2, n-1
c						local average variances
         al = (tmp(i)/float(i)) - zm
         vl = al*al
         do 30 j = 2, n-i+1
            al = al + (z(1+(j+i-2)*ns)-z(1+(j-2)*ns))/float(i)
            vl = vl + al*al
  30     continue
         if( lbias ) then
            vl = vl/(zv*float(n-i+1))
         else
            vl = vl/(zv*float(n-i))
         endif
c						store the results
         if( nit .eq. 1 ) then
            vf(i)  = vl
         else
            vf(i) = vf(i) + vl
            ni(i) = ni(i) + 1
            if( vl .gt. vmax(i) ) vmax(i) = vl
            if( vl .lt. vmin(i) ) vmin(i) = vl
         endif
  40  continue
c					normalize var func averages
      if( (nen .eq. nit) .and. (nit .gt. 1) ) then
         do 50 i = 2, nmax - 1
            if( ni(i) .gt. 0 ) vf(i) = vf(i)/float(ni(i))
  50     continue
      endif
c					all done
      return
      end
