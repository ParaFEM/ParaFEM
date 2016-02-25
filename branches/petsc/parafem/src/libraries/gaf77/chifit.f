c  *********************************************************************
c  *                                                                   *
c  *                         function chifit                           *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, 1992
c  Latest Update: Aug 27, 2000
c
c  PURPOSE   to estimate the Chi-squared goodness-of-fit critical significance
c            level (p-value)
c
c  This routine compares an observed histogram against a predicted histogram
c  (based on the assumed distribution), computes a Chi-squared statistic and
c  then the significance level, or p-value, corresponding to the statistic.
c  The function returns the critical significance level (p-value). Large
c  values of the critical significance level imply a better fit to the assumed
c  distribution. The hypothesis that the data comes from the assumed
c  distribution is not rejected if the target significance level is less than
c  the critical significance (p-value). Note that histogram buckets where the
c  predicted number of occurrences is less than PMIN are ignored.
c
c  Arguments to this routine are as follows;
c
c    xb    real vector of length at least NBKT+1 containing the values of the
c          assumed cumulative distribution at the boundaries of each bucket,
c          starting from the left edge of the first bucket and proceeding to
c          the right edge of the last bucket. (input)
c
c    freq  real vector of length at least NBKT containing the observed
c          frequencies in each bucket. (input)
c
c    npt   integer giving the total number of samples. (input)
c
c    nbkt  integer giving the number of buckets in the histogram. (input)
c
c    pmin  real value denoting the threshold of the predicted frequency
c          below which the bucket is ignored. (input)
c
c    nest  the number of parameters of the assumed distribution which are
c          estimated from the original data. (input)
c
c    chisq real value containing the computed chi-squared statistic. (output)
c
c  REVISION HISTORY:
c  1.2	revised above description to reflect `p-value' (Aug 27, 2000)
c------------------------------------------------------------------------------
      real function chifit(xb,freq,npt,nbkt,pmin,nest,chisq)
      dimension xb(*), freq(*)
      real*8 vpt, vmin, chi, epr, err, dble, zero, one, two
      real*8 dgmdst
      data zero/0.d0/, one/1.d0/, two/2.d0/

   1  format(a)
   2  format(a,i3,a)
   3  format(5f13.2)
   4  format(a,e12.5)

      vmin = dble(pmin)
      vpt  = dble(npt)
      nb0  = 0
      chi  = zero
      do 120 k = 1, nbkt
         epr = vpt*(dble(xb(k+1)) - dble(xb(k)))
         if( epr .gt. vmin ) then
            nb0 = nb0 + 1
            err = epr - dble(freq(k))
            chi = chi + err*err/epr
         endif
 120  continue

      chisq  = chi
      vpt    = dble(nb0 - nest - 1)
      if( vpt .lt. one ) then
         chifit = 0.0
      else
         chifit = sngl( one - dgmdst(chi,vpt/two,two) )
      endif

      return
      end
