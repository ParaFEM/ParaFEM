c  ********************************************************************
c  *                                                                  *
c  *                       Subroutine cov1d                           *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  calculates the unbiased covariance function of a zero-mean 1-D
c           process
c
c  Calculates the covariance structure (unbiased) of a 1-D process.
c  The true mean of the process is assumed to be zero so that we can
c  divide by N rather than (N-1) and still be unbiased. This allows the
c  estimate of covariance based on a single observation (at the maximum
c  lag) to be determined without blowing up.
c  Arguments are described as follows;
c
c      c    real vector of size at least n which on ouput will contain the
c           desired covariance structure up to a lag of (n-1) starting from
c           a lag of zero. (output)
c
c      z    real vector of size at least n that contains the input process
c           for which the covariance structure is to be calculated. (input)
c
c      n    number of points in the process. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine cov1d( c, z, n )
      dimension c(*), z(*)
      data zero/0.0/

      do 20 i = 1, n
         s = zero
         do 10 j = i, n
            s = s + z(j-i+1)*z(j)
  10     continue
         c(i) = s/float(n-i+1)
  20  continue

      return
      end
