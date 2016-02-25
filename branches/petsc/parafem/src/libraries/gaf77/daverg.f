c  ******************************************************************
c  *                                                                *
c  *                     Function Daverg                            *
c  *                                                                *
c  ******************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, Jan. 20, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  returns the average of a series of values X(i), i = 1,2, ..., n
c
c  Arguments to the routine are as follows;
c
c     X    real vector of length at least N containing the values to be
c          averaged. (input)
c
c     n    number of values in X to be averaged. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c------------------------------------------------------------------------------
      real*8 function daverg( X, n )
      implicit real*8 (a-h,o-z)
      dimension X(*)
      float(i) = dble(i)
c
      daverg = X(1)
      if( n .le. 1 ) return
      do 10 i = 2, n
         daverg = daverg + X(i)
  10  continue
c
      daverg = daverg/float(n)
c
      return
      end
