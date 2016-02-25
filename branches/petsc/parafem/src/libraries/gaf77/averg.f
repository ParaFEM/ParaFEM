c  ******************************************************************
c  *                                                                *
c  *                     Function Averg                             *
c  *                                                                *
c  ******************************************************************
c  Single Precision Version 1.01
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
      real function averg( X, n )
      dimension X(*)
c
      is    = 2
      averg = X(1)
      if( n .le. 1 ) return
      if( mod(n,2) .eq. 0 ) then
         is = 3
         averg = X(1) + X(2)
      endif
      do 10 i = is, n, 2
         averg = averg + X(i) + X(i+1)
  10  continue
c
      averg = averg/float(n)
c
      return
      end
