c  *********************************************************************
c  *                                                                   *
c  *                          function randu                           *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 2.0
c  Adapted from Press etal. "Numerical Recipes in Fortran", 2nd Ed.,
c  RAN1 program, pg 271
c  by Gordon A. Fenton, TUNS, Oct 14, 1996
c
c  PURPOSE  returns a pseudo-random number uniformly distributed on the
c           interval (0,1)
c
c  This routine uses a pair of multiplicative congruential generators to
c  produce a sequence of random numbers uniformly distributed on the interval
c  (0,1) which exceeds Park and Miller's ``Minimal Standard.'' The generator
c  has a period of about 2.3E+18 and uses a shuffling algorithm due to Bays
c  and Durham to remove low-order serial correlations. This routine is
c  about 2.0 times slower than the minimal routine `RAN0' presented in
c  Numerical Recipes. The endpoints of the interval (0,1) are excluded.
c
c  Arguments to the function are as follows;
c
c      jseed    integer seed. If jseed is zero, the next pseudo random number
c               in the sequence is returned. Otherwise if jseed is positive,
c               it is used to initialize the generator and the first value
c               in the sequence is returned. (input)
c
c  Notes:
c    1) if you always call this routine without initializing it (using
c       a positive seed) then you'll always get the same sequence of
c       pseudo-random numbers.
c
c  REVISION HISTORY:
c  2.0	implemented RAN2 from 2nd Edition of Numerical Recipes
c-----------------------------------------------------------------------------
      real function randu(jseed)
      integer jseed, idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      PARAMETER (IM1  = 2147483563,
     >           IM2  = 2147483399,
     >           AM   = 1./IM1,
     >           IMM1 = IM1-1,
     >           IA1  = 40014,
     >           IA2  = 40692,
     >           IQ1  = 53668,
     >           IQ2  = 52774,
     >           IR1  = 12211,
     >           IR2  = 3791,
     >           NTAB = 32,
     >           NDIV = 1+IMM1/NTAB,
     >           EPS  = 1.2e-7,
     >           RNMX = 1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum,idum2,ifirst
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/, ifirst/0/

      if( (jseed .gt. 0) .or. (ifirst .eq. 0) ) then
         ifirst = 1
         idum   = max0(jseed,1)
         idum2  = idum
         do 10 j = NTAB+8, 1, -1
            k    = idum/IQ1
            idum = IA1*(idum-k*IQ1) - k*IR1
            if( idum .lt. 0 ) idum  = idum + IM1
            if( j .le. NTAB ) iv(j) = idum
  10     continue
         iy = iv(1)
      endif

      k    = idum/IQ1
      idum = IA1*(idum-k*IQ1) - k*IR1
      if( idum .lt. 0 ) idum = idum + IM1
      k     = idum2/IQ2
      idum2 = IA2*(idum2-k*IQ2) - k*IR2
      if( idum2 .lt. 0 ) idum2 = idum2 + IM2
      j     = 1 + iy/NDIV
      iy    = iv(j) - idum2
      iv(j) = idum
      if( iy .lt. 1 ) iy = iy + IMM1
      randu = amin1( AM*iy, RNMX )

      return
      end
