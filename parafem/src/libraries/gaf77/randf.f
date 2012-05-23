c  *********************************************************************
c  *                                                                   *
c  *                          function randf                           *
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
c  This routine uses a simple multiplicative congruential generator to
c  produce a sequence of random numbers uniformly distributed on the interval
c  (0,1) which satisfies Park and Miller's ``Minimal Standard.'' The generator
c  has a period of about 2.1E+09 and uses a shuffling algorithm due to Bays
c  and Durham to remove low-order serial correlations. This routine is only
c  about 1.3 times slower than the minimal routine `RAN0' presented in
c  Numerical Recipes.
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
c  1.01	added comment re periodicity of the generator (old RAN1)
c  2.0	implemented RAN1 from 2nd Edition of Numerical Recipes
c-----------------------------------------------------------------------------
      real function randf(jseed)
      integer jseed, idum, IA, IM, IQ, IR, NTAB, NDIV
      real AM, EPS, RNMX
      parameter (IA   = 16807,
     >           IM   = 2147483647,
     >           AM   = 1./IM,
     >           IQ   = 127773,
     >           IR   = 2836,
     >           NTAB = 32,
     >           NDIV = 1+(IM-1)/NTAB,
     >           EPS  = 1.2e-7,
     >           RNMX = 1.-EPS)
      INTEGER j, k, iv(NTAB), iy
      save iv,iy, idum
      data iv /NTAB*0/, iy /0/

      if( (jseed .gt. 0) .or. (iy .eq. 0) ) then
         idum = max0(jseed,1)
         do 10 j = NTAB+8, 1, -1
            k    = idum/IQ
            idum = IA*(idum-k*IQ) - IR*k
            if( idum .lt. 0 ) idum  = idum + IM
            if( j .le. NTAB ) iv(j) = idum
  10     continue
         iy = iv(1)
      endif

      k     = idum/IQ
      idum  = IA*(idum-k*IQ) - IR*k
      if( idum .lt. 0 ) idum = idum + IM
      j     = 1 + iy/NDIV
      iy    = iv(j)
      iv(j) = idum
      randf = amin1( AM*iy, RNMX )

      return
      end
