c  *********************************************************************
c  *                                                                   *
c  *                       function dgamln                             *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Apr. 13, 1993
c
c  PURPOSE  returns the natural logarithm of the gamma function of the
c           supplied argument, z >= 0.
c
c  This routine calculates the definite integral
c
c                  inf
c    gamma(z) = int  t**(z-1)*exp(-t) dt
c                  0
c
c  which is the Gamma function, for argument z >= 0, and returns ln(gamma(z)).
c  The method used is an approximation derived by Lanczos as presented by
c  Press etal. in "Numerical Recipes in C", page 167. This algorithm supposedly
c  has a error of less than 2.E-10. See also GAMMA.f which has a relative
c  error of less than 1.E-14 and is based on a series expansion and can handle
c  negative arguments.
c  Arguments are as follows;
c
c     z    non-negative real value for which the gamma function is desired.
c          (input)
c
c  NOTES: 1) (MXFAC-1)! is the maximum integer factorial we can find without
c            overflowing
c         2) the integer factorial values are stored in fac(i) for later
c            (efficient) use. The first 7 are set explicitly, ie fac(1-7) = 
c            0!, 1!, ..., 6!.
c         3) set IDEBUG to 1 if you want overflow errors reported.
c--------------------------------------------------------------------------
      real*8 function dgamln(z)
      implicit real*8 (a-h,o-z)
      parameter (IDEBUG = 1, MXFAC = 171)
      dimension cof(6), fac(MXFAC)
      data cof/76.180091730d0,-86.50532033d0,24.01409822d0,
     >         -1.231739516d0,  .120858003d-2,-.536382d-5/
      data stp/2.50662827465d0/, eps/0.222045d-15/
      data zero,half,one,fpf/0.d0,0.5d0,1.0d0,5.5d0/
      data nfac/7/
      data (fac(i),i=1,7)/1.d0,1.d0,2.d0,6.d0,24.d0,120.d0,720.d0/
      float(i) = dble(i)
      int(x)   = idint(x)
      alog(x)  = dlog(x)

   1  format(a,e13.6,a)

      if( z .lt. zero ) then
         if( IDEBUG .eq. 1 )
     >      write(0,1)'Negative argument (',z,') in DGAMLN.'
         dgamln = zero
         return
      endif
c					is z an integer?
      q = z*(one + eps)
      n = int( q )
      f = z - float(n)
      if( f .lt. q*eps .and. n .le. MXFAC ) then
         if( n .gt. nfac ) then
            do 10 i = nfac, n-1
               fac(i+1) = float(i)*fac(i)
  10        continue
            nfac = n
         endif
         dgamln = alog(fac(n))
         return
      endif
c					z is non-integer (or is a big integer)
      if( z .lt. one ) then
         x   = z
         az  = alog(z)
      else
         x  = z - one
         az = zero
      endif
      tmp = x + fpf
      tmp = (x+half)*alog(tmp) - tmp
      ser = one
      do 20 j = 1, 6
         x   = x + one
         ser = ser + cof(j)/x
  20  continue
      dgamln = tmp + alog(stp*ser) - az
      return
      end
