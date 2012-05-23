c  ********************************************************************
c  *                                                                  *
c  *                          function dbeta                          *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 2.11
c  Written by Gordon A. Fenton, TUNS, Apr. 13, 1993
c  Latest Update: Feb 16, 2005
c
c  PURPOSE  returns the Beta distribution function for given x and parameters
c           alpha and beta.
c
c  This routine employs a continued fraction solution to compute the
c  cumulative probability associated with a variable following the Beta
c  distribution (also called the incomplete beta function),
c
c                G(a+b)     x
c     F(x) =   ---------  INT (1-t)**(b-1) * t**(a-1) dt,    0 <= x <= 1.0
c              G(a)*G(b)    0
c
c  where G(a) is the Gamma function ((a-1)! if `a' is an integer).
c  The algorithm follows that proposed by Press etal in "Numerical Recipes
c  in C", pg 178.
c  Arguments to the routine are as follows;
c
c     x   real value giving the position at which the cumulative distribution
c         is desired. (input)
c
c     a   the alpha parameter of the Beta distribution. (input)
c
c     b   the beta parameter of the Beta distribution. (input)
c
c  NOTES: 1) TOL is the acceptable relative error tolerance.
c         2) set IDEBUG = 1 if error messages should be written out.
c         3) IMAX is the maximum number of iterations allowed.
c
c  REVISION HISTORY:
c  2.1	changed `abs' to `dabs' and `float' to `dble' below (Dec 9/98)
c  2.2	added range on x in documentation above (Feb 16, 2005)
c---------------------------------------------------------------------------
      real*8 function dbeta(x,a,b)
      implicit real*8 (a-h,o-z)
      parameter (IDEBUG = 1, IMAX = 200, TOL = 1.d-10 )
      data zero/0.d0/, one/1.d0/, two/2.d0/
      exp(z)   = dexp(z)
      gamln(z) = dgamln(z)
      alog(z)  = dlog(z)
      abs(z)   = dabs(z)
      float(i) = dble(i)
c				check input data
      if( x .le. zero ) then
         dbeta = zero
         return
      elseif( x .ge. one ) then
         dbeta = one
         return
      endif

      ab = (one + a)/(a + b + two)
      bt = exp(gamln(a+b)-gamln(a)-gamln(b)+a*alog(x)+b*alog(one-x))
      if( x .lt. ab ) then
         xx = x
         aa = a
         bb = b
         st = zero
      else
         xx = one - x
         aa = b
         bb = a
         st = one
         bt = -bt
      endif
c				now use the continued fraction (Press etal.)
      am  = one
      bm  = one
      az  = one
      qab = aa + bb
      qap = aa + one
      qam = aa - one
      bz  = one - qab*xx/qap
      do 10 m = 1, IMAX
         em   = float(m)
         tem  = em + em
         d    = em*(bb-m)*xx/((qam+tem)*(aa+tem))
         ap   = az + d*am
         bp   = bz + d*bm
         d    =-(aa+em)*(qab+em)*xx/((aa+tem)*(qap+tem))
         app  = ap + d*az
         bpp  = bp + d*bz
         aold = az
         am   = ap/bpp
         bm   = bp/bpp
         az   = app/bpp
         bz   = one
         if( abs((az-aold)/az) .lt. TOL )go to 20
  10  continue
      if( IDEBUG .eq. 1 )
     >   write(0,'(a)')'DBETA couldn''t converge to a solution.'

  20  dbeta = st + bt*az/aa
      return
      end
