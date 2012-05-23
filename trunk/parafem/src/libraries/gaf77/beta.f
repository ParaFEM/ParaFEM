c  ********************************************************************
c  *                                                                  *
c  *                          function beta                           *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 2.0
c  Written by Gordon A. Fenton, TUNS, Apr. 13, 1993
c
c  PURPOSE  returns the Beta distribution function for given x and parameters
c           alpha and beta.
c
c  This routine employs a continued fraction solution to compute the
c  cumulative probability associated with a variable following the Beta
c  distribution (also called the incomplete beta function),
c
c                G(a+b)     x
c     F(x) =   ---------  INT (1-t)**(b-1) * t**(a-1) dt
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
c---------------------------------------------------------------------------
      real function beta(x,a,b)
      parameter (IDEBUG = 1, IMAX = 200, TOL = 1.e-6 )
      data zero/0.0/, one/1.0/, two/2.0/
   1  format(a)
c				check input data
      if( x .le. zero ) then
         beta = zero
         return
      elseif( x .ge. one ) then
         beta = one
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
     >   write(0,1)'BETA couldn''t converge to a solution.'

  20  beta = st + bt*az/aa
      return
      end
