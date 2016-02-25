c  ********************************************************************
c  *                                                                  *
c  *                         function dgmdst                          *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 2.0
c  Written by Press etal., Numerical Recipes, with some modifications by
c  Gordon A. Fenton, TUNS, 1994
c  Latest Update: May 9, 1994
c
c  PURPOSE  returns the Gamma distribution function for given x and parameters
c           alpha and beta.
c
c  This routine employs series or partial fractions solutions to compute the
c  cumulative probability associated with a variable following the Gamma
c  distribution,
c
c                   1       x
c     F(x) =   ---------  INT exp{-t/b}*t**(a-1) dt
c              G(a)*b**a    0
c
c  where G(.) is the Gamma function (G(a) = (a-1)! if `a' is an integer).
c  When b = 1, this is the incomplete Gamma function. When b = 2 and a = k/2,
c  this is the Chi-squared distribution with k degrees of freedom.
c  This routine is derived from Press etal. in "Numerical Recipes",
c  2nd Edition, pg. 218.
c  Arguments to the routine are as follows;
c
c     x   real value giving the position at which the cumulative distribution
c         is desired. (input)
c
c     a   the alpha parameter of the Gamma distribution (a > 0). (input)
c
c     b   the beta parameter of the Gamma distribution (b > 0). (input)
c---------------------------------------------------------------------------
      real*8 function dgmdst( x, a, b )
      implicit real*8 (a-h,o-z)
      parameter (ITMAX = 400)
      data zero/0.d0/, one/1.d0/, two/2.d0/, tol/1.d-12/, small/1.d-30/
      gamln(z) = dgamln(z)
      alog(z)  = dlog(z)
      exp(z)   = dexp(z)
      float(i) = dble(i)
      abs(z)   = dabs(z)

   1  format(a)
c-------------------------------- start executable statements ------------

      if( (x .le. zero) .or. (a .le. zero) .or. (b .lt. zero) ) then
         dgmdst = zero
         return
      endif

      xb = x/b
      gl = gamln(a)
      ta = exp(a*alog(xb) - xb - gl)
      if( xb .lt. (a + one) ) then
c						use series representation
         ap = a
         d  = one/a
         s  = d
         do 10 i = 1, ITMAX
            ap = ap + one
            d  = d*xb/ap
            s  = s + d
            if( d .lt. s*tol ) then
               dgmdst = s*ta
               return
            endif
  10     continue
      else
c						use continued fraction rep.
         ap = xb + one - a
         c  = one/small
         d  = one/ap
         h  = d
         do 20 i = 1, ITMAX
            an = -float(i)*(float(i) - a)
            ap = ap + two
            d  = an*d + ap
            if( abs(d) .lt. small ) d = small
            c  = ap + an/c
            if( abs(c) .lt. small ) c = small
            d  = one/d
            e  = d*c
            h  = h*e
            if( abs(e-one) .lt. tol ) then
               dgmdst = one - h*ta
               return
            endif
  20     continue
      endif

      write(6,1)'Gamma Distribution (dgmdst) unable to converge.'
      dgmdst = zero
      return
      end
