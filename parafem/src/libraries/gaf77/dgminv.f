c  ********************************************************************
c  *                                                                  *
c  *                         function dgminv                          *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 2.1
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: May 10, 1994
c
c  PURPOSE  returns the inverse Gamma distribution function for given
c           probability p and parameters alpha and beta.
c
c  This routine employs root finding techniques to obtain the inverse
c  Gamma distribution. The Gamma distribution is defined by
c
c                   1       x
c     F(x) =   ---------  INT exp{-t/b}*t**(a-1) dt
c              G(a)*b**a    0
c
c  so this routine returns x given p = F(x). Note that if p = 1 on input
c  the value 1.e+30 is returned. Relative error tolerance is fixed at 1.e-10.
c  Arguments to the routine are as follows;
c
c     p   real value giving the probability at which the corresponding x value
c         is desired. (input)
c
c     a   the alpha parameter of the Gamma distribution. (input)
c
c     b   the beta parameter of the Gamma distribution. (input)
c
c  ierr   integer error flag which is 0 if all goes well, -1 if convergence
c         not achieved.
c---------------------------------------------------------------------------
      real*8 function dgminv( p, a, b, ierr )
      parameter (MAXIT = 100)
      implicit real*8 (a-h,o-z)
      data zero/0.d0/,one/1.d0/,half/0.5d0/,tol/1.d-10/,small/1.d-20/
      data pt1/0.1d0/, opt5/1.5d0/, aone/0.999d0/, vbig/1.d+30/
      gmdst(z1,z2,z3) = dgmdst(z1,z2,z3)
      gamln(z) = dgamln(z)
      exp(z)   = dexp(z)
      abs(z)   = dabs(z)
      alog(z)  = dlog(z)

c-------------------------------- start executable statements ------------
c					eliminate the obvious
      if( p .eq. zero ) then
         dgminv = zero
         return
      elseif( p .eq. one ) then
         dgminv = vbig
         return
      endif
c					for Chi-Square, use bisection sooner
      if( a .lt. one ) then
         bg = one - pt1/a
      elseif( a .eq. one ) then
         dgminv = -b*alog(one - p)
         return
      else
         bg = aone
      endif
c					for small p, use bisection
      ierr = 0
      if( p .lt. small ) then
         xo = a*b
         do 20 j = 1, MAXIT
            x1 = half*xo
            if( gmdst(x1,a,b) .lt. p ) then
               do 10 i = 1, MAXIT
                  dgminv = half*(xo + x1)
                  if( abs((dgminv - xo)/dgminv) .lt. tol ) return
                  if( gmdst(dgminv,a,b) .lt. p ) then
                     x1 = dgminv
                  else
                     xo = dgminv
                  endif
  10           continue
               ierr = -1
               return
            endif
            xo = x1
  20     continue
c					for p almost 1.0, use bisection also
      elseif( p .gt. bg ) then
         xo = a*b
         do 40 j = 1, MAXIT
            x1 = opt5*xo
            if( gmdst(x1,a,b) .gt. p ) then
               do 30 i = 1, MAXIT
                  dgminv = half*(xo + x1)
                  if( abs((dgminv - xo)/dgminv) .lt. tol ) return
                  if( gmdst(dgminv,a,b) .gt. p ) then
                     x1 = dgminv
                  else
                     xo = dgminv
                  endif
  30           continue
               ierr = -1
               return
            endif
            xo = x1
  40     continue
      endif
c					otherwise use Newton-Raphson (modified)
      xo = a*b
      bg = gamln(a) + alog(b)
      bi = one/b
      c  = a - one
      do 60 i = 1, MAXIT
         xb = bi*xo
         fp = exp(c*alog(xb) - xb - bg)
         pg = p - gmdst(xo,a,b)
         dgminv = xo + pg/(fp + pg*(c/xo - bi))
         if( abs((dgminv - xo)/dgminv) .lt. tol ) return
         xo = dgminv
  60  continue
c					convergence not achieved
      ierr = -1
      return
      end
