c  ********************************************************************
c  *                                                                  *
c  *                         function studi                           *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Nov. 4, 1993
c  Latest Update: Dec 5, 1996
c
c  PURPOSE  returns the inverse Student-t distribution function for given
c           probability p and parameter k
c
c  This routine employs the inverse of the incomplete beta function `dbetnv'
c  to compute the inverse distribution function, ie the value of x
c  corresponding to a given probability p, where X follows a Student-t
c  distribution with k degrees of freedom. The Student-t distribution is
c
c                  G((k+1)/2)       x              1
c     F(x) =   -----------------  INT ------------------------- dt
c              G(k/2)*sqrt{pi*k} -inf  [ 1 + t**2/k]**((k+1)/2)
c
c  where G(a) is the Gamma function ((a-1)! if `a' is an integer).
c  This algorithm uses a Newton-Raphson iteration with first guess as
c  provided by Abramowitz and Stegun, pg. 949.
c
c  Note: if p = 0 or 1, then this function returns +/- 1.e+20.
c
c  Arguments to the routine are as follows;
c
c     p   real value giving the probability at which the inverse distribution
c         is desired. (input)
c
c     k   the number of degrees of freedom. (input)
c
c  REVISION HISTORY:
c  1.1	eliminated unused variable `two' (Dec 5/96)
c---------------------------------------------------------------------------
      real function studi(p,k)
      parameter (MAXIT = 50)
      data zero/0.0/, half/0.5/, one/1.0/, vbig/1.e+20/
      data pi/3.1415926535897932384/, tol/1.e-6/
      g1(x) = 0.25*(x**3 + x)
      g2(x) = (5.0*x**5 + 16.0*x**3 + 3.0*x)/96.0
      g3(x) = (3.0*x**7 + 19.0*x**5 + 17.0*x**3 - 15.0*x)/384.0

c					eliminate the obvious
      if( p .le. zero ) then
         studi = -vbig
         return
      elseif( p .ge. one ) then
         studi = vbig
         return
      elseif( p .eq. half ) then
         studi = zero
         return
      endif
c					use Newton-Raphson (modified)
c						first get initial guess
      v  = one/float(k)
      xp = phinv(p)
      xp = xp + v*(g1(xp) + v*(g2(xp) + v*g3(xp)))

      a  = half*float(k+1)
      b  = half*float(k)
      bg = gamma(a)/(gamma(b)*sqrt(pi*float(k)))
      xo = xp
      do 60 i = 1, MAXIT
         gx = p - dstudt(xo,k)
         gp = ((one + xo*xo*v)**a)/bg
         studi = xo + gx*gp
         if( abs((studi - xo)/studi) .lt. tol ) return
         xo = studi
  60  continue
c					convergence not achieved
      studi = xp
      return
      end
