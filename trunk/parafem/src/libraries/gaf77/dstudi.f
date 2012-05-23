c  ********************************************************************
c  *                                                                  *
c  *                         function dstudi                          *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Nov. 4, 1993
c  Latest Update: May 16, 1997
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
c  Note: if p = 0 or 1, then this function returns +/- 1.e+30.
c
c  Arguments to the routine are as follows;
c
c     p   real value giving the probability at which the inverse distribution
c         is desired. (input)
c
c     k   the number of degrees of freedom. (input)
c
c  REVISION HISTORY:
c  1.1	eliminated unused local variable `two' (Dec 5/96)
c  1.2	employ ln(Gamma) function dgamln to handle big arguments. (May 16/97)
c---------------------------------------------------------------------------
      real*8 function dstudi(p,k)
      parameter (MAXIT = 50)
      implicit real*8 (a-h,o-z)
      common/dbgstd/ i
      data zero/0.d0/, half/0.5d0/, one/1.d0/, vbig/1.d+30/
      data pi/3.1415926535897932384d0/, tol/1.d-10/
      g1(x) = 0.25d0*(x**3 + x)
      g2(x) = (5.d0*x**5 + 16.d0*x**3 + 3.d0*x)/96.d0
      g3(x) = (3.d0*x**7 + 19.d0*x**5 + 17.d0*x**3 - 15.d0*x)/384.d0
      float(i) = dble(i)
      sqrt(x)  = dsqrt(x)
      phinv(x) = dphinv(x)
      gamln(x) = dgamln(x)
      exp(x)   = dexp(x)

   1  format(a,f9.7,a,i10)
c					eliminate the obvious
      if( p .le. zero ) then
         dstudi = -vbig
         return
      elseif( p .ge. one ) then
         dstudi = vbig
         return
      elseif( p .eq. half ) then
         dstudi = zero
         return
      endif
c					use Newton-Raphson (modified)
c						first get initial guess
      xp = phinv(p)
c						just return for very large k
      if( k .gt. 1000 ) then
         dstudi = xp
         return
      endif
c						otherwise iterate
      v  = one/float(k)
      xp = xp + v*(g1(xp) + v*(g2(xp) + v*g3(xp)))

      a  = half*float(k+1)
      b  = half*float(k)
      bg = exp( gamln(a) - gamln(b) )/sqrt(pi*float(k))
      xo = xp
      do 60 i = 1, MAXIT
         gx = p - dstudt(xo,k)
         gp = ((one + xo*xo*v)**a)/bg
         dstudi = xo + gx*gp
         if( abs((dstudi - xo)/dstudi) .lt. tol ) return
         xo = dstudi
  60  continue
c					convergence not achieved

      write(6,1)'dstudi: convergence not achieved for p = ',p,', k = ',k
      write(6,1)'         Returning inverse normal distribution.'
      dstudi = xp
      return
      end


