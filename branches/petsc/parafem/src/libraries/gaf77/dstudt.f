c  ********************************************************************
c  *                                                                  *
c  *                         function dstudt                          *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Nov. 4, 1993
c
c  PURPOSE  returns the Student-t distribution function for given x and
c           parameter k
c
c  This routine employs the incomplete beta function `dbeta' to compute the
c  cumulative probability associated with a variable following the Student-t
c  distribution,
c
c                  G((k+1)/2)       x              1
c     F(x) =   -----------------  INT ------------------------- dt
c              G(k/2)*sqrt{pi*k} -inf  [ 1 + t**2/k]**((k+1)/2)
c
c  where G(a) is the Gamma function ((a-1)! if `a' is an integer).
c  Arguments to the routine are as follows;
c
c     x   real value giving the position at which the cumulative distribution
c         is desired. (input)
c
c     k   the number of degrees of freedom. (input)
c
c---------------------------------------------------------------------------
      real*8 function dstudt(x,k)
      implicit real*8 (a-h,o-z)
      data zero/0.d0/, half/0.5d0/, one/1.d0/
      beta(z1,z2,z3) = dbeta(z1,z2,z3)

      if( x .eq. zero ) then
         dstudt = half
         return
      endif

      v = dble(k)
      y = v/(v + x*x)

      dstudt = half*beta( y, half*v, half )

      if( x .gt. zero ) dstudt = one - dstudt

      return
      end
