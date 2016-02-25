c  ********************************************************************
c  *                                                                  *
c  *                         function studt                           *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.0
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
      real function studt(x,k)
      data zero/0.0/, half/0.5/, one/1.0/

      if( x .eq. zero ) then
         studt = half
         return
      endif

      v = float(k)
      y = v/(v + x*x)

      studt = half*beta( y, half*v, half )

      if( x .gt. zero ) studt = one - studt

      return
      end
