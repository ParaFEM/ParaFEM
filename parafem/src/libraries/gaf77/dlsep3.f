c  *******************************************************************
c  *                                                                 *
c  *                     Function dlsep3                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.31
c  Written by Gordon A. Fenton, Mar. 24, 1994
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points in a 3-D random field
c           having a separable Markovian covariance function.
c
c  This function returns the covariance between two points in the field
c  separated by lag vector {X,Y,Z}. The covariance function for this field
c  is given by
c
c      B(X,Y,Z) = var*exp{-(2|X|/dthx)}*exp{-(2|Y|/dthy)}*exp{-(2|Z|/dthz)}
c
c  where var is the point variance, and dthx, dthy, and dthz are the scales of
c  fluctuation in the x-, y-, and z-directions, respectively. These parameters
c  are brought in through the common block /dparam/.
c
c  NOTE: this function is NOT isotropic even if dthx = dthy = dthz.
c
c  If var < 0, then this function computes the variance of a local average
c  of the random field, averaged over a domain having side dimensions {X,Y,Z}.
c  For a separable covariance function, the variance function is also
c  separable and can be written as the product of three 1-D variance
c  functions corresponding to the Markovian correlation function
c
c         V(X,Y,Z) = V(X)*V(Y)*V(Z)
c
c  where individual functions are based on the 1-D analytical model
c
c                     dthx^2
c         V(X) = 2 * ------- [(2|X|/dthx) + exp{-2|X|/dthx} - 1]
c                     4*X^2
c
c  The variance of the local average is then given by var*V(X,Y,Z). See page
c  186 of "Random Fields" by E. Vanmarcke and G.A. Fenton's thesis for more
c  details.
c
c  The data value `eund' is designed to avoid underflow errors arising from
c  exp(-a) for a > eund. We quite happily accept a zero in this case.
c
c  REVISION HISTORY
c  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (May 1/00)
c  1.11	minor changes to documentation above (Jul 14/00)
c  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Apr 9/01)
c  1.3	reversed default - now return covariances if var > 0. (Apr 11/01)
c  1.31	revised above documentation to reflect revision 1.3 (May 9/01)
c---------------------------------------------------------------------------
      real*8 function dlsep3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, half/0.5d0/, one/1.d0/, two/2.d0/, eight/8.d0/
      data eund/300.d0/
      exp(s) = dexp(s)
      abs(s) = dabs(s)

      if( var .lt. zero ) then			! return variance function
c							x-direction var func
         a = two*abs(X)/dthx
         if( X .eq. zero ) then
            dsepa = half
         elseif( a .lt. eund ) then
            dsepa = (a + exp(-a) - one)/(a*a)
         else
            dsepa = (a - one)/(a*a)
         endif
c							y-direction var func
         b = two*abs(Y)/dthy
         if( Y .eq. zero ) then
            dsepb = half
         elseif( b .lt. eund ) then
            dsepb = (b + exp(-b) - one)/(b*b)
         else
            dsepb = (b - one)/(b*b)
         endif

c							z-direction var func
         c = two*abs(Z)/dthz
         if( Z .eq. zero ) then
            dsepc = half
         elseif( c .lt. eund ) then
            dsepc = (c + exp(-c) - one)/(c*c)
         else
            dsepc = (c - one)/(c*c)
         endif
c							final variance function
         dlsep3 = -eight*var*dsepa*dsepb*dsepc

      else					! return covariance
         if( dthx .eq. zero ) then
            if( X .eq. zero ) then
               a = zero
            else
               dlsep2 = zero
               return
            endif
         else
            a = two*abs(X)/dthx
         endif
         if( dthy .eq. zero ) then
            if( Y .eq. zero ) then
               b = zero
            else
               dlsep2 = zero
               return
            endif
         else
            b = two*abs(Y)/dthy
         endif
         if( dthz .eq. zero ) then
            if( Z .eq. zero ) then
               c = zero
            else
               dlsep2 = zero
               return
            endif
         else
            b = two*abs(Z)/dthz
         endif
         abc = a + b + c
         if( abc .gt. eund ) then
            dlsep3 = zero
         else
            dlsep3 = var*exp(-abc)
         endif
      endif

      return
      end
