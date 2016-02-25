c  *******************************************************************
c  *                                                                 *
c  *                     Function dlspx2                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.71
c  Written by Gordon A. Fenton, Sept. 1989
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points in a 2-D random field
c           with Gaussian type covariance function (separable and mean-square
c           differentiable).
c
c  This function returns the covariance between two points, separated
c  by lag vector {X,Y}, in a 2-D random field having separable (and mean-
c  square differentiable) covariance function
c
c      B(x,y) = var * exp{ -pi*(|x|/thx)^2 } * exp{ -pi*(|y|/thy)^2 }
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(X,Y), averaged over the domain X x Y.
c  Thus the variance function is also separable and can be written
c  as the product of two one-d variance functions corresponding to the
c  exponential correlation function
c
c         V(X,Y) = V(X)*V(Y)
c
c  where individual functions are based on the 1-D analytical model
c
c         V(X) = (1/a*a) * [pi*(|X|/thx)*erf(a) + exp(-a*a) - 1]
c
c  where a = (|X|*sqrt(pi))/thx and erf(.) is the error function,
c  dthx is the directional scale of fluctuation, and var is the
c  point variance. The variance of a local average is given by
c  var*V(X,Y).
c
c  See page 186 of "Random Fields" by E. Vanmarcke and G.A. Fenton's
c  thesis for more details.
c
c  The parameters var, dthx, and dthy are brought into this
c  routine through the common block /dparam/.
c
c  The arguments to this routine, X and Y, are the components of the lag
c  vector between the two points (or the dimensions of the averaging
c  region, if var < 0).
c
c  The data value `eund' is designed to avoid underflow errors arising from
c  exp(-a) for a > eund. We quite happily accept a zero in this case. This
c  tends to be compiler dependent.
c
c  REVISION HISTORY:
c  1.3	eliminated unused local variables `two' and `four' (Dec 5/96)
c  1.4	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (May 31/00)
c  1.5	eliminated lvarfn -- now return covariances only if var < 0 (Mar 27/01)
c  1.6	properly handled sign on var for covariances (Apr 5/01)
c  1.7	reversed default - now return covariances if var > 0. (Apr 11/01)
c  1.71	revised above documentation to reflect revision 1.7 (May 9/01)
c---------------------------------------------------------------------------
      real*8 function dlspx2( X, Y )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, dthx, dthy, dthz
      data zero/0.d0/, one/1.d0/
      data rtpi/1.77245385090551588d0/, eund/300.d0/
      data pi/3.1415926535897932384d0/
      exp(s)   = dexp(s)
      abs(s)   = dabs(s)
      r_erf(s) = derf(s)

      if( var .lt. zero ) then			! return variance function
         rx = abs(X)/dthx
         ry = abs(Y)/dthy
         a  = rtpi*rx
         aa = a*a
         b  = rtpi*ry
         bb = b*b
c						in the x direction
         if( X .eq. zero ) then
            dsepa = one
         elseif( aa .lt. eund ) then
            dsepa = (pi*rx*r_erf(a) + exp(-aa) - one)/aa
         else
            dsepa = (pi*rx - one)/aa
         endif
c						and in the y direction
         if( Y .eq. zero ) then
            dsepb = one
         elseif( bb .lt. eund ) then
            dsepb  = (pi*ry*r_erf(b) + exp(-bb) - one)/bb
         else
            dsepb = (pi*ry - one)/bb
         endif
         dlspx2 = -var*dsepa*dsepb

      else					! var < 0, return covariance
c							in the x direction
         if( dthx .eq. zero ) then
            if( X .eq. zero ) then
               dsepa = one
            else
               dsepa = zero
            endif
         else
            rx = X/dthx
            px = pi*rx*rx
            if( px .lt. eund ) then
               dsepa = exp(-px)
            else
               dsepa = zero
            endif
         endif
c							in the y direction
         if( dthy .eq. zero ) then
            if( Y .eq. zero ) then
               dsepb = one
            else
               dsepb = zero
            endif
         else
            ry = Y/dthy
            py = pi*ry*ry
            if( py .lt. eund ) then
               dsepb = exp(-py)
            else
               dsepb = zero
            endif
         endif
         dlspx2 = var*dsepa*dsepb
      endif

      return
      end
