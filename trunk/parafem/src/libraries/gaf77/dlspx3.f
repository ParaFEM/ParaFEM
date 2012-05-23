c  *******************************************************************
c  *                                                                 *
c  *                     Function dlspx3                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.42
c  Written by Gordon A. Fenton, Mar. 25, 1994
c  Latest Update: Sep 21, 2001
c
c  PURPOSE  returns the covariance between two points in a 3-D random field
c           having an exponential (Gaussian) type covariance function
c           (separable and mean-square differentiable).
c
c  This function returns the covariance between two points separated
c  by lag vector {X,Y,Z}. The process has a separable and mean-square
c  differentiable covariance function,
c
c      B(X,Y,Z) = var * exp{ -pi*[(|X|/thx)^2 + (|Y|/thy)^2 + (|Z|/thz)^2] }
c
c  where var is the point variance, thx, thy, and thz are the directional
c  scales of fluctuation and (X,Y,Z) are elements of the lag vector
c  separating the two points of interest.
c
c  The parameters var, dthx, dthy, and dthz are brought into this
c  routine through the common block /dparam/.
c
c  If var < 0, then this function computes the variance of a local average
c  of the random field, averaged over a domain having side dimensions {X,Y,Z}.
c  The variance function is also separable and can be written
c  as the product of three one-d variance functions
c
c         V(X,Y,Z) = V(X)*V(Y)*V(Z)
c
c  each corresponding to the 1-D Gaussian type correlation function,
c
c         p(x) = exp{-pi*(|x|/thx)^2}
c
c  The individual variance functions are based on the 1-D analytical model
c
c         V(X) = (1/a*a) * [pi*(|X|/thx)*erf(a) + exp(-a*a) - 1]
c
c  where a = (|X|*sqrt(pi))/thx and erf(.) is the error function.
c  See page 186 of "Random Fields" by E. Vanmarcke and G.A. Fenton's thesis
c  for more details. If var < 0, this function returns |var|*V(X,Y,Z).
c
c  The arguments to this routine, X, Y, and Z, are the components of the
c  lag vector between the two points of interest (or the dimensions of the
c  averaging volume, if var < 0)
c
c  The data value `eund' is designed to avoid underflow errors arising from
c  exp(-a) for a > eund. We quite happily accept a zero in this case. This
c  tends to be compiler dependent.
c
c  Requires:
c    1) from Fortran lib: d_erf	(error function), define the following
c	r_erf(s) = d_erf(s) (on SunOS 4.1.3)
c	r_erf(s) = derf(s)  (on DEC OSF/1)
c
c  REVISION HISTORY:
c  1.1	eliminated unused local variables `two' and `four' (Dec 5/96)
c  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (May 31/00)
c  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Apr 9/01)
c  1.4	reversed default - now return covariances if var > 0. (Apr 11/01)
c  1.41	revised above documentation to reflect revision 1.4 (May 9/01)
c  1.42	replaced 'dthx' with 'thx' (etc) below (Sep 21/01)
c---------------------------------------------------------------------------
      real*8 function dlspx3( X, Y, Z )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, thx, thy, thz
      data zero/0.d0/, one/1.d0/
      data rtpi/1.77245385090551588d0/, eund/300.d0/
      data pi/3.1415926535897932384d0/
      exp(s)   = dexp(s)
      abs(s)   = dabs(s)
      r_erf(s) = derf(s)

      if( var .lt. zero ) then			! return variance of local ave
         rx = abs(X)/thx
         a  = rtpi*rx
         aa = a*a
         ry = abs(Y)/thy
         b  = rtpi*ry
         bb = b*b
         rz = abs(Z)/thz
         c  = rtpi*rz
         cc = c*c
c							in the x direction
         if( X .eq. zero ) then
            dsepa = one
         elseif( aa .lt. eund ) then
            dsepa = (pi*rx*r_erf(a) + exp(-aa) - one)/aa
         else
            dsepa = (pi*rx - one)/aa
         endif
c							in the y direction
         if( Y .eq. zero ) then
            dsepb = one
         elseif( bb .lt. eund ) then
            dsepb = (pi*ry*r_erf(b) + exp(-bb) - one)/bb
         else
            dsepb = (pi*ry - one)/bb
         endif
c							and in the z direction
         if( Z .eq. zero ) then
            dsepc = one
         elseif( cc .lt. eund ) then
            dsepc = (pi*rz*r_erf(c) + exp(-cc) - one)/cc
         else
            dsepc = (pi*rz - one)/cc
         endif
         dlspx3 = -var*dsepa*dsepb*dsepc

      else					! return covariance
c							in the x direction
         if( thx .eq. zero ) then
            if( X .eq. zero ) then
               dsepa = one
            else
               dsepa = zero
            endif
         else
            rx = X/thx
            dsepa = exp(-pi*rx*rx)
         endif
c							in the y direction
         if( thy .eq. zero ) then
            if( Y .eq. zero ) then
               dsepb = one
            else
               dsepb = zero
            endif
         else
            ry = Y/thy
            dsepb = exp(-pi*ry*ry)
         endif
c							in the y direction
         if( thz .eq. zero ) then
            if( Z .eq. zero ) then
               dsepc = one
            else
               dsepc = zero
            endif
         else
            rz = Z/thz
            dsepc = exp(-pi*rz*rz)
         endif
         dlspx3 = var*dsepa*dsepb*dsepc
      endif


      return
      end


