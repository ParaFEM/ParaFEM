c  *******************************************************************
c  *                                                                 *
c  *                     Function dlspx1                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.41
c  Written by Gordon A. Fenton, Mar. 25, 1994
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points along a 1-D random
c           process having an exponential (Gaussian) type covariance function
c
c  This function returns the covariance between two points, separated by
c  distance T, along a random process having Gaussian type covariance function
c
c      B(T) = var * exp{ -pi*[(|T|/dthx)^2] }
c
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(T), averaged over the distance T, where the
c  variance function is given by
c
c         V(T) = (1/a*a) * [pi*(|T|/dthx)*erf(a) + exp(-a*a) - 1]
c
c  where a = (|T|*sqrt(pi))/dthx, erf(.) is the error function,
c  dthx is the directional scale of fluctuation, and var is the point
c  variance. The variance of a local average is given by var*V(X).
c
c  See page 186 of "Random Fields" by E. Vanmarcke and G.A. Fenton's thesis
c  for more details.
c
c  The parameters var and dthx are brought in via the common block
c  /dparam/.
c
c  The argument to this routine, T, is the distance between the two points
c  for which the covariance is desired (or the distance over which the local
c  averaging is performed, if var < 0).
c
c  REVISION HISTORY:
c  1.1	fixed argument list and dthx = 0 case, improved limit handling
c	(Jun 19/96)
c  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (May 1/00)
c  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
c  1.4	reversed default - now return covariances if var > 0. (Apr 11/01)
c  1.41	revised above documentation to reflect revision 1.4 (May 9/01)
c---------------------------------------------------------------------------
      real*8 function dlspx1( T )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, one/1.d0/, six/6.d0/
      data rtpi/1.77245385090551588d0/, pt1/0.1d0/
      data   pi/3.1415926535897932384d0/, big/1.d6/
      exp(s)   = dexp(s)
      abs(s)   = dabs(s)
      r_erf(s) = derf(s)

      if( var .lt. zero ) then			! return variance function
         aT = abs(T)
         if( dthx .ge. big*aT ) then
            if( dthx .eq. zero ) then
               dlspx1 = -var
            else
               rx = aT/dthx
               dlspx1 = -var*(one - pi*rx*rx/six)
            endif
         elseif( dthx .le. pt1*aT ) then
            rx = dthx/aT
            dlspx1 = -var*rx*(one - rx/pi)
         else
            rx = aT/dthx
            a  = rtpi*rx
            aa = a*a
            dlspx1 = -var*(pi*rx*r_erf(a) + exp(-aa) - one)/aa
         endif
      else					! return covariance
         if( dthx .eq. zero ) then
            if( T .eq. zero ) then
               dlspx1 = var
            else
               dlspx1 = zero
            endif
         else
            at = T/dthx
            dlspx1 = var * exp(-pi*at*at)
         endif
      endif

      return
      end

