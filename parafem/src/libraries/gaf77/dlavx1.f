c  *******************************************************************
c  *                                                                 *
c  *                     Real*8 Function dlavx1                      *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.62
c  Written by Gordon A. Fenton, 1989
c  Latest Update: Jul 10, 2002
c
c  PURPOSE  returns the covariance between two points along a 1-D Markov
c           process.
c
c  Returns the covariance between two points, separated by distance T, in
c  a 1-D Markov process. The covariance function is given by
c
c       B(T) = var * exp( -2|T|/dthx )
c
c  where var is the point variance and dthx is the scale of fluctuation.
c  var and dthx are brought into this routine through the common
c  block DPARAM.
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(T), averaged over the distance T, where V(T) is
c  the variance function (see E.H. Vanmarcke, "Random Fields: Analysis and
c  Synthesis", MIT Press, 1984, Chapter 5).
c                                  T
c       V(T) = (2/T) * INT  (1-s/T)*r(s) ds
c                                  0
c  where r(s) is the correlation function: r(s) = B(s)/var. The integration
c  gives us
c
c       V(T) = (2/a*a) * [ a + exp(-a) - 1]
c
c  where a = 2*T/dthx. When dthx >> T, a third order Taylor's series expansion
c  gives the numerically improved approximation
c
c	V(T) = [1 - 2*|T|/(3*dthx) ],	dthx >> T
c
c  and when dthx << T, the exp(-a) term disappears because `a' becomes very
c  large, leaving us with
c
c	V(T) = (dthx/T) * [ 1 - (dthx/2T) ]
c
c
c  REVISION HISTORY:
c  1.1	explicitly accomodate the case when dthx = 0.
c  1.2	corrected bug introduced in 1.1 (missed multiplying by var)
c  1.3	improved limit handling
c  1.4	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (Apr 29/00)
c  1.5	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
c  1.6	reversed default - now return covariances if var > 0 (Apr 11/01)
c  1.61	revised above documentation to reflect revision 1.6 (May 9/01)
c  1.62	corrected erroneous documentation above re variance func (Jul 10/02)
c =======================================================================
      real*8 function dlavx1( T )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, half/0.5d0/, one/1.d0/, two/2.d0/, three/3.d0/
      data small/0.0067d0/, big/1.d6/
      abs(x) = dabs(x)
      exp(x) = dexp(x)

      aT = abs(T)
      if( var .lt. zero ) then			! return variance function
         if( dthx .ge. big*aT ) then
            if( dthx .eq. zero ) then
               dlavx1 = -var
            else
               dlavx1 = -var*(one - two*aT/(three*dthx))
            endif
         elseif( dthx .le. small*aT ) then
            rx = dthx/aT
            dlavx1 = -var*rx*(one - half*rx)
         else
            a = dthx/abs(T)
            b = half*a
            c = exp(-one/b) - one
            dlavx1 = -var*a*(one + b*c)
         endif
      else					! return covariance
         if( dthx .eq. zero ) then
            if( T .eq. zero ) then
               dlavx1 = var
            else
               dlavx1 = zero
            endif
         else
            dlavx1 = var*exp( -two*aT/dthx )
         endif
      endif

      return
      end
