c  *********************************************************************
c  *                                                                   *
c  *                      real*8 function dlfsc1                       *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.51
c  Written by Gordon A. Fenton, TUNS, Nov 19, 1997
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns covariance between two points in a Fractal-Scale random
c           field.
c
c  DESCRIPTION
c  This function returns the covariance between
c  two points, separated by distance T, in a Fractal-Scale random field.
c  The Fractal-Scale process has covariance function
c
c                           var
c                B(T) = ----------------
c                       ( 1 + |T/pb| )^F
c
c  where F is the fractal-scale parameter and pb is a non-dimensionalizing
c  length scale. When 0 <= F <= 1 we get fractal type behaviour (correlation
c  decays as 1/T^F), in which case pb can be interpreted as a length below
c  which self-similar behaviour no longer holds, while for F > 1 we get
c  finite-scale behaviour with scale of fluctuation 2*pb/(F-1).
c  The parameters var, pb, and F are brought into this routine through
c  the common block /dparam/.
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(T), averaged over the distance T.
c  The variance function corresponding to the above covariance function
c  is
c
c                        2*pb^F      1        2-F     2-F        2-F
c                V(T) = --------- [ --- {(pb+T)   - pb   } - T*pb    ]
c                       T*T*(1-F)   2-F
c
c  with special cases for T = 0, F = 1, and F = 2.
c
c  ARGUMENTS
c
c	T	real value giving the distance between points (or the
c		size of the averaging domain, if var < 0). (input)
c
c  REVISION HISTORY:
c  1.1	correctly non-dimensionalized rho using additional parameter pb
c	(Jan 8/98)
c  1.2	corrected to divide by tmF (not tF!) on last assign line (Sep 19/98)
c  1.3	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (May 1/00)
c  1.4	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
c  1.5	reversed default - now return covariances if var > 0. (Apr 11/01)
c  1.51	revised above documentation to reflect revision 1.5 (May 9/01)
c-------------------------------------------------------------------------
      real*8 function dlfsc1(T)
      implicit real*8 (a-h,o-z)
      common/dparam/ var, pb, F, da, db
      data zero/0.d0/, one/1.d0/, two/2.d0/
      abs(x)  = dabs(x)
      alog(x) = dlog(x)

      aT = abs(T)
      if( var .lt. zero ) then			! return variance function
         if( T .eq. zero ) then
            dlfsc1 = -var
         elseif( F .eq. one ) then
            ama    = alog(pb+aT) - alog(pb)
            dlfsc1 = -var*two*pb*( (pb+aT)*ama - aT )/(aT*aT)
         elseif( F .eq. two ) then
            ama    = alog(pb+aT) - alog(pb)
            dlfsc1 = -var*two*pb*( aT - pb*ama )/(aT*aT)
         else
            omF    = one - F
            tmF    = two - F
            ama    = (pb + aT)**tmF - pb**tmF
            pbf    = pb**F
            pbf1   = pb**omF
            dlfsc1 = -var*two*pbf*( (ama/tmF) - aT*pbf1 )/(aT*aT*omF)
         endif
      else					! return covariance
         dlfsc1 = var/(one + (aT/pb)**F)
      endif

      return
      end
