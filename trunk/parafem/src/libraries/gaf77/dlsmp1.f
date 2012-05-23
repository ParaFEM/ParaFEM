c  *******************************************************************
c  *                                                                 *
c  *                     Real*8 Function dlsmp1                      *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.31
c  Written by Gordon A. Fenton, May 13, 1996
c  Latest Update: May 9, 2001
c
c  PURPOSE  returns the covariance between two points along a 1-D simple
c           variance function process 
c
c  This routine returns the covariance between two points, separated by
c  the distance T, along a process having the simple covariance function
c
c                  var*dthx^3
c	  B(T) = --------------
c                (dthx + |T|)^3
c
c  where var is the point variance, and dthx is the scale of fluctuation
c  of the process. Both parameters, var and dthx, are brought into this
c  function via the common block /dparam/.
c
c  If var < 0, then this function returns the variance of a local average
c  of the process, |var|*V(T), averaged over the distance T.
c  The random process considered has a particularly simple variance
c  function,
c
c                     dthx
c	    V(T) = ----------
c                  dthx + |T|
c
c  where `dthx' is the scale of fluctuation of the process. This function
c  returns var*V(T), where `var' is the point variance of the process.
c  For more details, see page 204 of E.H. Vanmarcke's ``Random Fields:
c  Analysis and Synthesis.''
c
c  REVISION HISTORY:
c  1.1	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (May 1/00)
c  1.2	eliminated lvarfn -- now return covariances only if var < 0 (Mar 2/01)
c  1.3	reversed default - now return covariances if var > 0. (Apr 11/01)
c  1.31	revised above documentation to reflect revision 1.3 (May 9/01)
c =======================================================================
      real*8 function dlsmp1( T )
      implicit real*8 (a-h,o-z)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/
      abs(x) = dabs(x)

      if( var .lt. zero ) then			! return variance function
         dlsmp1 = -var*dthx/(dthx + abs(T))
      else					! return covariance
         at = dthx + abs(T)
         if( at .eq. zero ) then
            dlsmp1 = var
         else
            bt = dthx/at
            dlsmp1 = var*bt*bt*bt
         endif
      endif

      return
      end
