c  *******************************************************************
c  *                                                                 *
c  *                     Function dlavx3                             *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.53
c  Written by Gordon A. Fenton, Sept. 1989
c  Latest Update: Jul 2, 2003
c
c  PURPOSE  returns the covariance between two points in a 3-D volume for
c           a random field having a Markovian covariance function.
c
c  This function returns the covariance between two
c  points in the field separated by lag vector {X,Y,Z}.
c  The 3-D Markovian process has covariance function
c
c   B(X,Y,Z) = var* exp{ -sqrt[(2|X|/dthx)^2 + (2|Y|/dthy)^2 + (2|Z|/dthz)^2] }
c
c  where var is the point variance of the process and dthx, dthy, and
c  dthz are the directional scales of fluctuation. The parameters var,
c  dthx, dthy, and dthz are brought in via the common block /dparam/.
c
c  Note that if var < 0, as provided in the common block /dparam/,
c  then this function returns the variance of a local average of the
c  process over a volume of side dimensions {X,Y,Z}. The `local average'
c  variance is obtained by Gaussian quadrature integration of the above
c  covariance function (the 6-D integral is collapsed to a 3-D integral by
c  taking advantage of the fact that the above function is quadrant
c  symmetric.
c  See "Random Fields" by E. Vanmarcke and G.A. Fenton's Ph.D. thesis for
c  more details on local average theory.
c
c  The arguments to this routine are just the elements of the lag vector
c  between points for which the covariance is desired.
c
c  REVISION HISTORY:
c  1.1	brought the pi/2 factor into the function directly, rather than
c	making it a parameter (was `dum1').
c  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
c	true if the variance function is to be returned. Otherwise this
c	function returns the covariance. (Apr 29/00)
c  1.3	covariance version now returns isotropic Markovian model (Jun 22/00)
c  1.4	eliminated lvarfn -- now return covariances only if var < 0 (Apr 9/01)
c  1.5	reversed default - now return covariances if var > 0.
c	Variance function values now computed by Gauss quadrature (Apr 11/01)
c  1.51	revised the above description to reflect Revision 1.5 (May 9/01)
c  1.52	now use 16 pt Gauss quadrature for variance calcs (Apr 8/03)
c  1.53	now use 20 pt Gauss quadrature for variance calcs (Jul 2/03)
c---------------------------------------------------------------------------
      real*8 function dlavx3( X, Y, Z )
      parameter (NG = 20)
      implicit real*8 (a-h,o-z)
      dimension w(NG), q(NG)
      common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, one/1.d0/, two/2.d0/
      data eighth/0.125d0/, half/0.5d0/
c			Gauss weights and points for NG = 5
c     data w/ .236926885056189d0,.478628670499366d0,.568888888888889d0,
c    >        .478628670499366d0,.236926885056189d0/
c     data q/-.906179845938664d0,-.538469310105683d0,.000000000000000d0,
c    >        .538469310105683d0, .906179845938664d0/
c			these are for NG = 16
c     data w/0.027152459411754094852d0, 0.062253523938647892863d0,
c    >       0.095158511682492784810d0, 0.124628971255533872052d0,
c    >       0.149595988816576732081d0, 0.169156519395002538189d0,
c    >       0.182603415044923588867d0, 0.189450610455068496285d0,
c    >       0.189450610455068496285d0, 0.182603415044923588867d0,
c    >       0.169156519395002538189d0, 0.149595988816576732081d0,
c    >       0.124628971255533872052d0, 0.095158511682492784810d0,
c    >       0.062253523938647892863d0, 0.027152459411754094852d0/
c     data q/-.989400934991649932596d0, -.944575023073232576078d0,
c    >       -.865631202387831743880d0, -.755404408355003033895d0,
c    >       -.617876244402643748447d0, -.458016777657227386342d0,
c    >       -.281603550779258913230d0, -.095012509837637440185d0,
c    >       0.095012509837637440185d0, 0.281603550779258913230d0,
c    >       0.458016777657227386342d0, 0.617876244402643748447d0,
c    >       0.755404408355003033895d0, 0.865631202387831743880d0,
c    >       0.944575023073232576078d0, 0.989400934991649932596d0/
c			these are for NG = 20
      data w/0.017614007139152118312d0, 0.040601429800386941331d0,
     >       0.062672048334109063570d0, 0.083276741576704748725d0,
     >       0.101930119817240435037d0, 0.118194531961518417312d0,
     >       0.131688638449176626898d0, 0.142096109318382051329d0,
     >       0.149172986472603746788d0, 0.152753387130725850698d0,
     >       0.152753387130725850698d0, 0.149172986472603746788d0,
     >       0.142096109318382051329d0, 0.131688638449176626898d0,
     >       0.118194531961518417312d0, 0.101930119817240435037d0,
     >       0.083276741576704748725d0, 0.062672048334109063570d0,
     >       0.040601429800386941331d0, 0.017614007139152118312d0/
      data q/-.993128599185094924786d0, -.963971927277913791268d0,
     >       -.912234428251325905868d0, -.839116971822218823395d0,
     >       -.746331906460150792614d0, -.636053680726515025453d0,
     >       -.510867001950827098004d0, -.373706088715419560673d0,
     >       -.227785851141645078080d0, -.076526521133497333755d0,
     >       0.076526521133497333755d0, 0.227785851141645078080d0,
     >       0.373706088715419560673d0, 0.510867001950827098004d0,
     >       0.636053680726515025453d0, 0.746331906460150792614d0,
     >       0.839116971822218823395d0, 0.912234428251325905868d0,
     >       0.963971927277913791268d0, 0.993128599185094924786d0/

      abs(qq)  = dabs(qq)
      exp(qq)  = dexp(qq)
      sqrt(qq) = dsqrt(qq)

      thx = dthx
      thy = dthy
      thz = dthz
      aX  = abs(X)
      aY  = abs(Y)
      aZ  = abs(Z)

      if( dthx .eq. zero ) then
         if( X .eq. zero ) then
            thx = one
         else
            dlavx3 = zero
            return
         endif
      endif
      if( dthy .eq. zero ) then
         if( Y .eq. zero ) then
            thy = one
         else
            dlavx3 = zero
            return
         endif
      endif
      if( dthz .eq. zero ) then
         if( Z .eq. zero ) then
            thz = one
         else
            dlavx3 = zero
            return
         endif
      endif

      if( var .lt. zero ) then			! return variance function
         rx = half*aX
         ry = half*aY
         rz = half*aZ
         tx = two/thx
         ty = two/thy
         tz = two/thz
         d1 = zero
         do 30 i = 1, NG
            xi = rx*(one + q(i))
            a1 = tx*xi
            d2 = zero
            do 20 j = 1, NG
               yj = ry*(one + q(j))
               a2 = ty*yj
               d3 = zero
               do 10 k = 1, NG
                  zk = rz*(one + q(k))
                  a3 = tz*zk
                  zp = sqrt(a1*a1 + a2*a2 + a3*a3)
                  d3 = d3 + w(k)*(one-q(k))*exp(-zp)
  10           continue
               d2 = d2 + w(j)*(one-q(j))*d3
  20        continue
            d1 = d1 + w(i)*(one-q(i))*d2
  30     continue
         dlavx3 = -eighth*var*d1
      else					! return covariance
         a1 = two*aX/thx
         a2 = two*aY/thy
         a3 = two*aZ/thz
         zp = sqrt(a1*a1 + a2*a2 + a3*a3)
         dlavx3 = var*exp( -zp )
      endif

      return
      end
