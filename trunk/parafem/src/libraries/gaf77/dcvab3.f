c  *********************************************************************
c  *                                                                   *
c  *                           function dcvab3                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.5
c  Written by Gordon A. Fenton, TUNS, July 17, 1992
c  Latest Update: Jul 2, 2003
c
c  PURPOSE  returns the covariance between two 3-D local averages, one
c           having eight times the volume of the other. Used by LAS3G.
c
c  This function evaluates the covariance between two local averages in
c  3-dimensional space. One local average is assumed to be of size Dx x Dy x Dz
c  and the other of size 2*Dx x 2*Dy x 2*Dz and they are separated in space by
c  a central distance having components Tx = C1*Dx, Ty = C2*Dy, and
c  Tz = C3*Dz.
c
c  The covariance is obtained by a 16-pt Gaussian quadrature of a
c  6-D integral collapsed (by assuming the covariance function to be
c  quadrant symmetric) to a 3-D integral.
c
c  NOTE: if the covariance function is not quadrant symmetric, this
c        routine will return erroneous results. Quadrant symmetry means
c        that cov(X,Y,Z) = cov(-X,Y,Z) = cov(X,-Y,Z), etc, where
c        cov(X,Y,Z) is a function returning the covariance between
c        points separated by lag (X,Y,Z), as discussed next.
c
c  The variance/covariance function is referenced herein using a call of the
c  form
c
c          V = cov( X, Y, Z )
c
c  where X, Y, and Z are the lag distances between the points in the
c  field.
c  Parameters, such as var, pb, dthx, dthy, and dthz are passed to the
c  variance/covariance function and this routine via the common block
c  'dparam'.
c
c  Arguments to this function are as follows
c
c    cov       external function provided by the user. On each invocation,
c              this routine calls cov up to 27*(5^3) = 3375 times.
c
c    Dx        x-dimension of the smaller local average cell. (input)
c
c    Dy        y-dimension of the smaller local average cell. (input)
c
c    Dz        z-dimension of the smaller local average cell. (input)
c
c    C1        x-direction distance between local average centers is C1*Dx.
c              (input)
c
c    C2        y-direction distance between local average centers is C2*Dy.
c              (input)
c
c    C3        z-direction distance between local average centers is C3*Dz.
c              (input)
c
c  REVISION HISTORY:
c  1.1	now including a Gaussian Quadrature integration option as an
c	alternative to the variance function approach to evaluate
c	covariances between local averages. (Jun 16/00)
c  1.2	use Gauss quadrature unless intervals overlap, eliminated lvarfn
c	(Apr 5/01)
c  1.3	now use Gauss quadrature for everything. (May 8/01)
c  1.4	now use 16 point Gauss quadrature (Mar 20/03)
c  1.5	simplification below only correct for volume completely enclosed.
c	Now using 20-point Gauss quadrature (Jul 2/03)
c---------------------------------------------------------------------------
      real*8 function dcvab3( cov, Dx, Dy, Dz, C1, C2, C3 )
      parameter (NG = 20)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)

      external cov
      common/dparam/ var, pb, dthx, dthy, dthz

      data zero/0.d0/, half/0.5d0/, one/1.d0/, onept5/1.5d0/, two/2.d0/
      data eighth/0.125d0/, fv12th/0.001953125d0/	!(= 1/512)
      data halfe/0.500000000001d0/			! half plus a wee bit
c			these are for NG = 3
c     data w/0.555555555555556d0,.888888888888889d0,.555555555555556d0/
c     data z/-.774596669241483d0,.000000000000000d0,.774596669241483d0/
c			and these are for NG = 5
c     data w/ .236926885056189d0,.478628670499366d0,.568888888888889d0,
c    >        .478628670499366d0,.236926885056189d0/
c     data z/-.906179845938664d0,-.538469310105683d0,.000000000000000d0,
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
c     data z/-.989400934991649932596d0, -.944575023073232576078d0,
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
      data z/-.993128599185094924786d0, -.963971927277913791268d0,
     >       -.912234428251325905868d0, -.839116971822218823395d0,
     >       -.746331906460150792614d0, -.636053680726515025453d0,
     >       -.510867001950827098004d0, -.373706088715419560673d0,
     >       -.227785851141645078080d0, -.076526521133497333755d0,
     >       0.076526521133497333755d0, 0.227785851141645078080d0,
     >       0.373706088715419560673d0, 0.510867001950827098004d0,
     >       0.636053680726515025453d0, 0.746331906460150792614d0,
     >       0.839116971822218823395d0, 0.912234428251325905868d0,
     >       0.963971927277913791268d0, 0.993128599185094924786d0/

c					enclosed volume: GQ simplifies
c					to variance of larger volume

      if( (C1.lt.halfe).and.(C2.lt.halfe).and.(C3.lt.halfe) ) then
         d1 = zero
         do 30 i = 1, NG
            xi = Dx*(one + z(i))
            d2 = zero
            do 20 j = 1, NG
               yj = Dy*(one + z(j))
               d3 = zero
               do 10 k = 1, NG
                  zk = Dz*(one + z(k))
                  d3 = d3 + w(k)*(one - z(k))*cov(xi,yj,zk)
  10           continue
               d2 = d2 + w(j)*(one-z(j))*d3
  20        continue
            d1 = d1 + w(i)*(one-z(i))*d2
  30     continue
         dcvab3 = eighth*d1
         return
      endif
c					otherwise, non-enclosed volume
      r1 = half*Dx
      r2 = half*Dy
      r3 = half*Dz

      s2 = two*C1
      s1 = s2 - two
      s3 = s2 + two
      v2 = two*C2
      v1 = v2 - two
      v3 = v2 + two
      u2 = two*C3
      u1 = u2 - two
      u3 = u2 + two

      d1 = zero
      do 60 i = 1, NG
         x1  = r1*(z(i) + s1)
         x2  = r1*(z(i) + s2)
         x3  = r1*(z(i) + s3)
         d21 = zero
         d22 = zero
         d23 = zero
         do 50 j = 1, NG
            y1 = r2*(z(j) + v1)
            y2 = r2*(z(j) + v2)
            y3 = r2*(z(j) + v3)
            d31 = zero
            d32 = zero
            d33 = zero
            d34 = zero
            d35 = zero
            d36 = zero
            d37 = zero
            d38 = zero
            d39 = zero
            do 40 k = 1, NG
               z1 = r3*(z(k) + u1)
               z2 = r3*(z(k) + u2)
               z3 = r3*(z(k) + u3)
               dp = one + z(k)
               dm = one - z(k)
               d31 = d31 + w(k)*(dp*cov(x1,y1,z1) + two*cov(x1,y1,z2)
     >                          +dm*cov(x1,y1,z3))
               d32 = d32 + w(k)*(dp*cov(x1,y2,z1) + two*cov(x1,y2,z2)
     >                          +dm*cov(x1,y2,z3))
               d33 = d33 + w(k)*(dp*cov(x1,y3,z1) + two*cov(x1,y3,z2)
     >                          +dm*cov(x1,y3,z3))
               d34 = d34 + w(k)*(dp*cov(x2,y1,z1) + two*cov(x2,y1,z2)
     >                          +dm*cov(x2,y1,z3))
               d35 = d35 + w(k)*(dp*cov(x2,y2,z1) + two*cov(x2,y2,z2)
     >                          +dm*cov(x2,y2,z3))
               d36 = d36 + w(k)*(dp*cov(x2,y3,z1) + two*cov(x2,y3,z2)
     >                          +dm*cov(x2,y3,z3))
               d37 = d37 + w(k)*(dp*cov(x3,y1,z1) + two*cov(x3,y1,z2)
     >                          +dm*cov(x3,y1,z3))
               d38 = d38 + w(k)*(dp*cov(x3,y2,z1) + two*cov(x3,y2,z2)
     >                          +dm*cov(x3,y2,z3))
               d39 = d39 + w(k)*(dp*cov(x3,y3,z1) + two*cov(x3,y3,z2)
     >                          +dm*cov(x3,y3,z3))
  40        continue
            dp = one + z(j)
            dm = one - z(j)
            d21 = d21 + w(j)*(dp*d31 + two*d32 + dm*d33)
            d22 = d22 + w(j)*(dp*d34 + two*d35 + dm*d36)
            d23 = d23 + w(j)*(dp*d37 + two*d38 + dm*d39)
  50     continue
         d1 = d1 + w(i)*((one+z(i))*d21 + two*d22 + (one-z(i))*d23)
  60  continue
      dcvab3 = fv12th*d1

      return
      end
