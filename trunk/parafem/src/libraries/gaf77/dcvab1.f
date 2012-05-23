c  *********************************************************************
c  *                                                                   *
c  *                           function dcvab1                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.4
c  Written by Gordon A. Fenton, DalTech, Jun 9, 2000
c  Latest Update: Apr 8, 2001
c
c  PURPOSE  returns the covariance between two 1-D local averages, one of
c           twice the length of the other. Used by LAS1G.
c
c  This function evaluates the covariance between two local averages in
c  1-dimensional space. The local averages are assumed to be of size Dx
c  and 2*Dx, and separated in space by the lag Tx=C1*Dx, center-to-center.
c
c  The covariance is obtained by a 5-pt Gaussian quadrature of a
c  1-dimensional integral of the externally provided point covariance
c  function of the process
c
c                   1    (C1-0.5)*Dx
c  Cov[X_a, X_b] = --- [ int  {x - (C1-1.5)Dx}*cov(x) dx
c                 2Dx^2 (C1-1.5)Dx
c
c			(C1+0.5)Dx
c			+ int Dx*cov(x) dx
c			(C1-0.5)*Dx
c
c			(C1+1.5)Dx
c			+ int {(C1+1.5) - Dx}*cov(x) dx ]
c			(C1+0.5)*Dx
c
c
c  Where X_a and X_b are the local averages of the two areas, of size
c  Dx and 2*DX, respectively, and separated center-to-center by lag C1*Dx.
c
c  The covariance function is referenced herein using a call of the form
c
c          V = cov( X )
c
c  where X is the lag distances between the points in the field.
c
c  Arguments to this function are as follows
c
c    cov       external variance/covariance function provided by the user.
c
c    Dx        x-dimension of the smaller local average. (input)
c
c    C1        x-direction distance between local average centers is C1*Dx.
c              (input)
c
c  REVISION HISTORY:
c  1.1	eliminated lvarfn - now choose integration method to minimize
c	errors. (Feb 27/01)
c  1.2	now choose Gauss Quadrature except when intervals overlap. (Apr 5/01)
c  1.3	now use Gauss Quadrature for everything (Apr 16/01)
c  1.4	now use 16-pt Gauss quadrature for increased accuracy (Apr 8/02)
c---------------------------------------------------------------------------
      real*8 function dcvab1( cov, Dx, C1 )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)

      external cov

      data zero/0.d0/, eighth/0.125d0/, half/0.5d0/
      data one/1.d0/, onept5/1.5d0/, two/2.d0/

c			these are for NG = 3
c     data w/0.555555555555556d0,.888888888888889d0,.555555555555556d0/
c     data z/-.774596669241483d0,.000000000000000d0,.774596669241483d0/
c			and these are for NG = 5
c     data w/ .236926885056189d0,.478628670499366d0,.568888888888889d0,
c    >        .478628670499366d0,.236926885056189d0/
c     data z/-.906179845938664d0,-.538469310105683d0,.000000000000000d0,
c    >        .538469310105683d0, .906179845938664d0/
c			these are for NG = 16
      data w/0.027152459411754094852d0, 0.062253523938647892863d0,
     >       0.095158511682492784810d0, 0.124628971255533872052d0,
     >       0.149595988816576732081d0, 0.169156519395002538189d0,
     >       0.182603415044923588867d0, 0.189450610455068496285d0,
     >       0.189450610455068496285d0, 0.182603415044923588867d0,
     >       0.169156519395002538189d0, 0.149595988816576732081d0,
     >       0.124628971255533872052d0, 0.095158511682492784810d0,
     >       0.062253523938647892863d0, 0.027152459411754094852d0/
      data z/-.989400934991649932596d0, -.944575023073232576078d0,
     >       -.865631202387831743880d0, -.755404408355003033895d0,
     >       -.617876244402643748447d0, -.458016777657227386342d0,
     >       -.281603550779258913230d0, -.095012509837637440185d0,
     >       0.095012509837637440185d0, 0.281603550779258913230d0,
     >       0.458016777657227386342d0, 0.617876244402643748447d0,
     >       0.755404408355003033895d0, 0.865631202387831743880d0,
     >       0.944575023073232576078d0, 0.989400934991649932596d0/

      abs(zz) = dabs(zz)
c					if intervals overlap, use special form
      if( abs(C1) .lt. one ) then	! NOTE: this really just means C1 = 0.5
         d1 = zero
         do 10 i = 1, NG
            xi = Dx*(one + z(i))
            d1 = d1 + w(i)*(one-z(i))*cov(xi)
  10     continue
         dcvab1 = half*d1
         return
      endif
c					otherwise, use GQ for non-overlapping
c					intervals
      r1 = half*Dx
      s2 = two*C1
      s1 = s2 - two
      s3 = s2 + two
      d1 = zero
      do 20 i = 1, NG
         xi1 = r1*(z(i) + s1)
         xi2 = r1*(z(i) + s2)
         xi3 = r1*(z(i) + s3)
         q1  = (one + z(i))*cov(xi1)
         q2  = two*cov(xi2)
         q3  = (one - z(i))*cov(xi3)
         d1  = d1 + w(i)*(q1 + q2 + q3)
  20  continue
      dcvab1 = eighth*d1

      return
      end
