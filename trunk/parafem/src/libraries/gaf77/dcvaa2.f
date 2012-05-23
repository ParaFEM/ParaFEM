c  *********************************************************************
c  *                                                                   *
c  *                           function dcvaa2                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.5
c  Written by Gordon A. Fenton, DalTech, July 17, 1992
c  Latest Update: Jul 2, 2003
c
c  PURPOSE  returns the covariance between two 2-D local averages of equal
c           area. Used by LAS2G.
c
c  This function evaluates the covariance between two local averages in
c  2-dimensional space. The local averages are assumed to be of equal
c  size, Dx x Dy, and separated in space by the lags Tx=C1*Dx and Ty=C2*Dy.
c
c  The covariance is obtained by a 16-pt Gaussian quadrature of a
c  4-D integral collapsed (by assuming the covariance function to be
c  quadrant symmetric) to a 2-D integral.
c  NOTE: if the covariance function is not quadrant symmetric, this
c        routine will return erroneous results. Quadrant symmetry means
c        that cov(X,Y) = cov(-X,Y) = cov(X,-Y) = cov(-X,-Y), where
c        cov(.,.) is the function returning the covariance between
c        points separated by lag (X,Y), as discussed next.
c
c  The covariance function is referenced herein using a call of the
c  form
c
c          V = cov( X, Y )
c
c  where X and Y are the lag distances between the points in the field.
c
c  Arguments to this function are as follows
c
c    cov       external covariance function provided by the user.
c
c    Dx        x-dimension of each local average. (input)
c
c    Dy        y-dimension of each local average. (input)
c
c    C1        x-direction distance between local average centers is C1*Dx.
c              (input)
c
c    C2        y-direction distance between local average centers is C2*Dy.
c              (input)
c
c  REVISION HISTORY:
c  1.1	now including a Gaussian Quadrature integration option as an
c	alternative to the variance function approach to evaluate
c	covariances between local averages. (Jun 16/00)
c  1.2	now choose between Gauss Quadrature and var function approaches
c	to minimize numerical errors. Eliminated lvarfn. (Mar 27/01)
c  1.3	use Gauss quadrature for everything (Apr 11/01)
c  1.4	now use 16-pt Gauss quadrature for increased accuracy (Apr 8/02)
c  1.5	simplification below only correct for areas completely overlapping
c	(Jul 2/03)
c---------------------------------------------------------------------------
      real*8 function dcvaa2( cov, Dx, Dy, C1, C2 )
      parameter (NG = 16)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)

      external cov

      data zero/0.d0/, sx/0.0625d0/, quart/0.25d0/, half/0.5d0/
      data one/1.d0/, two/2.d0/

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

      r1  = half*Dx
      r2  = half*Dy
c					if areas overlap, GQ simplifies

      if( (C1 .eq. zero) .and. (C2 .eq. zero) ) then
         
         d1 = zero
         do 20 i = 1, NG
            xi = r1*(one + z(i))
            d2 = zero
            do 10 j = 1, NG
               yi = r2*(one + z(j))
               d2 = d2 + w(j)*(one-z(j))*cov(xi,yi)
  10        continue
            d1 = d1 + w(i)*(one-z(i))*d2
  20     continue
         dcvaa2 = quart*d1
         return

      endif
c					otherwise, non-overlapping areas
      s1 = two*C1 - one
      s2 = two*C1 + one
      v1 = two*C2 - one
      v2 = two*C2 + one
      d1 = zero
      do 40 i = 1, NG
         xi1 = r1*(z(i) + s1)
         xi2 = r1*(z(i) + s2)
         d21 = zero
         d22 = zero
         do 30 j = 1, NG
            yi1 = r2*(z(j) + v1)
            yi2 = r2*(z(j) + v2)
            dz1 = one + z(j)
            dz2 = one - z(j)
            d21 = d21 + w(j)*(dz1*cov(xi1,yi1) + dz2*cov(xi1,yi2))
            d22 = d22 + w(j)*(dz1*cov(xi2,yi1) + dz2*cov(xi2,yi2))
  30     continue
         d1 = d1 + w(i)*((one+z(i))*d21 + (one-z(i))*d22)
  40  continue
      dcvaa2 = sx*d1

      return
      end
