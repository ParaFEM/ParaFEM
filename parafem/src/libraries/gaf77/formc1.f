c  *********************************************************************
c  *                                                                   *
c  *                         subroutine formc1                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Thu Jul 18 11:09:26 1996
c
c  PURPOSE  forms a covariance matrix for a 1-D process
c
c  This routine forms a Toeplitz covariance matrix for a stationary
c  1-D point process, where the points are equispaced. A variety of
c  correlation functions are supported;
c
c	'Gauss-Markov'	for an exponentially decaying correlation
c			of the form exp{-2*x/scale}, where `x' is
c			the physical distance between the two points,
c			and `scale' is the scale of fluctuation.
c
c	'Gaussian'	for a correlation function which decays in the
c			same was as a Gaussian distribution:
c			exp{-pi*(x/scale)^2}
c
c	'Fractal'	for a fractal correlation function with Hurst
c			parameter, 1/2 < H < 1. H = 1/2 gives a set
c			of independent samples, while H = 1 corresponds
c			to complete correlation over all lags.
c
c	'Simple'	for a simple correlation function of the form
c			[scale/(scale + x)]^3
c
c  Arguments to this routine are as follows;
c
c	C	real array of size at least n x n, which on output will
c		contain the desired covariance matrix. (output)
c
c      ic	leading dimension of the array C as prescribed in the
c		calling routine. (input)
c
c	n	the number of points in the 1-D process, and the size of
c		the array C. (input)
c
c     var	the process point variance. (input)
c
c   scale	the scale of fluctuation (or Hurst parameter) of the 1-D
c		process. (input)
c
c      dx	the spacing between points in the process. (input)
c
c    type	character string containing the name of the correlation
c		function to use in forming C. Possible values of type are;
c			'Gauss-Markov'
c			'Gaussian'
c			'Fractal'
c			'Simple'
c		See discussion above for descriptions of these types. (input)
c
c    ierr	integer error flag which is zero if all goes well. Ierr is
c		set to -1 and an error message issued to standard output if
c		`type' is not one of the above, or if type = 'Fractal' and
c		`scale' lies outside the interval [0.5, 1]. (output)
c-------------------------------------------------------------------------
      subroutine formc1( C, ic, n, var, scale, dx, type, ierr)
      implicit real*8 (a-h,o-z)
      dimension C(ic,*)
      character*(*) type
      data zero/0.d0/, half/0.5d0/, one/1.d0/, two/2.d0/
      data pi/3.1415926535897932384d0/
      exp(z)   = dexp(z)
      float(i) = dble(i)

   1  format(5a)
c					assume the worst
      ierr = -1
c					form essential elements of Toeplitz C
      C(1,1) = var
      if( type(1:12) .eq. 'Gauss-Markov' ) then
         div = -two/scale
         do 10 i = 2, n
            x = div*float(i-1)*dx
            C(i,1) = var*exp(x)
  10     continue
      elseif( type(1:8) .eq. 'Gaussian' ) then
         div = -pi/(scale*scale)
         do 20 i = 2, n
            x = float(i-1)*dx
            C(i,1) = var*exp(div*x*x)
  20     continue
      elseif( type(1:7) .eq. 'Fractal' ) then
         if( (scale .lt. half) .or. (scale .gt. one) ) then
            write(6,1)'formc1: Can''t form covariance matrix for fractal
     > process with H outside interval [0.5, 1].'
            return
         endif
         if( scale .eq. one ) then
            do 30 i = 2, n
               C(i,1) = var
  30        continue
         elseif( scale .eq. half ) then
            do 40 i = 2, n
               C(i,1) = zero
  40        continue
         else
            th  = two*scale
            div = var/(two*dx**th)
            xm  = zero
            x   = dx
            do 50 i = 2, n
               xp = float(i)*dx
               C(i,1) = div*(xp**th - two*x**th + xm**th)
               xm = x
               x  = xp
  50        continue
         endif
      elseif( type(1:6) .eq. 'Simple' ) then
         do 60 i = 2, n
            x = float(i-1)*dx
            th = scale/(scale + x)
            C(i,1) = var*th*th*th
  60     continue
      else
         lt = lnblnk(type)
         write(6,1)'formc1: Unknown correlation function type: ',
     >              type(1:lt)
         write(6,1)'        Known types: Gauss-Markov, Gaussian, Fractal
     >, Simple.'
         return
      endif
c					fill in remainder of C
      do 70 i = 2, n
         C(i,i) = C(1,1)
         C(1,i) = C(i,1)
  70  continue
      do 80 i = 2, n-1
      do 80 j = i+1, n
         C(j,j-i+1) = C(i,1)
         C(j-i+1,j) = C(i,1)
  80  continue
c					all done
      ierr = 0
      return
      end
