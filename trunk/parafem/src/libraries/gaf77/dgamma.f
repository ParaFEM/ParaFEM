c  *******************************************************************
c  *                                                                 *
c  *                          function dgamma                        *
c  *                                                                 *
c  *******************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, 1991
c
c  PURPOSE  to return the Gamma (factorial) function
c
c  This routine calculates the definite integral
c
c                  inf
c    dgamma(z) = int  t**(z-1)*exp(-t) dt
c                  0
c
c  which is the Gamma function, for argument z. The method used
c  is a series expansion followed by a recursion relationship. The
c  coefficients used in the series expansion are obtained from Abramowitz
c  and Stegun ("Handbook of Mathematical Functions", US Dept. of Commerce,
c  10th printing, 1972), pg 256. The result should have a relative accuracy
c  of better than 1.e-14. If gamma(zz) would result in a value larger than
c  the maximum double precision number, "big", then "big" is returned. This
c  occurs also if zz is a negative integer.
c  Arguments are as follows;
c
c     zz    real value for which the gamma function is desired. (input)
c
c  NOTES: 1) "eps" is the machine epsilon (smallest number which can be added
c            to 1.0 without being lost)
c         2) "big" is the biggest floating point number
c         3) set IDEBUG = 1 if error messages should be reported
c--------------------------------------------------------------------------
      real*8 function dgamma( zz )
      implicit real*8 (a-h,o-z)
      parameter (IDEBUG = 1)
      dimension c(26)
      data c/1.0000000000000000d0,  0.5772156649015329d0,
     >      -0.6558780715202538d0, -0.0420026350340952d0,
     >       0.1665386113822915d0, -0.0421977345555443d0,
     >      -0.0096219715278770d0,  0.0072189432466630d0,
     >      -0.0011651675918591d0, -0.0002152416741149d0,
     >       0.0001280502823882d0, -0.0000201348547807d0,
     >      -0.0000012504934821d0,  0.0000011330272320d0,
     >      -0.0000002056338417d0,  0.0000000061160950d0,
     >       0.0000000050020075d0, -0.0000000011812746d0,
     >       0.0000000001043427d0,  0.0000000000077823d0,
     >      -0.0000000000036968d0,  0.0000000000005100d0,
     >      -0.0000000000000206d0, -0.0000000000000054d0,
     >       0.0000000000000014d0,  0.0000000000000001d0/
c					eps is the machine epsilon
c					big is the largest possible number
      data zero/0.d0/, eps/0.222045d-15/, one/1.d0/, big/1.797693d+308/
      data pi/3.141592653589793238462643d0/
      float(i) = dble(i)
      int(x)   = idint(x)
      abs(x)   = dabs(x)
      sin(x)   = dsin(x)

   1  format(a,e13.6,a)
c----------------------------- start executable statements
      z = abs(zz)

      q = z + eps
      n = int( q )
      f = z - float(n)
c					is z an integer?
      if( f .lt. q*eps ) then
         if( zz .le. zero ) go to 20
         n = n - 1
         f = one
         a = one
      else
c					series expansion for f < 1.0
         a = f*(c(22) + f*(c(23) + f*(c(24) + f*(c(25) + f*c(26)))))
         a = f*(c(17) + f*(c(18) + f*(c(19) + f*(c(20) + f*(c(21)+a)))))
         a = f*(c(12) + f*(c(13) + f*(c(14) + f*(c(15) + f*(c(16)+a)))))
         a = f*(c( 7) + f*(c( 8) + f*(c( 9) + f*(c(10) + f*(c(11)+a)))))
         a = f*(c( 2) + f*(c( 3) + f*(c( 4) + f*(c( 5) + f*(c( 6)+a)))))
         a = one/(f*(c( 1) + a))
      endif
c					apply recursion and watch for overflow
      if( n .gt. 0 ) then
         a = a*f
         do 10 i = 1, n-1
            q = f + float(i)
            if( a .gt. (big/q) ) go to 20
            a = a*q
  10     continue
      endif

      if( zz .lt. zero ) then
         dgamma = -pi/(z*a*sin(pi*z))
      else
         dgamma = a
      endif

      return
c				overflow error
  20  if( IDEBUG .eq. 1 )
     >   write(0,1)'Overflow calculating Gamma(',zz,') in DGAMMA.'
      dgamma = big
      return
      end
