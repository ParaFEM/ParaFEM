c  *********************************************************************
c  *                                                                   *
c  *                      real*8 function dfdinv                       *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Dec  2, 1998
c  Latest Update: Dec  2, 1998
c
c  PURPOSE  returns the inverse F-distribution function for given probability
c           p and degrees-of-freedom v1 and v2.
c
c  DESCRIPTION
c  This function computes the point x such that P[ F < x ] = p, where F is
c  F-distributed with degrees-of-freedom v1 and v2. The inverse is computed
c  using the inverse of the Beta distribution, since
c
c	P[ F < x ] = 1 - beta(y,a,b)
c
c  where y = (v2/(v2 + v1*x)), a = v2/2, and b = v1/2.
c  For large degrees-of-freedoms, v1 or v2, three asymptotic approximations
c  exist;
c
c	1) when v2 =infinity >> v1, P[ F < x ] is given by a Chi-Squared
c		distribution
c		P[ C < y ] with a = v1/2, b = 2, and y = v1*x
c	   I invoked this whenever v2 > 500 + 10*v1
c
c	2) when v1 = infty >> v2, P[ F < x ] is given by a Chi-Squared
c		distribution
c		1 - P[ C < y ] with a = v2/2, b = 2, and y = v2/x
c          I invoke this whenever v1 > 500 + 10*v2
c
c	3) when both v1 and v2 are large, P[ F < x ] is approximated by
c
c			phi( (x - a1)/(a1*a2) )
c
c		where a1 = v2/(v2-2) and
c                     a2 = sqrt( 2*(v1+v2-2)/(v1*(v2-4)) )
c          I invoke this when (v1+v2) > 500, and (v1,v2) do NOT satisfy
c	   items (1) and (2)
c
c  I have used 500 as the basic boundary between large and small. This still
c  results in a small `jump' as one moves from one approximation to
c  another.
c
c  ARGUMENTS
c	p	real value giving the desired probability level such that
c		p = P[ F < x], where x is returned. (input)
c
c      v1	real value containing the first number of degrees-of-freedom.
c		(input)
c
c      v2	real value containing the second number of degrees-of-freedom.
c		(input)
c
c    ierr	integer flag which is zero if the inverse is successfully
c		computed, -1 if not. (output)
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      real*8 function dfdinv(p,v1,v2,ierr)
      implicit real*8 (a-h,o-z)
      data half/0.5d0/, one/1.d0/, two/2.d0/, four/4.d0/
      data big/500.d0/, s/10.d0/

c	NOTE: Linux's f77 compiler doesn't pass args in the following by
c	      address... thus `i' can't be returned in `ierr'...
c	      we have to specifically call dbetnv below.
c     betnv(z1,z2,z3,i) = dbetnv(z1,z2,z3,i)
      sqrt(z) = dsqrt(z)
      phinv(z) = dphinv(z)


      ierr = 0
c					check for limiting cases in v1,v2
      if( v1+v2 .le. big ) then
         y = dbetnv( one-p, half*v2, half*v1, ierr )
         if( ierr .ne. 0 ) return
         dfdinv = (v2*(one-y))/(v1*y)
      elseif( v2 .gt. (big + s*v1) ) then
         dfdinv = dgminv(p,half*v1,two,ierr)/v1
      elseif( v1 .gt. (big + s*v2) ) then
         dfdinv = v2/dgminv(one-p,half*v2,two,ierr)
      else
         a1 = v2/(v2-two)
         a2 = sqrt( two*(v1+v2-two)/(v1*(v2-four)) )
         dfdinv = a1*(one + a2*phinv(p))
      endif

      return
      end
