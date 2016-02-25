c  *********************************************************************
c  *                                                                   *
c  *                          function fdist                           *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, DalTech, Dec  2, 1998
c  Latest Update: Dec  2, 1998
c
c  PURPOSE  returns the F-distribution function for given x and
c           degrees-of-freedom v1 and v2.
c
c  DESCRIPTION
c  This function computes P[ F < x ], where F is F-distributed with
c  degrees-of-freedom v1 and v2. It does so by recognizing that this
c  probability is a special case of the Beta distribution, namely
c
c	P[ F < x ] = 1 - beta(y,a,b)
c
c  where y = (v2/(v2 + v1*x)), a = v2/2, and b = v1/2.
c
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
c	x	real value giving the point at which the cumulative
c		distribution is desired. (input)
c
c      v1	real value containing the first number of degrees-of-freedom.
c		(input)
c
c      v2	real value containing the second number of degrees-of-freedom.
c		(input)
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      real*4 function fdist(x,v1,v2)
      data half/0.5/, one/1.0/, two/2.0/, four/4.0/
      data big/500.0/, s/10.0/
c					check for limiting cases in v1,v2
      if( v1+v2 .le. big ) then
         fdist = one - beta((v2/(v2+v1*x)), half*v2, half*v1)
      elseif( v2 .gt. (big + s*v1) ) then
         fdist = gmdst( v1*x, half*v1, two )
      elseif( v1 .gt. (big + s*v2) ) then
         fdist = one - gmdst( v2/x, half*v2, two )
      else
         a1 = v2/(v2-two)
         a2 = sqrt( two*(v1+v2-two)/(v1*(v2-four)) )
         fdist = phi( (x - a1)/(a1*a2) )
      endif

      return
      end
