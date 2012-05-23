c  *********************************************************************
c  *                                                                   *
c  *                        subroutine gdate                           *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.0
c  Written by Fliegel, H. F. and van Flandern, T. C. (1968).
c  Communications of the ACM, Vol. 11, No. 10 (October, 1968).
c  Latest Update by Gordon A. Fenton: May 5, 2006
c
c  PURPOSE  computes the year, month, day of the Gregorian date given
c           its Julian day number.
c
c  DESCRIPTION
c  This function computes and returns the Gregorian date given the current
c  Julian Date.
c
c  The Julian date is the the number of days since noon on January 1,
c  -4712, i.e., January 1, 4713 BC (Seidelmann 1992). It was proposed by
c  Julius J. Scaliger in 1583, so the name for this system derived from
c  Julius Scaliger, not Julius Caesar. Scaliger defined Day One as a day
c  when three calendrical cycles converged.
c
c  This function assumes that the current dates are given according to
c  the Gregorian calendar, which was adopted in England and its colonies
c  on Sept 2, 1752, and which you are most likely using.
c
c  The algorithm was obtained from
c  http://aa.usno.navy.mil/faq/docs/JD_Formula.html
c
c  ARGUMENTS
c      jd	integer containing the Julian date (input)
c
c     iyr	integer containing the current year (e.g. 2006). (output)
c
c     imn	integer containing the month number, where Jan = 1 and
c		Dec = 12. (output)
c
c     idy	integer containing the day number (1 - 31). (output)
c
c  REVISION HISTORY:
c-------------------------------------------------------------------------
      subroutine gdate(jd,iyr,imn,idy)
      integer jd,iyr,imn,idy,i,j,l,n

      l   = jd+68569
      n   = 4*l/146097
      l   = l-(146097*n+3)/4
      i   = 4000*(l+1)/1461001
      l   = l-1461*i/4+31
      j   = 80*l/2447
      idy = l-2447*j/80
      l   = j/11
      imn = j+2-12*l
      iyr = 100*(n-49)+i+l

      return
      end
