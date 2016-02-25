c  *********************************************************************
c  *                                                                   *
c  *                          function julian                          *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.11
c  Written by Gordon A. Fenton, Dalhousie University, Apr 21, 2006
c  Latest Update: May 5, 2006
c
c  PURPOSE  computes the Julian day number given the year, month, and day
c
c  DESCRIPTION
c  This function computes and returns the Julian date given the current
c  year, month, and day of the month. It rounds to the nearest whole day
c  assuming that the time is around noon in the east of North America (this
c  corresponds to afternoon UT, so we add one day to the following
c  algorithm.
c
c  The Julian date is the the number of days since noon on January 1,
c  -4712, i.e., January 1, 4713 BC (Seidelmann 1992). It was proposed by
c  Julius J. Scaliger in 1583, so the name for this system derived from
c  Julius Scaliger, not Julius Caesar. Scaliger defined Day One as a day
c  when three calendrical cycles converged.
c
c  This function assumes that the current date is given according to
c  the Gregorian calendar, which was adopted in England and its colonies
c  on Sept 2, 1752, and which you are most likely using. The returned
c  Julian Date is valid for all AD dates in the Gregorian calendar.
c
c  The algorithm was obtained from
c  http://scienceworld.wolfram.com/astronomy/JulianDate.html
c
c  ARGUMENTS
c     iyr	integer containing the current year (e.g. 2006). (input)
c
c     imn	integer containing the month number, where Jan = 1 and
c		Dec = 12. (input)
c
c     idy	integer containing the day number (1 - 31). (input)
c
c  REVISION HISTORY:
c  1.1	replaced the algorithm http://www.jguru.com/faq/view.jsp?EID=14092
c  by Robert Baruch (May 29, 2000), which is broken. (Apr 27/06)
c  1.11	added comment that this is valid for all AD dates (May 5/06)
c-------------------------------------------------------------------------
      integer function julian(iyr,imn,idy)

      i1 = 7.*(iyr + int(float(imn+9)/12.))/4.
      i2 = 3.*(int((iyr+float(imn-9)/7.)/100.)+1.)/4.
      i3 = 275.*imn/9.
      i4 = 1721028 + idy + 1
      julian = 367*iyr - i1 - i2 + i3 + i4

      return
      end
