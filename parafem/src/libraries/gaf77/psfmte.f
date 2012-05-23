c  *********************************************************************
c  *                                                                   *
c  *                         subroutine psfmte                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Wed Jun  4 22:44:03 1997
c  Latest Update: Oct 2, 2001
c
c  PURPOSE  writes a number using e format to an internal character string
c
c  This routine writes the real value `val' to an internal character string
c  using the format
c
c	d.ddde+ee
c
c  where `ddd' expands into as many digits as necessary, or as specified
c  by the iw and id arguments, to the left of the decimal place, and `ee'
c  represents a power of ten. For example, if the number 123.4567 is to be
c  printed using the format e12.5 (see the iw.id description below), then the
c  printed value will appear as
c
c	1.23456e+02
c
c  in scientific notation.
c
c  If both iw and id are negative, then the actual format used is e11.5 if
c  val is positive, and e12.5 if val is negative.
c
c  If either iw or id are non-negative, then the format used to represent the
c  number is one of the following, where the variable y takes values
c
c			y = 6 if val is non-negative
c			y = 7 if val is negative
c
c
c	iw.id	(both iw and id are non-negative) This specifies that a
c		field width of iw is to be used with id digits to the
c		right of the decimal place. iw must be at least (id+y)
c		to accomodate the leading digit, the decimal, and the
c		trailing 'e+dd' characters. If iw < (id+y) it is internally
c		set to (id+y). If iw is greater than (id+y), the number is
c		right justified with leading blanks.
c
c	iw	(id is negative) This specifies the field width only. The
c		number of decimal digits is determined as id = (iw-y), so
c		iw must be at least y (if it isn't, it is set to y). Since
c		this is a single precision routine, values of iw greater than
c		12 result in a 12 character right justified number
c		(that is, at most 7 significant digits are shown).
c
c	.id	(iw is negative) This specifies the number of digits to
c		display to the right of the decimal place. iw is computed
c		as (id+y). If id is greater than 6, it is set internally
c		to 6 to provide at most 7 significant digits.
c
c
c  Arguments to this routine are as follows;
c
c     val	real value containing the number to print. (input)
c
c	s	character string to which the formatted number is to be
c		written starting at the l'th character. (input/output)
c
c       l	integer which indexes the position in the string `s' at
c		which the current number is to start. On output, l points
c		at the character just beyond the written string. (input/output)
c
c      iw	integer containing the minimum total field width of the number
c		to be printed. If iw is insufficient to contain the number, it
c		is increased accordingly. If iw < 0, then it is assumed to not
c		be set and a minimum value is computed internally. (input)
c
c      id	integer containing the maximum number of digits to show to the
c		right of the decimal point. If id < 0, then it is assumed to
c		not be set and the minimum number of digits to the right of
c		the decimal point required to show the number is computed
c		internally. (input)
c
c  REVISION HISTORY:
c  1.1	calling routine now provides iw.id format (Oct 2/01)
c-------------------------------------------------------------------------
      subroutine psfmte(val,s,l,iw,id)
      character*(*) s
      character d(10)*1
      logical lround
      data d/'0','1','2','3','4','5','6','7','8','9'/
c					transfer iw and id to temp values
      jw = iw
      jd = id
c					get absolute value of val
      aval = abs(val)
c					check specification; iw > id + 6
      iy = 6
      if( val .lt. 0 ) iy = 7
c					adjust jw and jd if necessary
      if( (jw .gt. 0) .and. (jd .gt. 0) ) then
         if( jw .lt. jd+iy ) jw = jd + iy
      elseif( jw .gt. 0 ) then
         if( jw .lt. iy ) jw = iy
         jd = jw - iy
      elseif( jd .gt. 0 ) then
         if( jd .gt. 6 ) jd = 6
         jw = jd + iy
      else
         jw = 5 + iy
         jd = 5
      endif
c					special case: val = 0 ---------------
      if( val .eq. 0.0 ) then
         do 10 i = jd+7, jw
            s(l:l) = ' '
            l = l + 1
  10     continue
         s(l:l+1) = '0.'
         l = l + 2
         do 20 i = 1, jd
            s(l:l) = '0'
            l = l + 1
  20     continue
         s(l:l+3) = 'e+00'
         l = l + 4
         return
      endif
c					pointer into fstr
      m = 1
      lround = .false.
c					derive magnitude of val
  30  al = alog10(aval)
      im = int( al )
      if( (aval .lt. 1.0) .and. (aval .ne. 10.**im) ) im = im - 1
c					normalize
      aval = aval/(10.**im)
c					round val
      aval = aval + 0.5*10.**(-jd)
c					recheck val to see if it has changed
      ial = int( alog10(aval) )
      if( ial .ne. 0 ) then
         im = im + 1
         aval = aval/10.
      endif
c					output val
      do 40 i = jd+iy+1, jw
         s(l:l) = ' '
         l = l + 1
  40  continue
c					set sign
      if( val .lt. 0.0 ) then
         s(l:l) = '-'
         l = l + 1
      endif
      t = aval
      ij = int( t )
      s(l:l) = d(1+ij)
      s(l+1:l+1) = '.'
      l = l + 2
      t = t - float( ij )
      m = 1
  50  if( m .le. jd ) then
         t  = t*10.
         ij = int( t )
         s(l:l) = d(1+ij)
         l = l + 1
         m = m + 1
         t = t - float( ij )
         go to 50
      endif
      s(l:l) = 'e'
      l = l + 1
      if( im .lt. 0 ) then
         s(l:l) = '-'
      else
         s(l:l) = '+'
      endif
      l = l + 1
      iam = iabs(im)
      ij  = iam/10
      s(l:l) = d(1+ij)
      l = l + 1
      iam = iam - ij*10
      s(l:l) = d(1+iam)
      l = l + 1
c					all done
      return
      end
