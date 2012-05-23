c  *********************************************************************
c  *                                                                   *
c  *                         subroutine prfmte                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Jun  4, 1997
c  Latest Update: Oct 2, 2001
c
c  PURPOSE  prints a number using e format
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
c	k	output unit number to which the number is printed (without
c		concluding newline character). (input)
c
c  REVISION HISTORY:
c  1.1	corrected incorrect rendering of powers greater than 9. (Jul 28/97)
c  1.2	calling routine now provides iw.id format (Oct 2/01)
c-------------------------------------------------------------------------
      subroutine prfmte(val,iw,id,k)
      character fstr*256
      character d(10)*1
      data d/'0','1','2','3','4','5','6','7','8','9'/

   1  format(a,$)
c					transfer iw and id to temp vars
      jw = iw
      jd = id
c					get absolute value of val
      aval = abs(val)
c					check specification; iw > id + 6
      idm = jd
      iwm = jw
      iy  = 6
      if( val .lt. 0 ) iy = 7
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
         if( (idm .lt. 0) .and. (iwm .lt. 0 ) ) then
            write(k,1) '0.e+00'
            return
         endif
         do 10 i = jd+7, jw
            write(k,1) ' '
  10     continue
         write(k,1) '0.'
         do 20 i = 1, jd
            write(k,1) '0'
  20     continue
         write(k,1) 'e+00'
         return
      endif
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
c					pointer into fstr
      m = 1
c					set sign
      if( val .lt. 0.0 ) then
         fstr(m:m) = '-'
         m = m + 1
      endif
      t = aval
      ij = int( t )
      fstr(m:m) = d(1+ij)
      fstr(m+1:m+1) = '.'
      m = m + 2
      t = t - float( ij )
      mm = jd + m
  40  if( m .lt. mm ) then
         t  = t*10.
         ij = int( t )
         fstr(m:m) = d(1+ij)
         m  = m + 1
         t  = t - float( ij )
         go to 40
      endif
c					optionally eliminate trailing zeroes
      if( idm .lt. 0 ) then
         do 50 i = m-1, 1, -1
            if( fstr(i:i) .ne. '0' ) go to 60
            m  = m - 1
            jd = jd - 1
  50     continue
      endif
  60  fstr(m:m) = 'e'
      m = m + 1
      if( im .lt. 0 ) then
         fstr(m:m) = '-'
      else
         fstr(m:m) = '+'
      endif
      iam = iabs(im)
      if( iam .gt. 99 ) then		! this shouldn't happen in real*4
         if( iam .gt. 999 ) then
            fstr(m+1:m+3) = '***'
            m = m + 3
         else
            i1 = iam/100
            i2 = (iam - 100*i1)/10
            i3 = (iam - 100*i1 - 10*i2)
            m = m + 1
            fstr(m:m) = d(1+i1)
            m = m + 1
            fstr(m:m) = d(1+i2)
            m = m + 1
            fstr(m:m) = d(1+i3)
         endif
      elseif( iam .gt. 9 ) then
         i1 = iam/10
         i2 = iam - 10*i1
         m = m + 1
         fstr(m:m) = d(1+i1)
         m = m + 1
         fstr(m:m) = d(1+i2)
      else
         m = m + 1
         fstr(m:m) = '0'
         m = m + 1
         fstr(m:m) = d(1+iam)
      endif
c					output val
      if( (idm .lt. 0) .and. (iwm .lt. 0) ) jw = jd + iy
      do 70 i = jd+iy+1, jw
         write(k,1) ' '
  70  continue
      write(k,1) fstr(1:m)
c					all done
      return
      end
