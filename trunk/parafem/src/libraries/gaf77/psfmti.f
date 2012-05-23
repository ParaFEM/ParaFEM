c  *********************************************************************
c  *                                                                   *
c  *                         subroutine psfmti                         *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.1
c  Written by Gordon A. Fenton, TUNS, Sun Jun  8 21:34:21 1997
c  Latest Update: Oct 2, 2001
c
c  PURPOSE  writes an integer to an internal character string using
c           an optional field width
c
c  This routine writes the integer `ival' to an internal character string
c  using an optional field width specification.
c
c  Arguments to this routine are as follows;
c
c    ival	integer containing the number to print. (input)
c
c	s	character string to which the formatted number is to be
c		written starting at the l'th character. (input/output)
c
c       l	integer which indexes the position in the string `s' at
c		which the current number is to start. On output, l points
c		at the character just beyond the written string. (input/output)
c
c      iw	integer containing the desired field width. If iw < 0
c		then the minimum field width is computed internally (note
c		that integers larger than 256 digits cannot be printed).
c		If iw is larger than required, the number is right
c		justified. (input)
c
c  REVISION HISTORY:
c  1.1	calling routine now provides the iw format (Oct 2/01)
c-------------------------------------------------------------------------
      subroutine psfmti(ival,s,l,iw)
      character*(*) s
      character fstr*256, d(10)*1
c					basic digits
      data d/'0','1','2','3','4','5','6','7','8','9'/

c					get absolute value of val
      iaval = iabs(ival)
c					build character string of number
      m = 0
      n = 257				! work backwards in fstr
  20  m = m + 1
      n = n - 1
      i = iaval/10
      ii = 10*i
      j = iaval - ii
      fstr(n:n) = d(j+1)
      iaval = i
      if( iaval .gt. 0 ) go to 20
c					add `-' if negative
      if( ival .lt. 0 ) then
         m = m + 1
         n = n - 1
         fstr(n:n) = '-'
      endif
c					now print it
      do 30 i = m+1, iw
         s(l:l) = ' '
         l = l + 1
  30  continue
      ll = 256 - n + 1
      s(l:l+ll-1) = fstr(n:256)
      l = l + ll

      return
      end
