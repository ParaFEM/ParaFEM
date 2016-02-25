c  *********************************************************************
c  *                                                                   *
c  *                         subroutine prints                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.61
c  Written by Gordon A. Fenton, TUNS, Sat Sep  9 13:21:15 1995
c  Latest Update: Jun 25, 2003
c
c  PURPOSE  prints real/integer values according to a C style format to
c           an internal character string.
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, along with a vector of additional
c  arguments and prints them to the given string. Normally, this routine is
c  called by one of prins1, prins2, ... See their headers for usage
c  instructions.
c  Arguments to this routine are as follows;
c
c       s	character string to which output is directed. (input)
c
c     fmt	format string. (input)
c
c     val	real vector of length at least n containing the elements
c		to be printed. An element of val can be either real or integer,
c		it's type is determined by the presence of a corresponding
c	 	%f or %e (for real) or %i (for integer) descriptor. See notes
c		below on the optional format specification which may proceed
c		any of these format descriptors.
c
c       n	the number of elements in val. Currently this must be less
c		than 256. (input)
c
c   The following % formats are recognized by this routine; only the %f, %e,
c   and %i formats actually consume an argument (ie, make use of one of
c   the val? arguments).
c
c	   %<space> to produce a space
c
c	   %%	to produce a % character
c
c	   %b	to produce a '\' character
c
c	   %t	to produce a tab character
c
c	   %n	to produce a newline character
c
c	   %f	to format a real value in the form '12.236'. If the value
c		is outside the range (1.e-6, 1.e+7) then the %e format is
c		automatically used. The %f format can be optionally written
c		in the following forms;
c
c		%w.df	specifies that a field width of `w' is to be used
c			with d digits to the right of the decimal place.
c			The decimal place takes up one character in the
c			total field width. `w' is actually a minimum width
c			and will be increased if more digits are required to
c			the left of the decimal point at the expense of
c			digits to the right of the decimal, if any. If the
c			number to be printed occupies less than w digits,
c			it is right justified. Digits to the right of the
c			7'th significant digit are left blank. For example,
c			%6.2f implies that the number 12.345 will be
c			printed as ' 12.35'.
c
c		%wf	specifies a minimum field width. The number of
c			decimal places is determined as necessary. Smaller
c			numbers are right justified. At most 7 significant
c			digits will be displayed. For example, %8f implies
c			that the number 12.34500 will be printed as
c			'  12.345'.
c
c		%.df	specifies a maximum number of digits to display
c			to the right of the decimal place.
c			Trailing zeroes are not shown. For example, %.2f
c			implies that the number 12.345 will be printed
c			as '12.35'.
c
c		%.f	at the moment, this is only specifies that zero is
c			printed as '0.'. The default is '0'.
c
c	   %e	to format a real value in scientific notion. For example
c		the number 12.345 is written as 1.2345e+01 in scientific
c		notation (where the e+01 means 'times 10 to the power 1').
c		The %e format can be optionally written in the following
c		forms, where;
c
c			y = 6 if the value being printed is non-negative
c			y = 7 if the value being printed is negative
c
c		%w.de	specifies that a field width of `w' is to be used
c			with d digits to the right of the decimal place.
c			w must be at least (d+y) to accomodate the leading
c			digit, the decimal, and the trailing 'e+dd' characters.
c			If w < (d+y) it is internally set to (d+y). If w
c			is greater than (d+y), the number is right justified
c			with leading blanks. For example, %10.2e implies that
c			the number 12.345 will be printed as '  1.23e+01'.
c
c		%we	specifies the field width only. The number of decimal
c			digits is determined as d = (w-y), so w must be at
c			least y (if it isn't, it is set to y). Since this is
c			a single precision routine, values of w greater than
c			12 result in a 12 character right justified number
c			(at most 7 significant digits are shown). For example,
c			%10e implies that the number 12.345 will be printed
c			as '1.2340e+01'.
c
c		%.de	specifies the number of digits to display to the right
c			of the decimal place. w is computed as (d+y). If d is
c			greater than 6, it is set internally to 6 to provide
c			at most 7 significant digits. For example, %.2e implies
c			that the number 12.345 will be printed as '1.23e+01'.
c
c	   %i	to format an integer value. Normally, the minimum field width
c		required to represent the number is used. However, the %i
c		format can be optionally written in the following form;
c
c		%wi	specifies that a field width of at least 'w' is used.
c			For example, %10i implies that the number 12 will be
c			printed as '        12'.
c
c
c  REVISION HISTORY:
c  1.1	added %fw specification, where w is interpreted as the desired max.
c	number of significant digits to show on real output (Oct. 9, 1995)
c  1.2	revised %i, %f, and %e output formatting (Jun 9/97)
c  1.21	made lnblnk external to ensure use of GAF77's lnblnk. (Jul 19/99)
c  1.3	replaced \\ with char(92) for portability (Dec 1/99)
c  1.4	added blanks to end of string s (Jan 25/00)
c  1.5	replaced '\n' and ' ' with char(10) and char(32) for portability
c	(Sep 24/01)
c  1.6	replaced all \ escapes with % escapes (for compatibility with
c	Win98's Visual Fortran). '\' is now treated as an ordinary
c	character. (Oct 1/01)
c  1.61	now properly handle '% ' at end of format string (in fact '%' at
c	end of format string will always produce a space now). (Jun 25/03)
c-------------------------------------------------------------------------
      subroutine prints( s, fmt, val, n )
      parameter (MPX = 256)
      character*(*) fmt, s
      character*1 bslash, NL, CR, spc, tab
      real val(*), tmp(MPX)
      integer itmp(MPX)
      logical ldot, ldigit, LDOS
      equivalence (itmp(1),tmp(1))
      external lnblnk

   1  format(3a)

      LDOS   = .false.			! set this to .true. for Windows env
      tab    = char(9)
      NL     = char(10)			! also known as LF
      CR     = char(13)
      spc    = char(32)
      bslash = char(92)
      lmax   = len(s)
      MP     = n
      if( MP .gt. MPX ) MP = MPX
      do 5 i = 1, MP
         tmp(i) = val(i)
   5  continue

      jt = 1
      lf = lnblnk(fmt)
      i  = 0
      l  = 1
  10  i = i + 1
      if( i .gt. lf ) go to 20
      if( fmt(i:i) .eq. '%' ) then				! DESCRIPTOR
         i = i + 1
         if( i .gt. lf ) then
            s(l:l) = spc
            l = l + 1
            go to 20
         endif
         if( fmt(i:i) .eq. 'b' ) then
            s(l:l) = bslash
            l = l + 1
            if( l .gt. lmax ) go to 99
            go to 10
         elseif( fmt(i:i) .eq. 'n' ) then
            if( LDOS ) then
               s(l:l) = CR
               l = l + 1
            endif
            s(l:l) = NL
            l = l + 1
            if( l .gt. lmax ) go to 99
            go to 10
         elseif( fmt(i:i) .eq. 't' ) then
            s(l:l) = tab
            l = l + 1
            if( l .gt. lmax ) go to 99
            go to 10
         elseif( fmt(i:i) .eq. '%' ) then
            s(l:l) = '%'
            l = l + 1
            if( l .gt. lmax ) go to 99
            go to 10
         elseif( fmt(i:i) .eq. spc ) then
            s(l:l) = spc
            l = l + 1
            if( l .gt. lmax ) go to 99
            go to 10
         elseif( ldigit(fmt(i:i)) ) then
            call getfsp(fmt,lf,i,j,iw,id,ldot)
            if( fmt(j:j) .eq. 'f' ) then
               call psfmtf(tmp(jt),s,l,iw,id,ldot)
            elseif( fmt(j:j) .eq. 'e' ) then
               call psfmte(tmp(jt),s,l,iw,id)
            elseif( fmt(j:j) .eq. 'i' ) then
               call psfmti(tmp(jt),s,l,iw)
            else
               k = j - i + 1
               s(l:l+k) = fmt(i-1:j)
               l = l + k + 1
               if( l .gt. lmax ) go to 99
               i = j
               go to 10
            endif
            i = j
         elseif( fmt(i:i) .eq. 'f' ) then
            call psfmtf(tmp(jt),s,l,-1,-1,.false.)
         elseif( fmt(i:i) .eq. 'e' ) then
            call psfmte(tmp(jt),s,l,-1,-1)
         elseif( fmt(i:i) .eq. 'i' ) then
            call psfmti(tmp(jt),s,l,-1)
         else
            s(l:l+1) = fmt(i-1:i)
            l = l + 2
            if( l .gt. lmax ) go to 99
            go to 10
         endif
         if( l .gt. lmax ) go to 99
         jt = jt + 1
         if( jt .gt. MP ) jt = MP
         go to 10
      else							! CHARACTER
         s(l:l) = fmt(i:i)
         l = l + 1
         if( l .gt. lmax ) go to 99
      endif
      go to 10
c						add blanks to end of s
  20  do 30 j = l, lmax
         s(j:j) = spc
  30  continue
      return

  99  if( (l-1) .le. lmax ) return
      write(6,1)'*** Warning: GAFlib prints utility wrote past the',
     >                      ' end of an internal character string.'
      write(6,1)'             Data corruption may have occurred.'
      write(6,1)'             Lately writing the following format:'
      write(6,1) fmt(1:lf)
      return
      end
