c  *********************************************************************
c  *                                                                   *
c  *                         subroutine printv                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.61
c  Written by Gordon A. Fenton, TUNS, Sat Sep  9 13:21:15 1995
c  Latest Update: Jun 25, 2003
c
c  PURPOSE  prints real/integer values according to a C style format
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, along with a vector of additional
c  arguments and prints them to the given unit. Normally, this routine is
c  called by one of print1, print2, ... See their headers for usage
c  instructions.
c
c  Arguments to this routine are as follows;
c
c       k	unit number to which output is directed. (input)
c
c     fmt	format string. (input)
c
c     val	real vector of length at least MP containing the elements
c		to be printed. An element of val can be either real or integer,
c		it's type is determined by the presence of a corresponding
c	 	%f or %e (for real) or %i (for integer) descriptor. See notes
c		below on the optional format specification which may follow
c		any of these format descriptors.
c		(input)
c
c       n	the number of elements in val. Currently this must be less
c		than 256. (input)
c
c   The following % formats are recognized by this routine;
c
c	   %<space> to produce a space
c
c	   %%	to produce a % character
c
c	   %b	to produce a '\' character
c
c	   %t	to produce a tab character
c
c	   %n	to produce a newline character. Note that the newline
c		character is NOT added by default to the end of an ouput
c		string, as it is in regular fortran. It must be explicitly
c		asked for via the %n format (as in C).
c
c	   %f	to format a real value in the form '12.236'. If the value
c		is outside the range (1.e-6, 1.e+7) then the %e format is
c		automatically used (except if the value is zero). The %f
c		format can be optionally written in one of the following
c		forms;
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
c		format can be optionally written in the following forms;
c
c		%wi	specifies that a field width of at least 'w' is used.
c			For example, %10i implies that the number 12 will be
c			printed as '        12'.
c
c
c  REVISION HISTORY:
c  1.1	turned fw.0 into iw
c  1.2	added %fw specification, where w is interpreted as the desired max.
c	number of significant digits to show on real output (Oct. 9, 1995)
c  1.21	added comment about trailing blanks in fmt not being printed (10/14/96)
c  1.3	revised %i, %f, and %e output formatting (Jun 9/97)
c  1.4	replaced \\ with char(92) for portability (Dec 1/99)
c  1.5	replaced '\n' with char(10) for portability (Sep 24/01)
c  1.6	replaced all \ escapes with % escapes (for compatibility with
c	Win98's Visual Fortran). '\' is now treated as an ordinary character
c	(Oct 1/01)
c  1.61	now properly handle '% ' at end of format string (in fact '%' at
c	end of format string will always produce a space now). (Jun 25/03)
c-------------------------------------------------------------------------
      subroutine printv( k, fmt, val, n )
      parameter (MPX = 256)
      character*(*) fmt
      character*1 bslash, spc, tab
      real*4 val(*), tmp(MPX)
      integer*4 itmp(MPX)
      logical ldot, ldigit
      equivalence (itmp(1),tmp(1))
      external lnblnk

   1  format(a,$)
   3  format()

      tab    = char(9)
      spc    = char(32)
      bslash = char(92)

      MP     = n
      if( MP .gt. MPX ) MP = MPX
      do 5 i = 1, MP
         tmp(i) = val(i)
   5  continue

      jt = 1
      lf = lnblnk(fmt)
      i  = 0
  10  i = i + 1
      if( i .gt. lf ) return
      if( fmt(i:i) .eq. '%' ) then				! DESCRIPTOR
         i = i + 1
         if( i .gt. lf ) then
            write(k,1) spc
            return
         endif
         if( fmt(i:i) .eq. 'b' ) then
            write(k,1) bslash
            go to 10
         elseif( fmt(i:i) .eq. 'n' ) then
            write(k,3)
            go to 10
         elseif( fmt(i:i) .eq. 't' ) then
            write(k,1) tab
            go to 10
         elseif( fmt(i:i) .eq. '%' ) then
            write(k,1) '%'
            go to 10
         elseif( fmt(i:i) .eq. spc ) then
            write(k,1) spc
            go to 10
         elseif( ldigit(fmt(i:i)) ) then
            call getfsp(fmt,lf,i,j,iw,id,ldot)
            if( fmt(j:j) .eq. 'f' ) then
               call prfmtf(tmp(jt),iw,id,ldot,k)
            elseif( fmt(j:j) .eq. 'e' ) then
               call prfmte(tmp(jt),iw,id,k)
            elseif( fmt(j:j) .eq. 'i' ) then
               call prfmti(tmp(jt),iw,k)
            else
               write(k,1) fmt(i-1:j)
               i = j
               go to 10
            endif
            i = j
         elseif( fmt(i:i) .eq. 'f' ) then
            call prfmtf(tmp(jt),-1,-1,.false.,k)
         elseif( fmt(i:i) .eq. 'e' ) then
            call prfmte(tmp(jt),-1,-1,k)
         elseif( fmt(i:i) .eq. 'i' ) then
            call prfmti(tmp(jt),-1,k)
         else
            write(k,1) fmt(i-1:i)
            go to 10
         endif
         jt = jt + 1
         if( jt .gt. MP ) jt = MP
         go to 10
      else							! CHARACTER
         write(k,1) fmt(i:i)
      endif
      go to 10

      end
