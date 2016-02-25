c  *********************************************************************
c  *                                                                   *
c  *                         subroutine prins6                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Sat Sep  9 13:21:15 1995
c  Latest Update: Oct 4, 2001
c
c  PURPOSE  prints six values according to a C style format to an internal
c           character string.
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, along with six additional
c  arguments and prints them to the given string. For example, say you wish
c  to print the values of the variables 'i', `d', 'j',  'f', 'nsim', and 'ierr'
c  to the string sub so that the string appears as
c
c	d(3) = 16.49, f(12) = -0.03 (6'th realization, ierr = 0)
c
c  (assuming i = 3, d = 16.49, j = 12, f = -0.03, nsim = 6, and ierr = 0),
c  then this routine would be called as follows;
c
c	call prins6( sub,
c      >    'd(%i) = %f, f(%i) = %f (%i''th realization, ierr = %i)',
c      >             i, d, j, f, nsim, ierr )
c
c  If the resulting string exceeds the dimensioned length of the target
c  string, then an error message is issued (however, this routine may still
c  have written past the end of the string variable).
c
c  Arguments to this routine are as follows;
c
c       s	character string to which output is directed. (input)
c
c     fmt	format string. (input)
c
c    val1	the first input argument. val1 can be either real or integer,
c		it's type is determined by the presence of a corresponding
c	 	%f or %e (for real) or %i (for integer) descriptor. See notes
c		below on the optional format specification which may precede
c		any of these format descriptors.
c		(input)
c
c    val2	see val1. (input)
c
c    val3	see val1. (input)
c
c    val4	see val1. (input)
c
c    val5	see val1. (input)
c
c    val6	see val1. (input)
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
c  REVISION HISTORY
c  1.1	added %fw specification, where w is interpreted as the desired max.
c	number of significant digits to show on real output (Oct 9/95)
c  1.11	modified above writeup for new prints routines (Jun 9/97)
c  1.2	modified above writeup to reflect new formats and escapes (Oct 4/01)
c-------------------------------------------------------------------------
      subroutine prins6( s, fmt, val1, val2, val3, val4, val5, val6 )
      parameter (MP = 6)
      character*(*) fmt, s
      real val1, val2, val3, val4, val5, val6, tmp(MP)

      tmp(1) = val1
      tmp(2) = val2
      tmp(3) = val3
      tmp(4) = val4
      tmp(5) = val5
      tmp(6) = val6

      call prints( s, fmt, tmp, MP )

      return
      end
