c  *********************************************************************
c  *                                                                   *
c  *                         subroutine print8                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Sat Sep  9 13:21:15 1995
c  Latest Update: Oct 4, 2001
c
c  PURPOSE  prints eight values according to a C style format
c
c  This routine takes as input a unit number, a format, which is a character
c  string containing C style format descriptors, and eight additional
c  arguments and prints them to the given unit number. For example, say
c  you wish to print the values of the variables nx1, nx2, nx3, nye, xl,
c  yl, w1, and w2 to the string sub so that the string appears as
c
c	'n_e = (12 + 24 + 12 x 20), D = (1.4 x 3.6), w_1 = 0.8, w_2 = 0.4'
c
c  (assuming nx1 = 12, nx2 = 24, nx3 = 12, nye = 20, xl = 1.4, yl = 0.8,
c  w1 = 0.8, and w2 = 0.4), then this routine would be called as follows;
c
c        call print8(
c        'n_e = (%i + %i + %i x %i), D = (%f x %f), w_1 = %f, w_2 = %f%n',
c        nx1,nx2,nx3,nye,xl,yl,w1,w2)
c
c  Note that if the trailing '%n' characters are absent from the format string,
c  then the carriage return is inhibited in the output.
c
c  Arguments to this routine are as follows;
c
c       k	unit number to which output is directed. (input)
c
c     fmt	format string. (input)
c
c    val1	the first input argument. val1 can be either real or integer;
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
c    val7	see val1. (input)
c
c    val8	see val1. (input)
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
c	   %n	to produce a newline character. This is NOT added to the
c		end of the output string by default. It must be explicitly
c		added (as in C).
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
c  NOTES:
c   1)	if more than 8 occurrences of the formats %f, %e or %i appear in
c	the format string, the last argument value is repeatedly used (this may
c	lead to problems if it's type doesn't match the repeated format
c	descriptor).
c
c   2)	trailing white space in the format string is not printed in the output.
c	To actually get it in the output, it must be escaped, as in % % % ...
c
c  REVISION HISTORY:
c-------------------------------------------------------------------------
      subroutine print8(k,fmt,val1,val2,val3,val4,val5,val6,val7,val8)
      parameter (MP = 8)
      character fmt*(*)
      real val1, val2, val3, val4, val5, val6, val7, val8, tmp(MP)

      tmp(1) = val1
      tmp(2) = val2
      tmp(3) = val3
      tmp(4) = val4
      tmp(5) = val5
      tmp(6) = val6
      tmp(7) = val7
      tmp(8) = val8

      call printv( k, fmt, tmp, MP )

      return
      end
