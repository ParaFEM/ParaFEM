c  *********************************************************************
c  *                                                                   *
c  *                    integer function nscan1                        *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Oct 5, 1995
c  Latest Update: Jul 13, 2004
c
c  PURPOSE  read one real/integer value according to a C style `sscanf' format
c           from an internal character string.
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, along with one additional
c  argument and reads its value from the given string. For example, if
c  the integer '25' is to be read from the string
c
c	str = 'd = (25)'
c
c  then this routine might be called as follows;
c
c	n =  nscan1(str,'(%i)',ival)
c
c  The function returns the number of values actually read from the string
c  (which should be 1).
c
c  In this example, the format string '(%i)' is used to locate the integer as
c  follows;
c	1) the first occurrence of the string '(' is located in str,
c	2) the end of the number is located by finding the next ')' in str,
c	3) the number type, i.e. integer or real, is defined by the occurrence
c	   of %i or %f, respectively,
c	4) the number is read from str
c
c  Note that an error would (probably) result if the invocation were
c
c	n = nscan1(str,'(%i',ival)
c
c  since then the integer would be read from the substring '25)'. Note
c  also that the invocation
c
c	n = nscan1(str,'(%i)',ival)
c
c  would not work if str = 'd(j) = (25)', since in this case the integer
c  would be read from the substring 'j' (ie this is the first occurrence
c  of the (...) string). In this case, the call should be something like
c
c	n = nscan1(str,' (%i)',ival)
c
c  minimally. That is, the format string should uniquely define the
c  characters both preceding and following the numeric value to be read.
c
c  Arguments to this routine are as follows;
c
c       s	character string from which the real/integer is read. (input)
c
c     fmt	format string. (input)
c
c    val1	the output argument. val1 can be either real or integer,
c		it's type is determined by the presence of a corresponding
c	 	%f (for real) or %i (for integer) descriptor in fmt.
c		(output)
c
c  REVISION HISTORY:
c  1.1	converted to a function returning number of items read, changed
c	name from scans1 to nscan1. (Jul 13/04)
c-------------------------------------------------------------------------
      integer function nscan1( s, fmt, val1 )
      parameter (MP = 1)
      character*(*) fmt, s
      real val1, tmp(MP)

      nscan1 = nscanv( s, fmt, tmp, MP )

      if( nscan1 == 1 ) val1 = tmp(1)

      return
      end
