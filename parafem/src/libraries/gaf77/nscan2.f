c  *********************************************************************
c  *                                                                   *
c  *                         subroutine scans2                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Oct 5, 1995
c  Latest Update: Jul 13, 2004
c
c  PURPOSE  reads two real/integer values according to a C style `sscanf'
c           format from an internal character string.
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, along with two additional
c  arguments and reads their values from the given string.  The actual number
c  read is assumed bounded by the character immediately preceding the '%'
c  descriptor in the fmt and the character immediately following the 'i'
c  or 'f' portion of the descriptor in the fmt. For example, if
c  the integer '25' and real number '-12.2' are to be read from the string
c
c	str = 'd = (25 x -12.2)'
c
c  then this routine might be called as follows;
c
c	n = nscan2(str,'(%i x %f)',ival,val)
c
c  The function returns the number of values actually read from the string
c  (which should be 2).
c
c  In this example, the format substring '(%i ' is used to locate the integer
c  as follows;
c	1) the first occurrence of the string '(' is located in str,
c	2) the end of the number is located by finding the next ' ' in str,
c	3) the number type, i.e. integer or real, is defined by the occurrence
c	   of %i or %f, respectively,
c	4) the number is read from str
c
c  Note that an error would (probably) result if the invocation were
c
c	n = nscan2(str,'(%i x %f',ival,val)
c
c  since then the real would be read from the substring '-12.2)'. Note
c  also that the invocation
c
c	n = nscan2(str,'(%i x %f)',ival,val)
c
c  would not work if str = 'd(j) = (25 x -12.2)', since in this case the
c  integer would be read from the substring 'j) ' (ie this is the first
c  occurrence of the '(* ' string, where * is anything). In this case,
c  the call should be something like
c
c	n = nscan2(str,' (%i x %f)',ival,val)
c
c  minimally. That is, the format string should uniquely define the characters
c  both preceding and following the numeric value to be read.
c
c  Arguments to this routine are as follows;
c
c       s	character string from which the real/integer is read. (input)
c
c     fmt	format string. (input)
c
c    val1	the first output argument. val1 can be either real or integer,
c		it's type is determined by the presence of a corresponding
c	 	%f (for real) or %i (for integer) descriptor in fmt.
c		(output)
c
c    val2	the second output argument. (output)
c
c  REVISION HISTORY:
c  1.1	converted to a function returning number of items read, changed
c	name from scans2 to nscan2. (Jul 13/04)
c-------------------------------------------------------------------------
      integer function nscan2( s, fmt, val1, val2 )
      parameter (MP = 2)
      character*(*) fmt, s
      real val1, val2, tmp(MP)

      nscan2 = nscanv( s, fmt, tmp, MP )

      if( nscan2 .gt. 0 ) then
         val1 = tmp(1)
         if( nscan2 .gt. 1 ) then
            val2 = tmp(2)
         endif
      endif

      return
      end
