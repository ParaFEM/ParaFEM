c  *********************************************************************
c  *                                                                   *
c  *                         subroutine scans5                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Oct 5, 1995
c  Latest Update: Jul 13, 2004
c
c  PURPOSE  reads five real/integer values according to a C style `sscanf'
c           format from an internal character string.
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, along with five additional
c  arguments and reads their values from the given string.
c
c  In general, the format string should uniquely define the characters both
c  preceding and following the numeric value to be read. The actual number
c  read is assumed bounded by the character immediately preceding the '%'
c  descriptor in the fmt and the character immediately following the 'i'
c  or 'f' portion of the descriptor in the fmt. The following are example
c  problems encountered with scans2 and similar statements can be made about
c  this routine. The number of values actually read depends on the number
c  of descriptors encountered in the fmt string; if more descriptors
c  appear in the fmt string than arguments to this routine (5), then
c  the excess are ignored.
c
c  The function returns the number of values actually read from the string
c  (which should be 5).
c
c  EXAMPLE (2 argument case)
c  Consider trying to read the two values from the string
c
c	str = 'd = (25 x -12.2)'
c
c  One possible approach would be to invoke scans2 as follows;
c
c	n = nscan2( str, '(%i x %f)', ival, val )
c
c  In this example, the format substring '(%i ' is used to locate the integer
c  as follows;
c	1) the first occurrence of the string '(' is located in str,
c	2) the end of the number is located by finding the next ' ' in str,
c	3) the number type, i.e. integer or real, is defined by the occurrence
c	   of %i or %f, respectively,
c	4) the number is read from str
c
c  Note that an error would (probably) result if the invocation were written
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
c  minimally (note the leading space). That is, the format string should
c  uniquely define the characters both preceding and following the numeric
c  value to be read.
c
c  ARGUMENTS
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
c    val3	the third output argument. (output)
c
c    val4	the fourth output argument. (output)
c
c    val5	the fifth output argument. (output)
c
c  REVISION HISTORY:
c  1.1	converted to a function returning number of items read, changed
c	name from scans5 to nscan5. (Jul 13/04)
c-------------------------------------------------------------------------
      integer function nscan5( s, fmt, val1, val2, val3, val4, val5 )
      parameter (MP = 5)
      character*(*) fmt, s
      real val1, val2, val3, val4, val5, tmp(MP)

      nscan5 = nscanv( s, fmt, tmp, MP )

      if( nscan5 .gt. 0 ) then
         val1 = tmp(1)
         if( nscan5 .gt. 1 ) then
            val2 = tmp(2)
            if( nscan5 .gt. 2 ) then
               val3 = tmp(3)
               if( nscan5 .gt. 3 ) then
                  val4 = tmp(4)
                  if( nscan5 .gt. 4 ) then
                     val5 = tmp(5)
                  endif
               endif
            endif
         endif
      endif

      return
      end
