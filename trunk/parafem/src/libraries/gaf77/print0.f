c  *********************************************************************
c  *                                                                   *
c  *                         subroutine print0                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Sun Sep 8, 2002
c  Latest Update: Sep 8, 2002
c
c  PURPOSE  prints a string according to a C style format
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, with no numeric arguments.
c  This is primarily provided as a means to print special characters
c  such as a %, \, tab, or newline character. For example, to print
c
c	G(\w)
c
c  this routine would be called as
c
c	call print0( 6, 'G(%bw)%n' )
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
c	   %n	to produce a newline character. This is NOT added to the
c		end of the output string by default. It must be explicitly
c		added (as in C).
c
c  NOTES:
c   1)	trailing white space in the format string is not printed in the output.
c	To actually get it in the output, it must be escaped, as in % % % ...
c
c  REVISION HISTORY
c-------------------------------------------------------------------------
      subroutine print0( k, fmt )
      character*(*) fmt
      real val

      call printv( k, fmt, val, 0 )

      return
      end
