c  *********************************************************************
c  *                                                                   *
c  *                         subroutine prins0                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Sat Sep  9 13:21:15 1995
c  Latest Update: Oct 4, 2001
c
c  PURPOSE  prints a string according to a C style format to an internal
c           character string.
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, with no numeric arguments.
c  This is primarily provided as a means to place special characters
c  such as a %, \, tab, or newline character into another string.
c  For example, to produce the string
c
c	G(\w)
c
c  in the character string `str', this routine would be called as
c
c	call print0( str, 'G(%bw)' )
c
c  Arguments to this routine are as follows;
c
c       s	character string to which output is directed. (input)
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
c  REVISION HISTORY
c-------------------------------------------------------------------------
      subroutine prins0( s, fmt )
      character*(*) s, fmt
      real tmp

      call prints( s, fmt, tmp, 0 )

      return
      end
