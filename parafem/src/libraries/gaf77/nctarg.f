c  *********************************************************************
c  *                                                                   *
c  *                   integer function nctarg                         *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.1
c  Written by Gordon A. Fenton, TUNS, Mon Sep 26 20:45:03 1994
c  Latest Update: Sep 24, 2001
c
c  PURPOSE  returns number of white-space or comma separated arguments found
c           in a character string
c
c  This function takes one argument, `str', which contains the character
c  string to parse. It traverses the character string and counts the number
c  of non-blank arguments (each separated by white-space or a comma).
c
c  REVISION HISTORY:
c  1.1	replaced the space (' ') and tab character ('\t') with char(32) and
c	char(9) for portability (Sep 24/01)
c-------------------------------------------------------------------------
      integer function nctarg(str)
      character*(*) str
      character*1 tab, spc

      nctarg = 0
      l      = len(str)
      i      = 0
      spc    = char(32)
      tab    = char(9)

  10  i = i + 1
      if( i .gt. l ) return
      if( str(i:i) .eq. spc .or. str(i:i) .eq. tab
     >                      .or. str(i:i) .eq. ',' ) go to 10
      nctarg = nctarg + 1
  20  i = i + 1
      if( i .gt. l ) return
      if( str(i:i) .eq. spc .or. str(i:i) .eq. tab
     >                      .or. str(i:i) .eq. ',' ) go to 10
      go to 20

      end
