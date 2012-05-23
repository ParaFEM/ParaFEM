c  *********************************************************************
c  *                                                                   *
c  *                         subroutine getfsp                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Wed Jun  4 22:19:27 1997
c
c  PURPOSE  gets the format specification from the format string
c
c  This routine examines the string fmt(i+1:j) for a format specification
c  of the form `iw.id', where `iw' and `id' are whole numbers separated by
c  a decimal place, `.'. The index j is set to point at the next character
c  in fmt which is not either a digit in the range [0,9] or a decimal point.
c  The format specification can take one of four forms
c
c	1) both iw and id specified, separated by a `.'
c	2) only iw specified, id is set to -1
c	3) only .id specified, iw is set to -1
c	4) neither specified, both id and iw are set to -1
c
c  Arguments to this routine are as follows;
c
c     fmt	character string optionally containing the format specification
c		starting at index (i+1). (input)
c
c      lf	integer giving the non-blank length of the format string, fmt.
c		(input)
c
c       i	integer index into fmt at the point where the iw.id format
c		specification is searched for. (input)
c
c	j	integer index into fmt which on output points to the next
c		character in fmt beyond the iw.id format specification. If no
c		specification is found, j is set to i. (output)
c
c      iw	integer which on output will contain the left hand portion
c		of the format specification. If not provided, iw is set to
c		-1. (output)
c
c      id	integer which on output will contain the right hand portion
c		of the format specification. If not provided, id is set to
c		-1. (output)
c
c    ldot	logical flag which is set to true if a decimal point was
c		found in the format specification. This allows format
c		specifications of the form `iw.', without a trailing
c		id specification, for example. (output)
c-------------------------------------------------------------------------
      subroutine getfsp(fmt,lf,i,j,iw,id,ldot)
      character*(*) fmt
      character*1 c
      integer iq(2), il(2)
      logical ldot
c					process format descriptor
      ldot  = .false.			! true if decimal point found
      iq(1) = 0				! field width (later iw)
      iq(2) = 0				! no. decimal digits (later id)
      il(1) = -1			! default if iw not specified
      il(2) = -1			! default if id not specified
      m     = 1
      j     = i
  10  if( j .lt. lf ) then
         c     = fmt(j:j)
         if( c .eq. '0' ) then
            iq(m) = 10*iq(m)
            il(m) = 0
         elseif( c .eq. '1' ) then
            iq(m) = 10*iq(m) + 1
            il(m) = 0
         elseif( c .eq. '2' ) then
            iq(m) = 10*iq(m) + 2
            il(m) = 0
         elseif( c .eq. '3' ) then
            iq(m) = 10*iq(m) + 3
            il(m) = 0
         elseif( c .eq. '4' ) then
            iq(m) = 10*iq(m) + 4
            il(m) = 0
         elseif( c .eq. '5' ) then
            iq(m) = 10*iq(m) + 5
            il(m) = 0
         elseif( c .eq. '6' ) then
            iq(m) = 10*iq(m) + 6
            il(m) = 0
         elseif( c .eq. '7' ) then
            iq(m) = 10*iq(m) + 7
            il(m) = 0
         elseif( c .eq. '8' ) then
            iq(m) = 10*iq(m) + 8
            il(m) = 0
         elseif( c .eq. '9' ) then
            iq(m) = 10*iq(m) + 9
            il(m) = 0
         elseif( c .eq. '.' ) then
            if( ldot ) go to 20		! done, found specification
            ldot = .true.
            m    = 2
         else
            go to 20			! done, found specification
         endif
         j = j + 1
         go to 10
      endif
c					now assign the specification
  20  iw = iq(1) + il(1)
      id = iq(2) + il(2)

      return
      end
