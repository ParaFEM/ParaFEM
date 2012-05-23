c  *********************************************************************
c  *                                                                   *
c  *                   integer function nscanv                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.3
c  Written by Gordon A. Fenton, TUNS, Oct 5, 1995
c  Latest Update: Aug 20, 2004
c
c  PURPOSE  reads real/integer values according to a C style `sscanf' format
c           from an internal character string.
c
c  This routine takes as input a format, which is a character string
c  containing C style format descriptors, along with a vector of additional
c  arguments and reads their values from the given string. Normally, this
c  routine is called by one of scans1, scans2, ... See their headers for usage
c  instructions.
c  Arguments to this routine are as follows;
c
c       s	character string from which data is read. (output)
c
c     fmt	format string. (input)
c
c     val	real vector of length at least MP to contain the elements
c		read from s. An element of val can be either real or integer,
c		it's type is determined by the presence of a corresponding
c	 	%f (for real) or %i (for integer) descriptor in fmt.
c		(output)
c
c       n	the number of elements in val. Currently this must be less
c		than 256. (input)
c
c  NOTES:
c  1) '%' characters are not treated specially in fmt unless they are
c     immediately followed by either 'i' or 'f'.
c  2) there is no way to actually include the strings '%i' or '%f' as
c     part of the fmt without them being treated as descriptors.
c
c  REVISION HISTORY:
c  1.1	eliminated unused vector `iderr' (Dec 5/96)
c  1.2	converted to a function returning number of items read, changed
c	name from scansv to nscanv. (Jul 13/04)
c  1.3	corrected starting counter js to start from 0. (Aug 20/04)
c-------------------------------------------------------------------------
      integer function nscanv( s, fmt, val, n )
      parameter (MPX = 256)
      character*(*) fmt, s
      real*4 val(*), tmp(MPX)
      integer*4 itmp(MPX)
      character c*2, errstr(MPX)*16
      equivalence (itmp(1),tmp(1))

   1  format(3a)
c					the number of elements to obtain
      MP = n
      if( MP .gt. MPX ) MP = MPX
      ierr = 0
c					the non-blank string lengths
      ls = lnblnk(s)
      lf = lnblnk(fmt)
c					pointers into s and fmt
      js = 0
      jf = 0
      jt = 0

  10  jf = jf + 1
      if( jf .gt. lf ) go to 40
c					check to see if we have a descriptor
      c = fmt(jf:jf+1)
      if( c .eq. '%i' .or. c .eq. '%f' ) then
         js = js + 1
         jf = jf + 1
         jt = jt + 1
         je = js
         if( jf .lt. lf ) then
  20        if( s(je:je) .ne. fmt(jf+1:jf+1) ) then
               je = je + 1
               go to 20
            endif
         else
            je = ls + 1
         endif
         if( fmt(jf:jf) .eq. 'i' ) then
            if( nrdint(s(js:je-1),itmp(jt),1) .le. 0 ) then
               ierr = ierr + 1
               errstr(ierr) = s(js:je-1)
               itmp(jt) = 0
            endif
         else
            if( nrdfp(s(js:je-1),tmp(jt),1) .le. 0 ) then
               ierr = ierr + 1
               errstr(ierr) = s(js:je-1)
               tmp(jt) = 0.0
            endif
         endif
         jf = jf + 1
         js = je
         if( (jf.gt.lf) .or. (js.gt.ls) .or. (jt.eq.MP) ) go to 40
      endif
c					find point where fmt(jf) = s(js)
  30  if( s(js:js) .ne. fmt(jf:jf) ) then
         js = js + 1
         if( js .gt. ls ) go to 40
         go to 30
      endif
      go to 10

c 40  if( ierr .gt. 0 ) then
c        do 50 i = 1, ierr
c           write(6,1)'*** scans: error reading value',
c    >                ' from following substring: ',errstr(i)
c 50     continue
c        write(6,1)'Original string: ',s(1:ls)
c        write(6,1)'Reading format : ',fmt(1:lf)
c     endif

  40  nscanv = jt

      do 60 i = 1, jt
         val(i) = tmp(i)
  60  continue

      return
      end
