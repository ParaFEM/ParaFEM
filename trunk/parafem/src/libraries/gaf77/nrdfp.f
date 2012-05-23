c  *******************************************************************
c  *                                                                 *
c  *                    integer function nrdfp                       *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 2.4
c  Written by Gordon A. Fenton, TUNS, Dec. 2, 1993
c  Latest Update: Nov 26, 2004
c
c  PURPOSE   read a sequence of real variables from an internal character
c            string using free format with provision for sequences
c
c  This function performs a free format internal string read of real values
c  similar to what would be obtained by `read(str,*) fp1, fp2, ...' and
c  returns the actual number of real values read. Values must be separated
c  by whitespace (not commas). If an error is encountered during the reading
c  of the first MX or less real values, the function returns the negative
c  of the index of the white-space separated element at which the error
c  occurred. For example, trying to read the string
c
c       str = '2.2e+03  3.2  45  5.r2  6'
c
c  with the reference `n = nrdfp( str, v, 20 )' would result in n = -4
c  with v(1) = 2200.0, v(2) = 3.2, v(3) = 45.0. Note that the string (see
c  discussion below)
c
c       str = '2.2e+03 3.2-4.2,0.2 45 5.r2 6'
c
c  would also result in n = -4 since the `5.r2' element is the 4th (it is
c  the 9th element in the expanded sequence, as discussed next).
c
c  This function also can handle sequences of the form f1-f2,f3 (with no
c  intervening white space). This is interpreted as the sequence
c
c	f1, f1+f3, f1+2*f3, ..., f1+n*f3
c
c  where n = int[(f2 - f1)/f3]. The value of `n' must be greater than or
c  equal to 0. The function returns the number of elements in the constructed
c  sequence. If the `,f3' entry is absent, f3 defaults to 1.0. The string
c  `f1,f3' is illegal.
c  For example, the string
c
c	1-10,2	12-15 -1--5,-2 0.1-0.3,0.1
c
c  expands into the sequence of values;
c
c	1. 3. 5. 7. 9.  12. 13. 14. 15.  -1. -3. -5.	0.1 0.2 0.3
c
c  and nrdfp returns 15.
c
c  Arguments are as follows;
c
c    str    character string containing the data to be read. (input)
c
c    v      real vector of length at least MX which on output will
c	    contain the sequence of real variables read. (output)
c
c    MX     the maximum number of real values to read. If there
c	    are more values in the string, then nrdfp returns the
c           actual number found, but only stores the first MX in v.
c	    (input)
c
c  REVISION HISTORY:
c  2.1	changed `i .gt. ls' to `i .ge. ls' in block 20 below.
c  2.2	provided for sequences of the form f1-f2 or f1-f2,f3. No longer
c	allows for commas to separate numeric elements except in sequences.
c	(June 21, 1995)
c  2.3	use char(9) and char(32) for tab and space for portability (Jun 20/00)
c  2.4	make f(3) real*4 so we don't have to worry about specifying too
c	many significant digits and overrun an integer f(3). (Nov 26/04)
c--------------------------------------------------------------------------
      integer function nrdfp(str,v,MX)
      real v(*)
      real f(3)
      character*(*) str
      character*1 c, tab, spc
      logical lseq, linc

      nrdfp = 0
      tab = char(9)
      spc = char(32)

      ls   = lnblnk(str)
      i    = 1
      j    = 1
      k    = 1
      s    = 1.
      me   = 1
      f(1) = 0.
      f(2) = 0.
      f(3) = 0.
      f3   = 1.
      lseq = .false.
      linc = .false.
      i2   = 0
      i3   = 0
      jt   = 1

  10  if( i .gt. ls ) return
      c = str(i:i)
      if( (c .eq. spc) .or. (c .eq. tab) ) then		! skip leading space
         i = i + 1
         go to 10
      endif
      if( c .eq. '-' ) then				! set global sign
         s = -1.
         i = i + 1
         c = str(i:i)
      endif

  20  if( c .eq. '0' ) then
         f(k) = 10.*f(k)
      elseif( c .eq. '1' ) then
         f(k) = 10.*f(k) + 1.
      elseif( c .eq. '2' ) then
         f(k) = 10.*f(k) + 2.
      elseif( c .eq. '3' ) then
         f(k) = 10.*f(k) + 3.
      elseif( c .eq. '4' ) then
         f(k) = 10.*f(k) + 4.
      elseif( c .eq. '5' ) then
         f(k) = 10.*f(k) + 5.
      elseif( c .eq. '6' ) then
         f(k) = 10.*f(k) + 6.
      elseif( c .eq. '7' ) then
         f(k) = 10.*f(k) + 7.
      elseif( c .eq. '8' ) then
         f(k) = 10.*f(k) + 8.
      elseif( c .eq. '9' ) then
         f(k) = 10.*f(k) + 9.
      elseif( c .eq. 'e' .or. c .eq. 'E'
     >   .or. c .eq. 'd' .or. c .eq. 'D' ) then
         if( k .eq. 3 ) go to 999
         if( k .eq. 2 ) i3 = i
         k = 3
         i = i + 1
         c = str(i:i)
         if( c .eq. '+' .or. c .eq. '-' ) then
            if( c .eq. '-' ) me = -1
            i = i + 1
            c = str(i:i)
         endif
         go to 20
      elseif( c .eq. '.' ) then
         if( k .eq. 2 .or. k .eq. 3 ) go to 999
         k  = 2
         i2 = i+1
      elseif( c .eq. '-' ) then
         lseq = .true.
         if( f(3) .ne. 0. ) then
            s = s*10.**(float(me)*f(3))
         endif
         if( k .eq. 2 ) i3 = i
         fr   = f(2)/(10.**(i3-i2))
         f1   = s*(f(1) + fr)
         i    = i + 1
         s    = 1.
         me   = 1
         f(1) = 0.
         f(2) = 0.
         f(3) = 0.
         k    = 1
         i2   = 0
         i3   = 0
         go to 10
      elseif( c .eq. ',' ) then
         if( .not. lseq ) go to 999
         linc = .true.
         if( f(3) .ne. 0. ) then
            s = s*10.**(float(me)*f(3))
         endif
         if( k .eq. 2 ) i3 = i
         fr   = f(2)/(10.**(i3-i2))
         f2   = s*(f(1) + fr)
         i    = i + 1
         s    = 1.
         me   = 1
         f(1) = 0.
         f(2) = 0.
         f(3) = 0.
         k    = 1
         i2   = 0
         i3   = 0
         go to 10
      elseif( (c .eq. spc) .or. (c .eq. tab) .or. (i .ge. ls) ) then
         if( lseq ) then
            if( f(3) .ne. 0. ) then
               s = s*10.**(float(me)*f(3))
            endif
            if( k .eq. 2 ) i3 = i
            fr = f(2)/(10.**(i3-i2))
            if( linc ) then
               f3 = s*(f(1) + fr)
            else
               f2 = s*(f(1) + fr)
            endif
            if( f3 .eq. 0. ) go to 999
            n  = int( 1.000001*(f2 - f1)/f3 )
            if( n .lt. 0 ) go to 999
            do 30 L = 0, n
               if( j .le. MX ) v(j) = f1 + float(L)*f3
               j = j + 1
  30        continue
            nrdfp = j - 1
            lseq  = .false.
            linc  = .false.
            f3    = 1.0
         else
            if( j .le. MX ) then
               if( f(3) .ne. 0 ) then
                  s = s*10.**(float(me)*f(3))
               endif
               if( k .eq. 2 ) i3 = i
               fr = f(2)/(10.**(i3-i2))
               v(j) = s*(f(1) + fr)
            endif
            nrdfp = j
            j     = j + 1
         endif
         jt   = jt + 1
         i    = i + 1
         s    = 1.
         me   = 1
         f(1) = 0.
         f(2) = 0.
         f(3) = 0.
         k    = 1
         i2   = 0
         i3   = 0
         go to 10
      else
         go to 999
      endif
      i = i + 1
      if( i .gt. ls ) then
         c = spc
      else
         c = str(i:i)
      endif
      go to 20
c						read error
 999  if( j .gt. MX ) return
      nrdfp = -jt
      return
      end
