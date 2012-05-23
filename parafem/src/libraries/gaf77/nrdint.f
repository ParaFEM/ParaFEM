c  *******************************************************************
c  *                                                                 *
c  *                    integer function nrdint                      *
c  *                                                                 *
c  *******************************************************************
c  Integer Version 2.3
c  Written by Gordon A. Fenton, TUNS, Dec. 2, 1993
c  Latest Update: Jun 20, 2000
c
c  PURPOSE   read a sequence of integer variables from an internal character
c            string using free format with special provision for sequences
c
c  This function performs a free format internal string read of integers
c  similar to what would be obtained by `read(str,*) int1, int2, ...' and
c  returns the actual number of integers read. Numeric entries in the string
c  must be separated by white space (not commas).
c
c  If an error is encountered during the reading of the first MX or less
c  real values, the function returns the negative of the index of the
c  white-space separated element at which the error occurred. For example,
c  trying to read the string
c
c       str = '2 3 5.2 6'
c
c  with the reference `n = nrdint( str, nv, 20 )' would result in n = -3
c  with nv(1) = 2, nv(2) = 3. Note that the string (see discussion below)
c
c       str = '2 3-7,2 5.2 6'
c
c  would also result in n = -3 since the `5.2' element is the 3rd white-
c  space separated argument (it is the 5th element in the expanded sequence,
c  as discussed next).
c
c  This function also can handle sequences of the form i1-i2,i3 (with no
c  intervening white space). This is interpreted as the sequence
c
c	i1, i1+i3, i1+2*i3, ..., i1+n*i3
c
c  where n = int[(i2 - i1)/i3]. The value of `n' must be greater than or
c  equal to 0. The function returns the number of elements in the constructed
c  sequence. If the `,i3' entry is absent, i3 defaults to 1.
c  For example, the string
c
c	1-10,2	12-15 -1--5,-2
c
c  expands into the sequence of values;
c
c	1 3 5 7 9  12 13 14 15  -1 -3 -5
c
c  and nrdint returns 12.
c
c  The sequence `i1,i3' is interpreted slightly differently: it expands into
c  i3 identical values, all equal to i1. For example, the string `2,3' expands
c  into the sequence `2 2 2'.
c
c  Arguments are as follows;
c
c    str    character string containing the data to be read. (input)
c
c    nv     integer vector of length at least MX which on output will
c	    contain the sequence of integer variables read. (output)
c
c    MX     the maximum number of integer values to read. If there
c	    are more values in the string, then nrdint returns the
c           actual number found, but only stores the first MX in nv.
c	    (input)
c
c  REVISION HISTORY:
c  2.1	no longer reads beyond the end of the string.
c  2.2	provided for sequences of the form i1-i2 or i1-i2,i3. No longer
c	allows for commas to separate numeric elements except in sequences.
c	(June 21, 1995)
c  2.3	use char(9) and char(32) for tab and space for portability (Jun 20/00)
c--------------------------------------------------------------------------
      integer function nrdint(str,nv,MX)
      integer nv(*)
      character*(*) str
      character*1 c, tab, spc
      logical lseq, linc

      nrdint = 0
      tab = char(9)
      spc = char(32)

      ls   = lnblnk(str)
      i    = 1
      j    = 1
      jt   = 1
      k    = 0
      m    = 1
      i3   = 1
      lseq = .false.
      linc = .false.

  10  if( i .gt. ls ) return
      c = str(i:i)
      if( (c .eq. spc) .or. (c .eq. tab) ) then
         i = i + 1
         go to 10
      endif
      if( c .eq. '-' ) then
         m = -1
         i = i + 1
         c = str(i:i)
      endif
  20  if( c .eq. '0' ) then
         k = 10*k
      elseif( c .eq. '1' ) then
         k = 10*k + 1
      elseif( c .eq. '2' ) then
         k = 10*k + 2
      elseif( c .eq. '3' ) then
         k = 10*k + 3
      elseif( c .eq. '4' ) then
         k = 10*k + 4
      elseif( c .eq. '5' ) then
         k = 10*k + 5
      elseif( c .eq. '6' ) then
         k = 10*k + 6
      elseif( c .eq. '7' ) then
         k = 10*k + 7
      elseif( c .eq. '8' ) then
         k = 10*k + 8
      elseif( c .eq. '9' ) then
         k = 10*k + 9
      elseif( c .eq. '-' ) then
         lseq = .true.
         i1 = m*k
         i  = i + 1
         k  = 0
         m  = 1
         go to 10
      elseif( c .eq. ',' ) then
c        if( .not. lseq ) go to 999
         linc = .true.
         i2 = m*k
         i  = i + 1
         k  = 0
         m  = 1
         go to 10
      elseif( (c .eq. spc) .or. (c .eq. tab) .or. (i .ge. ls) ) then
         if( lseq ) then
            if( linc ) then
               i3 = m*k
            else
               i2 = m*k
            endif
            if( i3 .eq. 0 ) go to 999
            n  = (i2 - i1)/i3
            if( n .lt. 0 ) go to 999
            do 30 L = 0, n
               if( j .le. MX ) nv(j) = i1 + L*i3
               j = j + 1
  30        continue
            nrdint = j - 1
            lseq = .false.
            linc = .false.
            i3   = 1
         elseif( linc ) then
            i3 = m*k
            do 40 L = 1, i3
               if( j .le. MX ) nv(j) = i2
               j = j + 1
  40        continue
            nrdint = j - 1
            linc   = .false.
            i3     = 1
         else
            if( j .le. MX ) nv(j) = m*k
            nrdint = j
            j = j + 1
         endif
         jt = jt + 1
         i  = i + 1
         k  = 0
         m  = 1
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
      nrdint = -jt
      return
      end
