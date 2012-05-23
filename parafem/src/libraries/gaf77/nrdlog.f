c  *******************************************************************
c  *                                                                 *
c  *                    integer function nrdlog                      *
c  *                                                                 *
c  *******************************************************************
c  Integer Version 2.2
c  Written by Gordon A. Fenton, TUNS, Dec. 2, 1993
c  Latest Update: Jun 20, 2000
c
c  PURPOSE   read a sequence of logical variables from an internal character
c            string using free format
c
c  This function performs a free format internal string read of logical
c  similar to what would be obtained by `read(str,*) log1, log2, ...' and
c  returns the actual number of logicals read. The logical values are
c  determined purely by their first letters, `f' or `t', ignoring leading
c  decimals. Thus .false., false, and fool are all considered to be false.
c  If an error is encountered during the reading of the first MX or less
c  logicals, the function returns the negative of the index of the logical
c  at which the error occurred. For example, trying to read the string
c
c       str = 't f .true. hello false'
c
c  with the reference `n = nrdlog( str, lv, 20 )' would result in n = -4
c  with lv(1) = .true., lv(2) = .false., lv(3) = .true.
c  Arguments are as follows;
c
c    str    character string containing the data to be read. (input)
c
c    lv     logical vector of length at least MX which on output will
c	    contain the sequence of logical variables read. (output)
c
c    MX     the maximum number of logical values to read. If there
c	    are more values in the string, then nrdlog returns the
c           actual number found, but only stores the first MX in nv.
c	    (input)
c
c  NOTE: this function replaces rdlogs.f
c
c  REVISION HISTORY:
c  2.1	changed `i .le. ls' to `i .lt. ls' in block 20 below (July 15, 1994)
c  2.2	use char(9) and char(32) for tab and space for portability (Jun 20/00)
c--------------------------------------------------------------------------
      integer function nrdlog(str,lv,MX)
      logical lv(*), lf, lc
      character*(*) str
      character*1 c, tab, spc

      tab = char(9)
      spc = char(32)

      ls = lnblnk(str)
      nrdlog = 0
      i = 1
      j = 1
  10  if( i .gt. ls ) return
      c = str(i:i)
      if( (c .eq. spc) .or. (c .eq. tab)
     >.or.(c .eq. ',') .or. (c .eq. '.') ) then
         i = i + 1
         go to 10
      endif
      lf = ((c .eq. 't') .or. (c .eq. 'T'))
      if( lf .or. c .eq. 'f' .or. c .eq. 'F' ) then
         if( j .le. MX ) lv(j) = lf
         nrdlog = j
         j = j + 1
  20     i = i + 1
         c = str(i:i)
         lc =  ( (c .ne. spc) .and. (c .ne. tab) .and. (c .ne. ',') )
         if( lc .and. i .lt. ls ) go to 20
      else
         if( j .gt. MX ) return
         nrdlog = -j
         return
      endif
      go to 10

      end
