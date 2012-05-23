c  *******************************************************************
c  *                                                                 *
c  *                    integer function nrdarg                      *
c  *                                                                 *
c  *******************************************************************
c  Character Version 1.1
c  Written by Gordon A. Fenton, TUNS, July 1994
c  Latest Update: Jun 20, 2000
c
c  PURPOSE   read a sequence of whitespace separated character strings from
c            an internal character string
c
c  This function reads a sequence of character `arguments', separated by
c  whitespace, from an internal character string and returns the actual
c  number of arguments read. For example, the string
c
c      str = 'hello, how are you?'
c
c  would be parsed into the individual arguments
c
c      cv(1) = 'hello,  '
c      cv(2) = 'how     '
c      cv(3) = 'are     '
c      cv(4) = 'you?    '
c
c  assuming that each cv is a character variable of length 8. Note that
c  argument strings are padded with blanks on the right.
c
c  Arguments are as follows;
c
c     str	character string containing the data to be read. (input)
c
c      cv	character vector of length at least MX which on output will
c		contain the sequence of character arguments read. If the
c		character argument has more characters than can be stored
c		in the cv element, the excess characters (at the end) are
c		discarded. (output)
c
c      mx	the maximum number of character arguments to read. If there
c		are more arguments in str, then nrdarg returns the actual
c		number found, but only stores the first MX in cv. (input)
c
c  REVISION HISTORY:
c  1.2	use char(9) and char(32) for tab and space for portability (Jun 20/00)
c--------------------------------------------------------------------------
      integer function nrdarg(str,cv,mx)
      character*(*) str, cv(*)
      character*1 spc, tab

      tab = char(9)
      spc = char(32)

      ls     = lnblnk(str)
      lenc   = len(cv(1))
      nrdarg = 0
      i      = 1
      j      = 1
      k      = 1
  10  if( i .gt. ls ) return
      if( (str(i:i) .eq. spc) .or. (str(i:i) .eq. tab) ) then
         i = i + 1
         go to 10
      endif
  20  if( (j .le. mx) .and. (k .le. lenc) ) cv(j)(k:k) = str(i:i)
      i = i + 1
      k = k + 1
      if( (i.gt.ls) .or. (str(i:i).eq.spc) .or. (str(i:i).eq.tab) ) then
         nrdarg = j
         if( k .le. lenc ) cv(j)(k:) = spc
         j = j + 1
         i = i + 1
         k = 1
         go to 10
      endif
      go to 20

      end
