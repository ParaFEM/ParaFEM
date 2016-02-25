c  *********************************************************************
c  *                                                                   *
c  *                          subroutine getop                         *
c  *                                                                   *
c  *********************************************************************
c  Character Version 1.1
c  Written by R. Dan Lohnes, TUNS, Sep 1994.
c  Latest Update: Jul 27, 1999
c
c  PURPOSE   to return the next option letter and argument from the command
c            line.
c
c  This is modeled after "getopt(3)" except that a semicolon (;) 
c  following an option letter means that a following argument is
c  asscociated which occurs immediately, with no intervening white
c  space.  In this case, the argument is considered to be optional.  
c  A colon (:) following an option character denotes that a following 
c  argument is associated and that it follows after a white space.  
c  In this case, the following argument is mandatory.  In both cases, 
c  the option must appear on its own or at the end of a string of other 
c  options.
c
c  This routine returns one option at a time; eoa is set to true when the
c  argument list has been used up.
c
c  If the next argument in the command line is not an option, then c is
c  set to ' ' and optarg is set to the value of the argument.
c
c  Concatenated options are allowed, although options having arguments
c  must appear at the end of the option string. For example if
c  optstr = 'a;bcd:' (see below for more details), then the option string
c  `-bca12.2' is equivalent to `-b -c -a12.2'. This routine will first
c  return the -b option, then the -c option, then the -a option with its
c  optional argument string `12.2' in optarg on subsequent calls.
c
c  Arguments to this routine are as follows;
c
c       c	char string of length at least 1 which returns the option
c		letter. (output)
c
c  optarg	char string which returns the argument to the option if
c		present. (output)
c
c  optstr	char string of all possible options. In opstr, the characters
c		`;' and `:' have special meaning when following an option
c		letter. The `;' implies that the option is to be immediately
c		followed by an optional character string with no intervening
c		white space. If the string is missing, optarg is set blank.
c		The `:' implies that the option is to be followed by a
c		non-optional argument and will be separated from the option
c		letter by white space. For example if optstr = 'a;bcd:' then
c		options -a, -b, -c, and -d are recognized. The -a option may
c		be immediately followed by an argument string (as in -a12.2)
c		while the -d option must be followed after some white space
c		by an argument (as in -d hello). optstr should not be changed
c		after the first call to this routine. (input)
c
c     eoa	logical flag which is set to true if all arguments supplied on
c		the command line have been considered. (output)
c
c  REVISION HISTORY:
c  1.1	now getarg is an intrinsic (Jul 27/99)
c--------------------------------------------------------------------------
      subroutine getop( c, optarg, optstr, eoa )
      character*256 str
      character*(*) c, optarg, optstr
      logical eoa
      save narg, iarg, j, str, lenopt, lenstr
      data narg/-1/, j/1/, iarg/0/

   1  format(a,a,a)
c					initialize things
      if( narg .eq. -1 ) then
         narg   = iargc()
         lenopt = lnblnk(optstr)
         eoa    = .false.
      endif
c					clear optional argument string
      optarg = ' '
c					get next argument from command line
      if( j .eq. 1 ) then
         iarg = iarg + 1
         if( iarg .gt. narg ) then
            eoa = .true.
            return
         endif
         call getarg( iarg, str )
         lenstr = lnblnk(str)
      endif
c					check if an option
      itest = 0
      if( str(1:1).eq.'-' ) then
         j = j + 1
         do 20 i = 1, lenopt
            if( str(j:j) .eq. optstr(i:i) ) then
               itest = 1
               c     = str(j:j)
c					check for `:' after option

               if( optstr(i+1:i+1).eq.':' ) then
                  iarg = iarg + 1
                  if ( iarg .gt. narg ) then
                     write(6,1)'Option -',str(j:j),
     >                         ' is missing required argument.'
                     stop
                  else 
                     call getarg( iarg, optarg )
                     if( optarg(1:1) .eq. '-' ) then
                        write(6,1)'Option -',str(j:j),
     >                            ' is missing required argument.'
                        stop
                     endif
                  endif
                  j = 1 
                  return
c					check for `;' after option

               elseif( optstr(i+1:i+1).eq.';' ) then
                  if( lenstr .ge. (j+1) ) optarg = str(j+1:lenstr)
                  j = 1
                  return
               endif
               go to 30
            endif
  20     continue
c					did we find a valid option?
  30     if( itest. eq. 0 ) then
            write(6,1)'Unknown option: -',str(j:j)
            stop
         endif
      else
c					not an option, return it in optarg
         c = ' '
         optarg = str
      endif
c					reset j
      if( j .eq. lenstr ) j = 1

      return
      end
