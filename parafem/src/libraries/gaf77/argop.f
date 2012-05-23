c  *********************************************************************
c  *                                                                   *
c  *                       Subroutine argop                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.11
c  Written by Gordon A. Fenton, Princeton, Nov. 5, 1988.
c  Latest Update: Sep 8, 1999
c
c  PURPOSE  to obtain filenames and options from the command line
c
c  Obtains the optional argument which is assumed to be the name of the
c  input file.  This routine can be called a number of times and it will
c  get succeeding arguments and assign them to filein. Arguments are also
c  checked to see if they are proceeded by a `-' in which case they are
c  assumed to be options and the characters immediately following the `-'
c  are stored in the character vector `option'. `nopt' is the number of
c  options found. Note that the options are assumed to precede the file
c  names and that if any options appear within the list of file names,
c  they only are reported at the time they are encountered (this can be
c  used to change options before running certain files). If no further
c  arguments are found, the logical flag `lstop' is set to true.
c
c  NOTE: Options are assumed to be separate arguments (ie, don't combine
c  options after a single - character!).
c
c  Arguments to the routine are as follows;
c
c   filein    character string which on output will contain the input
c             file name as obtained from the command line arguments. (output)
c
c     lstop   logical flag which is set to true if no further filenames
c             are specified (by the user or on the command line). (output)
c
c    option   vector of character strings to contain the option strings
c             discovered in the command line arguments. (output)
c
c      nopt   number of options found. This number is cumulative from
c             call to call if nopt is not changed in the calling routine.
c             (output)
c
c  REVISION HISTORY:
c  1.1	now getarg is an intrinsic (Jul 27/99)
c  1.11	emphasized that options must appear separately above (Sep 8/99)
c-------------------------------------------------------------------------
      subroutine argop( filein, lstop, option, nopt )
      logical lstop
      integer nopt
      character*(*) filein, option(*)
      save narg
      data karg/0/

   1  format(a)
   2  format(a,$)

      lstop = .false.

      if( karg .eq. 0 ) then
         nopt = 0
         narg = iargc()
      endif

      karg = karg + 1
      if( narg .eq. 0 ) then
          if( karg .eq. 1 ) then
              write(6,2)'Input filename? '
              read(5,1) filein
          else
              go to 30
          endif
      else
  10      if( karg .gt. narg ) go to 30
          call getarg( karg, filein )
          if( filein(1:1) .eq. '-' ) then
              karg = karg + 1
              nopt = nopt + 1
              option(nopt) = filein(2:)
              go to 10
          endif
      endif

      return

  30  lstop = .true.
      return

      end
