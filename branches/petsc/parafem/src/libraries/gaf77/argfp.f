c  *********************************************************************
c  *                                                                   *
c  *                       Subroutine argfp                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, Princeton, Feb. 19, 1988.
c  Latest Update: Jul 27, 1999
c
c  PURPOSE  to obtain input file names and options and create output names
c
c  Obtains the optional argument which is assumed to be the name of the
c  input file. argfp then strips the name of its extension and adds the
c  extension 'ext' to form the output name fileout. If no arguments are
c  supplied, this routine prompts for a filename. This routine can
c  be called a number of times and it will get succeeding arguments and
c  assign them to filein. When the argument list is exhausted, the logical
c  flag `lstop' is set to true. Arguments are also checked
c  to see if they are proceeded by a `-' in which case they are assumed
c  to be options and the characters immediately following the `-' are
c  stored in the character vector `option'. `nopt' is the number of options
c  found. Note that the options are assumed to preceed the file names and
c  that if any options appear within the list of file names, they only
c  are reported at the time they are encountered (this can be used to change
c  options before running certain files). Arguments to the routine are
c  as follows;
c
c   filein    character string which on output will contain the input
c             file name as obtained from the command line arguments. (output)
c
c   fileout   character string which on output will contain the output
c             file name created by appending `ext' to the basename of
c             `filein'. (output)
c
c       ext   character string denoting the desired extension of the
c             output file name. `ext' should not include the leading `.'
c             (input)
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
c-------------------------------------------------------------------------
      subroutine argfp( filein, fileout, ext, lstop, option, nopt )
      logical lstop
      integer nopt
      character*(*) filein, fileout, ext, option(1)
      save narg
      data karg/0/

   1  format(a)
   2  format(a,$)

      lstop = .false.

      if( karg .eq. 0 ) then
          narg = iargc()
          nopt = 0
      endif

      karg = karg + 1
      if( narg .eq. 0 ) then
          if( karg .eq. 1 ) then
              write(6,2)'Input filename? '
              read(5,1) filein
          else
              go to 20
          endif
      else
   5      if( karg .gt. narg ) go to 20
          call getarg( karg, filein )
          if( filein(1:1) .eq. '-' ) then
              karg = karg + 1
              nopt = nopt + 1
              option(nopt) = filein(2:)
              go to 5
          endif
      endif

c                                 construct output filename
  10  i = lnblnk( filein )
      do 12 il = i, 1, -1
         if( filein(il:il) .eq. '.' ) go to 15
  12  continue
      il = i + 1
  15  il = il - 1
      fileout = filein(1:il)//'.'//ext

      return

  20  lstop = .true.
      return

      end
