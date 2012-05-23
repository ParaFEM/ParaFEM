c  *********************************************************************
c  *                                                                   *
c  *                       Subroutine argfn                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, Princeton, Feb. 18, 1988.
c  Latest Update: Jul 27, 1999
c
c  PURPOSE  to obtain input file argument and create output file name
c
c  Obtains the optional argument which is assumed to be the name of the
c  input file. argfn then strips the name of its extension and adds the
c  extension 'ext' to form the output filename. This routine can
c  be called a number of times and it will get succeeding arguments and
c  assign them to filein. When the argument list is exhausted, the user
c  is prompted for any further input files. The logical flag `lstop'
c  is set to true if the user does not specify a filename.
c  Note: `ext' should not include the `.' Arguments to this routine
c  are as follows;
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
c  REVISION HISTORY:
c  1.1	now getarg is an intrinsic (Jul 27/99)
c--------------------------------------------------------------------------
      subroutine argfn( filein, fileout, ext, lstop )
      logical lstop
      character*1 chr
      character*(*) filein, fileout, ext
      save narg
      data karg/0/

   1  format(a)
   2  format(a,$)

      lstop = .false.

      if( karg .eq. 0 ) narg = iargc()
      karg = karg + 1
      if( narg .eq. 0 ) then
          if( karg .eq. 1 ) then
              write(6,2)'Input filename? '
              read(5,1) filein
          else
              go to 20
          endif
      else
          if( karg .gt. narg ) go to 20
          call getarg( karg, filein )
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

  20  write(6,2)'Do you wish to read a new data file y/(n)? '
      read(5,1,end=30) chr
      if( chr .eq. 'y' .or. chr .eq.'Y' ) then
          write(6,2)'Input filename? '
          read(5,1) filein
          go to 10
      endif

  30  lstop = .true.
      return

      end
