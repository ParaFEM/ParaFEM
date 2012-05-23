c  *********************************************************************
c  *                                                                   *
c  *                       Subroutine argin                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, Princeton, July 22, 1988.
c  Latest Update: Jul 27, 1999
c
c  PURPOSE  to obtain a character string filename from the argument list
c
c  Obtains the optional argument which is assumed to be the name of the
c  input file. This routine can be called a number of times and it will get
c  succeeding arguments and assign them to filein. When the argument list
c  is exhausted, the logical flag `lstop' is set to true. Arguments to
c  this routine are as follows;
c
c    filein   character string which on output will contain the input
c             file name as obtained from the command line arguments. (output)
c
c     lstop   logical flag which is set to true if no further filenames
c             are specified (by the user or on the command line). (output)
c
c  REVISION HISTORY:
c  1.1	now getarg is an intrinsic (Jul 27/99)
c----------------------------------------------------------------------------
      subroutine argin( filein, lstop )
      logical lstop
      character*(*) filein
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
      return

  20  lstop = .true.
      return

      end
