c  ********************************************************************
c  *                                                                  *
c  *                       subroutine openfl                          *
c  *                                                                  *
c  ********************************************************************
c  Integer Version 1.1
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jul 27, 1999
c
c  PURPOSE  to open a data file given by an argument on the command line
c
c  The parameter `karg' is provided because not all compilers number
c  their command line arguments the same way (ie HP-UX uses karg = 1).
c
c  REVISION HISTORY:
c  1.1	now getarg is an intrinsic (Jul 27/99)
c---------------------------------------------------------------------------
      subroutine openfl( iunit )
      parameter (karg = 0)
      character*64 datfil
      logical found

   1  format(a)
   2  format(a,$)
c						get the data file name
      narg = iargc()
      if( narg .lt. (karg+1) ) then
         write(6,2)'Please enter the data file name: '
         read(5,1) datfil
      else
         call getarg( (karg+1), datfil )
      endif
c						open the data file
      inquire( file=datfil, exist=found)
      if( .not. found ) then
         write(6,1)'Data file not found: '//datfil
         write(6,1)'Please create or check spelling.'
         stop
      endif
      open( iunit, file = datfil )

      return
      end
