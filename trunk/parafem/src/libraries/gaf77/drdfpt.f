c  *********************************************************************
c  *                                                                   *
c  *                      Function drdfpt                              *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, Princeton, Feb. 18, 1988.
c
c  PURPOSE  reads a real value from standard input (with default)
c
c  Reads a real value from the terminal and assigns a default value if
c  no value is entered;
c
c      drdfpt     is the real value returned by the function
c      default    is the supplied default value
c
c---------------------------------------------------------------------------
      real*8 function drdfpt( default )
      implicit real*8 (a-h,o-z)
      character*20 input
c
   1  format(a)

      read(5,1) input

      if( lnblnk(input) .eq. 0 ) then
         drdfpt = default
      else
         read(input,*) drdfpt
      endif

      return
      end
