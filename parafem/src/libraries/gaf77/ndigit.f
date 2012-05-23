c  *********************************************************************
c  *                                                                   *
c  *                      integer function ndigit                      *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.0
c  Written by Gordon A. Fenton, TUNS, Thu Sep  7 16:02:32 1995
c
c  PURPOSE  returns the number of digits in the supplied integer value
c
c  This function returns the number of digits in the supplied integer value.
c  If the value is 0, the number of digits return is 1. For negative numbers,
c  the number of digits includes the `-' sign, thus -43 is interpreted as
c  having 3 digits. This function is especially useful in modifying Fortran
c  format statements to just include the field width necessary to print the
c  number. For example, if you are trying to print the integer variable ival
c  with just enough width, try
c
c      character fmt*16, dig(9)*1
c      data dig/'1','2','3','4','5','6','7','8','9'/
c       .
c       .
c      fmt = '(i9)'
c      nd  = ndigit(ival,'int')
c      if( nd .gt. 9 ) then
c         write(6,'(a)')'Can only print integers up to 9 digits in length'
c         stop
c      endif
c      fmt(3:3) = dig(nd)
c      write(6,fmt) ival
c
c  where the limitation of up to 9 digits is only in the above example. If
c  you desire the number of digits to the left of the decimal place in a
c  real number, try ndigit( int(val) ).
c-------------------------------------------------------------------------
      integer function ndigit(ival)
      integer ival

      if( ival .eq. 0 ) then
         ndigit = 1
      elseif( ival .lt. 0 ) then
         ndigit = 2 + int( alog10( 0.0001 + float(-ival) ) )
      else
         ndigit = 1 + int( alog10( 0.0001 + float( ival ) ) )
      endif

      return
      end
