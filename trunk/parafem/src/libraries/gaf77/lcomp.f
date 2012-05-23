c  *********************************************************************
c  *                                                                   *
c  *                    Logical Function lcomp                         *
c  *                                                                   *
c  *********************************************************************
c  Logical Version 1.0
c  Written by Gordon A. Fenton, Princeton, June 25, 1988.
c
c  PURPOSE compare two strings and return true if str1 is contained in str2.
c
c----------------------------------------------------------------------------
      logical function lcomp( str1, str2 )
      character*(*) str1, str2

      k1 = lnblnk( str1 ) - 1
      k2 = lnblnk( str2 ) - k1

      if( k2 .lt. 1 ) then
          lcomp = .false.
          return
      endif

      do 10 i = 1, k2
         if( str2(i:i+k1) .eq. str1 ) then
             lcomp = .true.
             return
         endif
  10  continue

      lcomp = .false.
      return
      end
