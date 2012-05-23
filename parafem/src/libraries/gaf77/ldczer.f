c  ***********************************************************************
c  *                                                                     *
c  *                            Function ldczer                          *
c  *                                                                     *
c  ***********************************************************************
c  Logical/Double Precision/Complex Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to determine if the complex argument is algorithmically zero
c
c  This function accepts a complex argument (treated here as two real values)
c  and determines if they (ie. real and imaginary parts) are effectively
c  zero. The function returns true if this is the case. The argument is
c  considered to be effectively zero if both its components are less than
c  `small' which has been chosen rather arbitrarily here.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      logical function ldczer( c )
      real*8 c(*), dabs, small
      data small/1.d-60/

      if( dabs(c(1)) .lt. small .and. dabs(c(2)) .lt. small ) then
         ldczer = .true.
      else
         ldczer = .false.
      endif

      return
      end
