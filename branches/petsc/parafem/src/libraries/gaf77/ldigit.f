c  *********************************************************************
c  *                                                                   *
c  *                  logical function ldigit                          *
c  *                                                                   *
c  *********************************************************************
c  Logical Version 1.0
c  Written by Gordon A. Fenton, Dalhousie University, Oct  1, 2001
c  Latest Update: Oct  1, 2001
c
c  PURPOSE  returns true if the provided character is a digit in the range
c           0-9 or a decimal '.'
c
c  DESCRIPTION
c
c  This function checks to see if the provided character, c, is a digit in
c  the range 0 to 9 or a decimal point (.). It does so by checking the
c  integer index of the character, which assumes that the character set
c  follows the ASCII indexing.
c
c  ARGUMENTS
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      logical function ldigit(c)
      character*1 c
      integer ic

      ic     = ichar(c)
      ldigit = ((ic .ge. 48) .and. (ic .le. 57)) .or. (ic .eq. 46)

      return
      end
