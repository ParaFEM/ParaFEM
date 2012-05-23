c  ***********************************************************************
c  *                                                                     *
c  *                            Function dsmall                          *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.0
c  Written by G.A. Fenton, TUNS, 1991
c
c  PURPOSE   return the machine epsilon (unit roundoff)
c
c  This routine computes the dsmallest number that can be added to 1.0
c  without getting lost (the machine epsilon). It takes no arguments.
c----------------------------------------------------------------------------
      real*8 function dsmall()
      implicit real*8 (a-h,o-z)
      data one/1.d0/, half/0.5d0/, two/2.d0/

      dsmall = one
10    if( (dsmall + one) .GT. one ) then
         dsmall = half*dsmall
         go to 10
      endif

      dsmall = two*dsmall

      return
      end

