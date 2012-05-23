c  ***********************************************************************
c  *                                                                     *
c  *                            Function small                             *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.0
c  Written by G.A. Fenton, TUNS, 1991
c
c  PURPOSE   return the machine epsilon (unit roundoff)
c
c  This routine computes the smallest number that can be added to 1.0
c  without getting lost (the machine epsilon). It takes no arguments.
c----------------------------------------------------------------------------
      real function small()
      data one/1.0/, half/0.5/, two/2.0/

      small = one
10    if( (small + one) .GT. one ) then
         small = half*small
         go to 10
      endif

      small = two*small

      return
      end

