c  ***********************************************************************
c  *                                                                     *
c  *                              Function cmplx                         *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, 1992
c
c  PURPOSE  to convert two real arguments into a complex value
c
c  Arguments to this routine are as follows;
c
c    DR     real value containing the real part of the number. (input)
c
c    DI     real value containing the imaginary part of the number. (input)
c
c-----------------------------------------------------------------------------
      complex function cmplx( DR, DI )
      complex cmp
      real DR, DI, rcmp(2)
      equivalence (cmp,rcmp(1))

      rcmp(1) = DR
      rcmp(2) = DI
      cmplx  = cmp

      return
      end
