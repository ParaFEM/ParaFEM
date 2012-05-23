c  ***********************************************************************
c  *                                                                     *
c  *                             Function dcmplx                         *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, 1991
c
c  PURPOSE  to convert two double precision arguments into a single double
c           precision complex value.
c
c  Arguments to this routine are as follows;
c
c    DR     real value containing the real part of the number. (input)
c
c    DI     real value containing the imaginary part of the number. (input)
c
c-----------------------------------------------------------------------------
      complex*16 function dcmplx( DR, DI )
      complex*16 cmp
      real*8 DR, DI, rcmp(2)
      equivalence (cmp,rcmp(1))

      rcmp(1) = DR
      rcmp(2) = DI
      dcmplx  = cmp

      return
      end
