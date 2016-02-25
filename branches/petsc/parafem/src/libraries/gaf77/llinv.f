c  ***********************************************************************
c  *                                                                     *
c  *                          Subroutine llinv                           *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the inverse of a lower triangular matrix
c
c  This routine takes a lower triangular matrix TL and computes its inverse
c  (which is also lower triangular) TI. Arguments to the routine are
c
c   TL   real array of size at least N x N which contains the elements of
c        the lower triangular matrix to be inverted (the elements above the
c        diagonal of TL are zeroed in this routine). (input/output)
c
c   TI   real array of size at least N x N which will contain the inverse
c        of TL in its lower triangle. Elements above the diagonal are set
c        to zero. (output)
c
c   IM   integer giving the leading dimension of TL and TI as defined in the
c        calling routine. (input)
c
c   N    integer giving the size of the matrices TL and TI. (input)
c
c ierr   integer flag which is set to zero if all goes well. If one of the
c        diagonal elements of TL is zero, then the matrix is singular and
c        ierr is set to -1. (output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine llinv( TL, TI, IM, N, ierr )
      dimension TL(IM,*), TI(IM,*)
      data zero/0.0/, one/1.0/
c						assume all will not go well
      ierr = -1
      do 30 i = 1, N
         if( TL(i,i) .eq. zero ) return
         TI(i,i) = one/TL(i,i)
         do 20 j = 1, i - 1
            s = TL(i,j)*TI(j,j)
            do 10 k = j + 1, i - 1
               s = s + TL(i,k)*TI(k,j)
  10        continue
            TI(i,j) = -s * TI(i,i)
            TI(j,i) = zero
            TL(j,i) = zero
  20     continue
  30  continue
c						all went well, set flag
      ierr = 0

      return
      end
