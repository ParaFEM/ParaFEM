c  *********************************************************************
c  *                                                                   *
c  *                       subroutine uuinv                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the inverse of an upper triangular matrix
c
c  This routine takes an upper triangular matrix TU and computes its inverse
c  (which is also upper triangular) TI. Arguments to the routine are
c
c   TU   real array of size at least N x N which contains the elements of
c        the upper triangular matrix to be inverted. (input)
c
c   TI   real array of size at least N x N which will contain the inverse
c        of TU in its upper triangle. Elements below the diagonal are set
c        to zero. (output)
c
c   IM   integer giving the leading dimension of TU and TI as defined in the
c        calling routine. (input)
c
c   N    integer giving the size of the matrices TU and TI. (input)
c
c ierr   integer flag which is set to zero if all goes well. If one of the
c        diagonal elements of TU is zero, then the matrix is singular and
c        ierr is set to -1. (output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine uuinv( TU, TI, IM, N, ierr )
      dimension TU(IM,*), TI(IM,*)
      data zero/0.0/, one/1.0/
c						assume all will not go well
      ierr = -1
      do 30 j = 1, N
         if( TU(j,j) .eq. zero ) return
         TI(j,j) = one/TU(j,j)
         do 20 i = j-1, 1, -1
            s = TU(i,j)*TI(j,j)
            do 10 k = i + 1, j - 1
               s = s + TU(i,k)*TI(k,j)
  10        continue
            TI(i,j) = -s * TI(i,i)
            TI(j,i) = zero
  20     continue
  30  continue
c						all went well, set flag
      ierr = 0

      return
      end
