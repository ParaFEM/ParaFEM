c  ******************************************************************
c  *                                                                *
c  *                        Subroutine tridg                        *
c  *                                                                *
c  ******************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  reduces a real symmetrix matrix to symmetric tridiagonal form
c
c  This routine reduces a real symmetric matrix to a symmetric
c  tridiagonal matrix using and accumulating orthogonal similarity
c  transformations [Q]. This routine was obtained from EISPACK.
c
c  Input parameters are discussed as follows:
c
c   A    is a real two-dimensional variable of size A(IA,N).  A contains
c        the symmetric matrix of order N to be reduced to tridiagonal form.
c        Only the full upper triangle of the matrix need be supplied.
c
c   Q    is a real two-dimensional output variable of dimension IA x N which
c        will contain the orthogonal transformation matrix produced in the
c        reduction to the tridiagonal form.
c
c   DIAG is a real output one-dimensional variable of dimension at least N in
c        which the diagonal elements of the reduced matrix will be placed.
c
c   SD   is a real output one-dimensional variable of dimension at least N in
c        which the subdiagonal elements of the reduced matrix will be placed.
c        Since there are only N-1 of these, the element SD(1) is set to zero.
c
c   IA   is an integer input variable set equal to the column dimension
c        of the two-dimensional arrays A and Q as specified in the
c        DIMENSION statement of the calling routine.
c
c   N    is an integer set equal to the order of the matrix A.  N must
c        not be greater than IA.
c
c
c   NOTES:
c          - if the matrix has elements of widely varying magnitudes, the
c            smaller ones should be in the top left-hand corner.
c
c          - the array Q is not needed. The transformation matrix can be easily
c            written over the array A.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c------------------------------------------------------------------------------
      Subroutine tridg( A, Q, diag, sd, ia, n )
      dimension a(ia,*), diag(*), sd(*), Q(ia,*)
      data zero/0.0/, one/1.0/
c                                                 initialize Q
      do 10 i = 1,n
      do 10 j = 1,i
         Q(i,j) = A(j,i)
  10  continue

      if ( n .gt. 1 ) then
         do 80 ii = 2,n
            i = n + 2 - ii
            l = i - 1
            h = zero
            scale = zero
            if ( l .lt. 2 ) then
               sd(i)   = Q(i,l)
               diag(i) = h
            else
c                                                  scale row
               do 20 k = 1,l
                  scale = scale + abs(Q(i,k))
  20           continue

               if ( scale .eq. zero ) then
                  sd(i)   = Q(i,l)
                  diag(i) = h
               else
                  do 30 k = 1,l
                     Q(i,k) = Q(i,k) / scale
                     h      = h + Q(i,k)*Q(i,k)
  30              continue

                  f = Q(i,l)
                  g = -sign( sqrt(h), f )
                  sd(i) = scale * g
                  h = h - f*g
                  Q(i,l) = f - g
                  f = zero

                  do 60 j = 1,l
                     Q(j,i) = Q(i,j) / h
                     g = zero
c                                               form element of A*U
                     do 40 k = 1,j
                        g = g + Q(j,k)*Q(i,k)
  40                 continue
                     do 50 k = j+1,l
                        g = g + Q(k,j)*Q(i,k)
  50                 continue
c                                               form element of P
                     sd(j) = g / h
                     f = f + sd(j)*Q(i,j)
  60              continue

                  hh = f / ( h + h )
c                                               form reduced A
                  do 70 j = 1,l
                     f = Q(i,j)
                     g = sd(j) - hh*f
                     sd(j) = g
                     do 70 k = 1,j
                        Q(j,k) = Q(j,k) - f*sd(k) - g*Q(i,k)
  70              continue
               endif
            endif
  80     continue
      endif
c                                      accumulation of transformation matrix
      diag(1) = zero
      sd(1) = zero
      do 120 i = 1,n
         l = i - 1
         if ( diag(i) .ne. zero ) then
            do 100 j = 1,l
               g = zero
               do 90 k = 1,l
                  g = g + Q(i,k)*Q(k,j)
  90           continue
               do 100 k = 1,l
                  Q(k,j) = Q(k,j) - g*Q(k,i)
 100        continue
         endif
         diag(i) = Q(i,i)
         Q(i,i)  = one
         if ( l .ge. 1 ) then
            do 110 j = 1,l
               Q(i,j) = zero
               Q(j,i) = zero
  110       continue
         endif
  120 continue

      return
      end
