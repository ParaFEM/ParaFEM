c  *******************************************************************
c  *                                                                 *
c  *                        subroutine sifa                          *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 08/14/78 1.01
c  Written by James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  factors a single precision symmetric matrix by elimination
c           with symmetric pivoting (LINPACK)
c
c     To solve  A*X = B ,          follow sifa by dsisl.
c     To compute  inverse(A)*C ,   follow sifa by dsisl.
c
c     On entry
c
c        A       single precision(lda,n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        A       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        ierr    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that sisl may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas axpy, aswap,iamax
c     fortran abs,amax1,sqrt
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-------------------------------------------------------------------------
      subroutine sifa(a,lda,n,kpvt,ierr)
      dimension a(lda,*)
      integer kpvt(*)
      logical swap
      data zero/0.0/, one/1.0/, eight/8.0/, svtn/17.0/
c						initialize
c						choose pivot block size.
      alpha = (one + sqrt(svtn))/eight
      ierr  = 0
c						loop on k, from n to 1.
      k = n
c						leave the loop if k=0 or k=1.
  10  if (k .eq. 0) return
      if (k .eq. 1) then
         kpvt(1) = 1
         if (a(1,1) .eq. zero) ierr = 1
         return
      endif
      km1 = k - 1
      absakk = abs(a(k,k))
c						largest element in col k
      imax = iamax(k-1,a(1,k),1)
      colmax = abs(a(imax,k))
      if (absakk .ge. alpha*colmax) then
         kstep = 1
         swap = .false.
      else
c						largest element in row imax
         rowmax = zero
         imaxp1 = imax + 1
         do 20 j = imaxp1, k
            rowmax = amax1(rowmax,abs(a(imax,j)))
  20     continue
         if (imax .ne. 1) then
            jmax = iamax(imax-1,a(1,imax),1)
            rowmax = amax1(rowmax,abs(a(jmax,imax)))
         endif
         if (abs(a(imax,imax)) .ge. alpha*rowmax) then
            kstep = 1
            swap = .true.
         else
            if (absakk .ge. alpha*colmax*(colmax/rowmax)) then
               kstep = 1
               swap = .false.
            else
               kstep = 2
               swap = imax .ne. km1
            endif
         endif
      endif
      if (amax1(absakk,colmax) .eq. zero) then
c						column k is zero.  set ierr
         kpvt(k) = k
         ierr = k
      elseif (kstep .eq. 1) then
c						1 x 1 pivot block.
         if (swap) then
            call aswap(imax,a(1,imax),1,a(1,k),1)
            do 30 jj = imax, k
               j = k + imax - jj
               t = a(j,k)
               a(j,k) = a(imax,j)
               a(imax,j) = t
  30        continue
         endif
c						perform the elimination.
         do 40 jj = 1, km1
            j = k - jj
            amk = -a(j,k)/a(k,k)
            t = amk
            call axpy(j,t,a(1,k),1,a(1,j),1)
            a(j,k) = amk
  40     continue
c	     				set the pivot array.
         kpvt(k) = k
         if (swap) kpvt(k) = imax
      else
c						2 x 2 pivot block.
         if (swap) then
            call aswap(imax,a(1,imax),1,a(1,k-1),1)
            do 50 jj = imax, km1
               j = km1 + imax - jj
               t = a(j,k-1)
               a(j,k-1) = a(imax,j)
               a(imax,j) = t
  50        continue
            t = a(k-1,k)
            a(k-1,k) = a(imax,k)
            a(imax,k) = t
         endif
c						perform the elimination.
         km2 = k - 2
         if (km2 .ne. 0) then
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            denom = one - ak*akm1
            do 60 jj = 1, km2
               j = km1 - jj
               bk = a(j,k)/a(k-1,k)
               bkm1 = a(j,k-1)/a(k-1,k)
               amk = (akm1*bk - bkm1)/denom
               amkm1 = (ak*bkm1 - bk)/denom
               t = amk
               call axpy(j,t,a(1,k),1,a(1,j),1)
               t = amkm1
               call axpy(j,t,a(1,k-1),1,a(1,j),1)
               a(j,k) = amk
               a(j,k-1) = amkm1
  60        continue
         endif
c            				set the pivot array.
         kpvt(k) = 1 - k
         if (swap) kpvt(k) = -imax
         kpvt(k-1) = kpvt(k)
      endif
      k = k - kstep
      go to 10

      end
