c  *******************************************************************
c  *                                                                 *
c  *                        subroutine sisl                          *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 08/14/78 1.01
c  Written by James Bunch, Univ. Calif. San Diego, Argonne Nat. Lab.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   solves the single precision symmetric system [A]{X} = {B}
c            using the factors computed by sifa (LINPACK)
c
c     On entry
c
c        A       single precision(lda,n)
c                the output from sifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from sifa.
c
c        B       single precision(n)
c                the right hand side vector.
c
c     on return
c
c        B       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  sico  has set rcond .eq. 0.0
c        or  sifa  has set ierr .ne. 0  .
c
c     subroutines and functions
c
c     blas axpy, addot
c     fortran iabs
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-------------------------------------------------------------------------
      subroutine sisl(a,lda,n,kpvt,b)
      real a(lda,*), b(*)
      real ak, akm1, bk, bkm1, denom, temp
      integer lda, n, kpvt(*)
      integer k, kp
c							work backwards
      k = n
   10 if (k .ne. 0) then
         if (kpvt(k) .ge. 0) then
c							1 x 1 pivot block.
            if (k .ne. 1) then
               kp = kpvt(k)
               if (kp .ne. k) then
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
               endif
               call axpy(k-1,b(k),a(1,k),1,b(1),1)
            endif
            b(k) = b(k)/a(k,k)
            k = k - 1
         else
c							2 x 2 pivot block.
            if (k .ne. 2) then
               kp = iabs(kpvt(k))
               if (kp .ne. k - 1) then
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
               endif
               call axpy(k-2,b(k),a(1,k),1,b(1),1)
               call axpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
            endif
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
         endif
         go to 10
      endif
c							now work forwards
      k = 1
   90 if (k .le. n) then
         if (kpvt(k) .ge. 0) then
c							1 x 1 pivot block.
            if (k .ne. 1) then
               b(k) = b(k) + addot(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .ne. k) then
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
               endif
            endif
            k = k + 1
         else
c							2 x 2 pivot block.
            if (k .ne. 1) then
               b(k) = b(k) + addot(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + addot(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .ne. k) then
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
               endif
            endif
            k = k + 2
         endif
         go to 90
      endif

      return
      end
