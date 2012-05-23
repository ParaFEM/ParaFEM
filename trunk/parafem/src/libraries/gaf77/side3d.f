c  *********************************************************************
c  *                                                                   *
c  *                        subroutine side3d                          *
c  *                                                                   *
c  *********************************************************************
c  Mixed Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Aug. 2, 1993
c  Latest Update: Feb. 22, 1994
c
c  PURPOSE  creates the parameter matrices required for the side cell
c           subdivisions of LAS3G.
c
c  Requires:
c   1) from libGAFsim:	DSIFA, DSISL, DCHOL2, DAXPY, DSWAP, DDOT, IDAMAX
c-----------------------------------------------------------------------------
      subroutine side3d(R,ir,B,ib,S,is,CS,AS,iout,tol)
      real CS(28,*), AS(18,7,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RS(18,18), DA(18), BB(7,7)
      integer ic(18,6), indx(18)
      data ic/10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
     >         4, 5, 6, 7, 8, 9,13,14,15,16,17,18,22,23,24,25,26,27,
     >         2, 3, 5, 6, 8, 9,11,12,14,15,17,18,20,21,23,24,26,27,
     >         1, 2, 4, 5, 7, 8,10,11,13,14,16,17,19,20,22,23,25,26,
     >         1, 2, 3, 4, 5, 6,10,11,12,13,14,15,19,20,21,22,23,24,
     >         1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18/

   1  format(a,a,a)
   2  format(a,e13.4)
c							for each side...
      do 70 ns = 1, 6
c							extract R
         do 10 j = 1, 18
            do 10 i = 1, j
               RS(i,j) = R(ic(i,ns), ic(j,ns))
  10     continue
c							factorize R
         call dsifa(RS,18,18,indx,ierr)
         if( ierr .ne. 0 ) then
            write(iout,1)'Error: unable to factorize side covariance mat
     >rix in SIDE3D.'
            stop
         endif
c							make a copy of S
         do 50 j = 1, 7
            do 20 i = 1, 18
               DA(i) = S(ic(i,ns),j)
  20        continue
c							and solve for A
            call dsisl(RS,18,18,indx,DA)
c							store in real*4
            do 30 i = 1, 18
               AS(i,j,ns) = DA(i)
  30        continue
c							update B
            do 40 i = 1, j
               BB(i,j) = B(i,j)
               do 40 k = 1, 18
                  BB(i,j) = BB(i,j) - S(ic(k,ns),i)*DA(k)
  40        continue
  50     continue
c							Cholesky Decomposition
         call dchol2( BB, 7, 7, rerr )
         if( rerr .gt. tol ) then
            write(iout,1)'Warning: Cholesky decomposition of side covari
     >ance matrix BB'
            write(iout,2)'         has maximum relative error of ',rerr
         endif
c							store in real*4
         ii = 0
         do 60 j = 1, 7
         do 60 i = 1, j
            ii = ii + 1
            CS(ii,ns) = BB(i,j)
  60     continue
  70  continue

      return
      end
