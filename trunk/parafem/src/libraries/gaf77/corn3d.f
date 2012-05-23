c  *********************************************************************
c  *                                                                   *
c  *                        subroutine corn3d                          *
c  *                                                                   *
c  *********************************************************************
c  Mixed Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Aug. 2, 1993
c  Latest Update: Feb. 22, 1994
c
c  PURPOSE  creates the parameter matrices required for the corner cell
c           subdivisions of LAS3G.
c
c  Requires:
c   1) from libGAFsim:	DSIFA, DSISL, DCHOL2, DAXPY, DSWAP, DDOT, IDAMAX
c-----------------------------------------------------------------------------
      subroutine corn3d(R,ir,B,ib,S,is,CC,AC,iout,tol)
      real CC(28,*), AC(8,7,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RC(8,8), DA(8), BB(7,7)
      integer ic(8,8), indx(8)
      data ic/14,15,17,18,23,24,26,27,
     >        13,14,16,17,22,23,25,26,
     >        11,12,14,15,20,21,23,24,
     >        10,11,13,14,19,20,22,23,
     >         5, 6, 8, 9,14,15,17,18,
     >         4, 5, 7, 8,13,14,16,17,
     >         2, 3, 5, 6,11,12,14,15,
     >         1, 2, 4, 5,10,11,13,14/

   1  format(a,a,a)
   2  format(a,e13.4)
c							extract R
      do 10 j = 1, 8
         do 10 i = 1, j
            RC(i,j) = R(ic(i,8), ic(j,8))
  10  continue
c							factorize R
      call dsifa(RC,8,8,indx,ierr)
      if( ierr .ne. 0 ) then
         write(iout,1)'Error: unable to factorize corner covariance matr
     >ix in CORN3D.'
         stop
      endif
c							for each corner...
      do 70 nc = 1, 8
c							make a copy of S
         do 50 j = 1, 7
            do 20 i = 1, 8
               DA(i) = S(ic(i,nc),j)
  20        continue
c							and solve for A
            call dsisl(RC,8,8,indx,DA)
c							store in real*4
            do 30 i = 1, 8
               AC(i,j,nc) = DA(i)
  30        continue
c							update B
            do 40 i = 1, j
               BB(i,j) = B(i,j)
     >                 - S(ic(1,nc),i)*DA(1) - S(ic(2,nc),i)*DA(2)
     >                 - S(ic(3,nc),i)*DA(3) - S(ic(4,nc),i)*DA(4)
     >                 - S(ic(5,nc),i)*DA(5) - S(ic(6,nc),i)*DA(6)
     >                 - S(ic(7,nc),i)*DA(7) - S(ic(8,nc),i)*DA(8)
  40        continue
  50     continue
c							Cholesky Decomposition
         call dchol2( BB, 7, 7, rerr )
         if( rerr .gt. tol ) then
            write(iout,1)'Warning: Cholesky decomposition of corner cova
     >riance matrix BB'
            write(iout,2)'         has maximum relative error of ',rerr
         endif
c							store in real*4
         ii = 0
         do 60 j = 1, 7
         do 60 i = 1, j
            ii = ii + 1
            CC(ii,nc) = BB(i,j)
  60     continue
  70  continue

      return
      end
