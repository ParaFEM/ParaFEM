c  *********************************************************************
c  *                                                                   *
c  *                        subroutine edge3d                          *
c  *                                                                   *
c  *********************************************************************
c  Mixed Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, Aug. 2, 1993
c  Latest Update: Feb. 22, 1994
c
c  PURPOSE  creates the parameter matrices required for the edge cell
c           subdivisions of LAS3G.
c
c  Requires:
c   1) from libGAFsim:	DSIFA, DSISL, DCHOL2, DAXPY, DSWAP, DDOT, IDAMAX
c-----------------------------------------------------------------------------
      subroutine edge3d(R,ir,B,ib,S,is,CE,AE,iout,tol)
      real CE(28,*), AE(12,7,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RS(12,12), DA(12), BB(7,7)
      integer ic(12,12), indx(12)
      data ic/13,14,15,16,17,18,22,23,24,25,26,27,
     >        11,12,14,15,17,18,20,21,23,24,26,27,
     >        10,11,13,14,16,17,19,20,22,23,25,26,
     >        10,11,12,13,14,15,19,20,21,22,23,24,
     >         5, 6, 8, 9,14,15,17,18,23,24,26,27,
     >         4, 5, 7, 8,13,14,16,17,22,23,25,26,
     >         2, 3, 5, 6,11,12,14,15,20,21,23,24,
     >         1, 2, 4, 5,10,11,13,14,19,20,22,23,
     >         4, 5, 6, 7, 8, 9,13,14,15,16,17,18,
     >         2, 3, 5, 6, 8, 9,11,12,14,15,17,18,
     >         1, 2, 4, 5, 7, 8,10,11,13,14,16,17,
     >         1, 2, 3, 4, 5, 6,10,11,12,13,14,15/

   1  format(a,a,a)
   2  format(a,e13.4)
c							for each edge...
      do 70 ne = 1, 12
c							extract R
         do 10 j = 1, 12
            do 10 i = 1, j
               RS(i,j) = R(ic(i,ne), ic(j,ne))
  10     continue
c							factorize R
         call dsifa(RS,12,12,indx,ierr)
         if( ierr .ne. 0 ) then
            write(iout,1)'Error: unable to factorize edge covariance mat
     >rix in EDGE3D.'
            stop
         endif
c							make a copy of S
         do 50 j = 1, 7
            do 20 i = 1, 12
               DA(i) = S(ic(i,ne),j)
  20        continue
c							and solve for A
            call dsisl(RS,12,12,indx,DA)
c							store in real*4
            do 30 i = 1, 12
               AE(i,j,ne) = DA(i)
  30        continue
c							update B
            do 40 i = 1, j
               BB(i,j) = B(i,j)
               do 40 k = 1, 12
                  BB(i,j) = BB(i,j) - S(ic(k,ne),i)*DA(k)
  40        continue
  50     continue
c							Cholesky Decomposition
         call dchol2( BB, 7, 7, rerr )
         if( rerr .gt. tol ) then
            write(iout,1)'Warning: Cholesky decomposition of edge covari
     >ance matrix BB'
            write(iout,2)'         has maximum relative error of ',rerr
         endif
c							store in real*4
         ii = 0
         do 60 j = 1, 7
         do 60 i = 1, j
            ii = ii + 1
            CE(ii,ne) = BB(i,j)
  60     continue
  70  continue

      return
      end
