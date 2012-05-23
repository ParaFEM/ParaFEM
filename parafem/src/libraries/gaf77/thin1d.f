c  *********************************************************************
c  *                                                                   *
c  *                        subroutine thin1d                          *
c  *                                                                   *
c  *********************************************************************
c  Mixed Precision Version 1.11
c  Written by Gordon A. Fenton, TUNS, Dec. 8, 1993
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  creates the parameter matrices required for the interior cell
c           subdivisions of LAS2G for the case k1 or k2 = 1 and for LAS3G
c           for k1=k2=1, k1=k3=1, or k2=k3=1.
c
c  Requires:
c   1) from libGAFsim:	DSIFA, DSISL, DCHOL2, DAXPY, DSWAP, DDOT, IDAMAX
c
c  i1, i2, and i3 are the indexes into the global covariance matrix R to be
c  extracted to form the local covariance matrix.
c  n is the number of cells in the subdivision less 1. This is equal to
c  3 for LAS2G and to 7 for LAS3G.
c
c  REVISION HISTORY:
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine thin1d(R,ir,B,ib,S,is,AS,ias,AI,iai,CS,CI,
     >                  i1,i2,i3,n,iout,tol)
      real AS(ias,n,*), AI(iai,*), CS(*), CI(*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RT(2,2), RI(3,3), DA(3), B1(7,7), B2(7,7)
      integer indT(2), indI(3)

   1  format(a,a,a)
   2  format(a,e13.4)
c							extract R
      RT(1,1) = R(i2,i2)
      RT(1,2) = R(i2,i3)
      RT(2,2) = R(i3,i3)

      RI(1,1) = R(i1,i1)
      RI(1,2) = R(i1,i2)
      RI(1,3) = R(i1,i3)
      RI(2,2) = R(i2,i2)
      RI(2,3) = R(i2,i3)
      RI(3,3) = R(i3,i3)
c							factorize RT (edges)
      call dsifa(RT,2,2,indT,ierr)
      if( ierr .ne. 0 ) then
         write(iout,1)'Error: unable to factorize covariance matrix RT i
     >n THIN1D.'
         stop
      endif
c							factorize RI (interior)
      call dsifa(RI,3,3,indI,ierr)
      if( ierr .ne. 0 ) then
         write(iout,1)'Error: unable to factorize covariance matrix RI i
     >n THIN1D.'
         stop
      endif

      do 30 j = 1, n
         DA(1) = S(i2,j)
         DA(2) = S(i3,j)
         call dsisl( RT, 2, 2, indT, DA )
         AS(1,j,1) = DA(1)
         AS(2,j,1) = DA(2)
         do 10 i = 1, j
            B1(i,j) = B(i,j) - S(i2,i)*DA(1) - S(i3,i)*DA(2)
  10     continue

         DA(1) = S(i1,j)
         DA(2) = S(i2,j)
         DA(3) = S(i3,j)
         call dsisl( RI, 3, 3, indI, DA )
         AI(1,j) = DA(1)
         AI(2,j) = DA(2)
         AI(3,j) = DA(3)
         do 20 i = 1, j
            B2(i,j) = B(i,j) - S(i1,i)*DA(1) - S(i2,i)*DA(2)
     >                       - S(i3,i)*DA(3)
  20     continue
  30  continue
c							Cholesky Decomposition
      call dchol2( B1, 7, n, rerr )
      if( rerr .gt. tol ) then
         write(iout,1)'Warning: Cholesky decomposition of THIN1D covaria
     >nce matrix B1'
         write(iout,2)'         has maximum relative error of ',rerr
      endif

      call dchol2( B2, 7, n, rerr )
      if( rerr .gt. tol ) then
         write(iout,1)'Warning: Cholesky decomposition of THIN1D covaria
     >nce matrix B2'
         write(iout,2)'         has maximum relative error of ',rerr
      endif
c							store in real*4
      ii = 0
      do 40 j = 1, n
      do 40 i = 1, j
         ii = ii + 1
         CS(ii) = B1(i,j)
         CI(ii) = B2(i,j)
  40  continue

      return
      end
