c  *********************************************************************
c  *                                                                   *
c  *                        Subroutine dprtmt                          *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 3.01
c  Written by Gordon A. Fenton, 1989
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to print a matrix using f format
c
c  This routine prints an n x m matrix to unit IOUT using the f format. See
c  prtme for a routine which uses the e format.
c  Arguments are as follows;
c
c iout   fortran unit number (assumed already opened). (input)
c
c    U   real array containing the matrix to be printed. (input)
c
c   iu   leading dimension of the array U as specified in the dimension
c        statement of the calling routine. (input)
c
c    n   column dimension of U (number of rows). (input)
c
c    m   row dimension of U (number of columns). (input)
c
c title  character string giving the title to be printed. (input)
c
c The parameter nCOL is the desired number of columns within which to format
c the rows of the matrix. The actual number of columns used will vary,
c depending on the number of columns in U.
c
c  REVISION HISTORY:
c  3.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c--------------------------------------------------------------------------
      subroutine dprtmt( iout, U, iu, n, m, title )
      parameter (nCOL = 80)
      implicit real*8 (a-h,o-z)
      dimension U(iu,*)
      character*(*) title
      character fmt*16, rfmt*26
      data zero/0.d0/, one/1.d0/
      abs(x)    = dabs(x)
      alog10(x) = dlog10(x)
      int(x)    = idint(x)

   1  format(a)
   2  format()

      write(iout,2)
      write(iout,1) title
c					check data to find maximum
      umax = abs(U(1,1))
      isyn = 0
      do 10 j = 1, m
      do 10 i = 1, n
         u0 = abs(U(i,j))
         if( u0 .gt. umax ) umax = u0
         if( U(i,j) .lt. zero ) isyn = 1
  10  continue
c					actual field width
      nwd = nCOL/m
c					available field width (less decimal)
      nfld = nwd - 2 - isyn
c					4 digits minimum
      if( nfld .lt. 4 ) then
         nfld = 4
         nwd  = 6 + isyn
      endif
c					required digits to left of decimal
      if( umax .ge. one ) then
         nl = 1 + int( alog10(umax) )
      else
         nl = 1
      endif
c					will our numbers fit?
      if( nl .gt. nfld ) then
         nfld = nl
         nwd  = nl + 2 + isyn
      endif
c					number of digits to right of decimal
      nr = nfld - nl
      if( nr .gt. 16 ) then
         nfld = nfld - (nr - 16)
         nwd  = nfld + 2 + isyn
         nr   = 16
      endif
c					now set up the format
      rfmt = '(''('',i2,''f'',i2,''.'',i2,'')'')'
      if( nwd .lt. 10 ) rfmt(14:14) = '1'
      if( nr  .lt. 10 ) rfmt(21:21) = '1'
      write(fmt,rfmt) m, nwd, nr
c					and print the matrix
      do 20 i = 1, n
         write(iout,fmt) ( u(i,j), j = 1, m )
  20  continue

      return
      end
