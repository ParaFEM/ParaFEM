c  *********************************************************************
c  *                                                                   *
c  *                         Subroutine prtmt                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 3.01
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
      subroutine prtmt( iout, U, iu, n, m, title )
      parameter (nCOL = 80)
      dimension U(iu,*)
      character*(*) title
      character*16 fmt
      data zero/0.0/, one/1.0/

   1  format(a)
   2  format()
   3  format('(',i2,'f',i1,'.',i1,')')
   4  format('(',i2,'f',i2,'.',i1,')')

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
c					3 digits minimum
      if( nfld .lt. 3 ) then
         nfld = 3
         nwd  = 5 + isyn
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
      if( nr .gt. 7 ) then
         nfld = nfld - (nr - 7)
         nwd  = nfld + 2 + isyn
         nr   = 7
      endif
c					now set up the format
      if( nwd .lt. 10 ) then
         write(fmt,3) m, nwd, nr
      else
         write(fmt,4) m, nwd, nr
      endif
c					and print the matrix
      do 20 i = 1, n
         write(iout,fmt) ( u(i,j), j = 1, m )
  20  continue

      return
      end
