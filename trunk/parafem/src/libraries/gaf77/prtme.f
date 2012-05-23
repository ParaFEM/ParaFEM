c  *********************************************************************
c  *                                                                   *
c  *                         Subroutine prtme                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 3.11
c  Written by Gordon A. Fenton, 1989
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to print a matrix using e format
c
c  This routine prints an n x m matrix to unit IOUT using the e format. See
c  prtmt for a routine which uses the f format.
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
c  3.1	eliminated unused variable `one' (Dec 5/96)
c  3.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c--------------------------------------------------------------------------
      subroutine prtme( iout, U, iu, n, m, title )
      parameter (nCOL = 80)
      dimension U(iu,*)
      character*(*) title
      character*16 fmt
      data zero/0.0/

   1  format(a)
   2  format()
   3  format('(',i2,'e',i1,'.',i1,')')
   4  format('(',i2,'e',i2,'.',i1,')')

      write(iout,2)
      write(iout,1) title
c					check data to find sign
      isyn = 0
      do 10 j = 1, m
      do 10 i = 1, n
         if( U(i,j) .lt. zero ) isyn = 1
  10  continue
c					actual field width
      nwd = nCOL/m
c					available decimal positions
      nr  = nwd - 7 - isyn
c					will our numbers fit (3 digits min)?
      if( nr .lt. 3 ) then
         nr  = 3
         nwd = 9 + isyn
      elseif( nr .gt. 7 ) then
         nr  = 7
         nwd = 13 + isyn
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
