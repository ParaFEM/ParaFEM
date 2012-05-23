c  *********************************************************************
c  *                                                                   *
c  *                     integer function nchoos                       *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.0
c  Written by Gordon A. Fenton, TUNS, Thu May 26 09:28:15 1994
c
c  PURPOSE  returns number of combinations: n choose m
c
c-------------------------------------------------------------------------
      integer function nchoos( n, m )
      data ibig/2147483647/
      data half/0.5/, one/1.0/, ierr/0/

   1  format(a,i5,a,i5,a)

      nn = max0( m, n-m )
c					compute numerator (and avoid overflow)
      it = 1
      id = 1
      do 10 i = n, nn+1, -1
         if( it .gt. ibig/i ) go to 20
         it = it*i
         id = id*(i-nn)
  10  continue
      nchoos = it/id
      return
c					too big, use floating point
  20  f = one
      do 30 i = n, nn+1, -1
         f = f*float(i)/float(i-nn)
  30  continue
c					at most 20 error messages
      if( f .gt. float(ibig) ) then
         if( ierr .lt. 20 ) then
         write(6,1)'Can''t represent choose(',n,',',m,') as an integer'
         write(6,1)'Returning nchoos = 1'
         ierr = ierr + 1
         endif
         nchoos = 1
      else
         nchoos = int( f + half )
      endif

      return
      end
