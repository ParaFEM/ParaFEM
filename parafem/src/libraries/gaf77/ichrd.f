c  *********************************************************************
c  *                                                                   *
c  *                    Integer Function ichrd                         *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.01
c  Written by Gordon A. Fenton, Princeton, 1989
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  performs a free format real or integer read from an internal string
c
c  Performs the equivalent of a free format read from an internal record
c  (ie. a character string). Some compilers do not allow you to do this
c  using the READ(STR,*) function directly.
c  This function takes four arguments;
c
c    STR     character string from which to read the data
c
c    ITYPE   integer flag to indicate type of data to be read (all the same)
c            = 1 for integer data
c            = 2 for single precision data
c
c    N       integer giving the expected number of data values to read in
c
c    VAL     vector of length at least N to store the values being read in
c
c  The function returns 0 if the read was successful, `i' if not, where `i'
c  is the index of the value being read at the time of an error.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      integer function ichrd( str, itype, n, val )
      character*(*) str
      character*7 ufmt
      integer itype, n
      real val(*)

      ufmt = '(f  .0)'
      if( itype .eq. 1 ) ufmt = '(i  )  '
      ichrd = 0
      len = lnblnk( str )
      i = 0
      k = 1
  10  i = i + 1
         if( i .gt. len ) go to 40
         if( str(i:i) .eq. ' ' .or. str(i:i) .eq. ',' ) go to 10
         ist = i
  20  i = i + 1
         if( i .gt. len ) go to 30
         if( str(i:i) .ne. ' '  .and. str(i:i) .ne. ',' ) go to 20

  30     ufmt(3:4) = '  '
         write(ufmt(3:4),'(i2)') i - ist
         read(str(ist:i-1),ufmt,err=40,end=40) val(k)
         if( k .eq. n ) return
         k = k + 1
         go to 10

  40  ichrd = k
      return

      end
