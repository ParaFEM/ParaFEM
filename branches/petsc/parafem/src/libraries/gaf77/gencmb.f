c  *********************************************************************
c  *                                                                   *
c  *                        subroutine gencmb                          *
c  *                                                                   *
c  *********************************************************************
c  Integer Version 1.1
c  Written by Gordon A. Fenton, TUNS, Wed May 25 16:52:41 1994
c  Latest Update: May 26, 1994
c
c  PURPOSE  generates sequential elements in a combination
c
c  This routine generates, one at a time, the sequence of combinations found
c  when m indistinguishable objects are placed in n slots (m <= n). The number
c  of such combinations is choose(n,m) = n!/m!(n-m)!. On each call to this
c  routine, the next combination is returned (if this routine is called more
c  than choose(n,m) times, the last combination is repeatedly returned).
c  Each individual combination is defined by n binary digits (0 or ib)
c  containing m ib's in the appropriate locations. For example, if the problem
c  is 5 choose 2, then the first 5 (out of 10) combinations returned are
c
c	1 1 0 0 0
c	1 0 1 0 0
c	1 0 0 1 0
c	1 0 0 0 1
c	0 1 1 0 0
c
c  if ib = 1, etc. Arguments to this routine are as follows;
c
c    icmb   integer vector of length at least n containing the next combination
c           on output. It should contain the current combination on input if
c           this is not the first call to gencmb. (input/output)
c
c    n      number of `slots' into which the m objects are to be placed.
c           If n changes on any call after the first, the new sequence is
c           started and the old abandoned. (input)
c
c    m      number of indistinguishable objects being placed. (input)
c           If m changes on any call after the first, the new sequence is
c           started and the old abandoned. (input)
c
c    ib     integer to which the object placement should be set. If ib is
c           0, then object placements are denoted by 1. On output ib is the
c           value actually used. (input/output)
c
c  PARAMETERS
c    MXM    maximum value of m allowed (currently 100)
c-------------------------------------------------------------------------
      subroutine gencmb( icmb, n, m, ib )
      parameter (MXM = 100)
      integer icmb(*), k(MXM)
      save k
      data nn/0/, mm/0/
c						initial combination?
      if( (n .ne. nn) .or. (m .ne. mm) ) then
         nn = n
         mm = m
         do 10 i = 1, m
            icmb(i) = ib
            k(i)    = i
  10     continue
         do 20 i = m+1, n
            icmb(i) = 0
  20     continue
         return
      endif
      if( ib .eq. 0 ) ib = 1
      nl = n
      ml = m
  30  if( k(ml) .eq. nl ) then
         if( ml .eq. 1 ) return
         if( k(ml-1) .eq. (nl-1) ) then
            ml = ml - 1
            nl = nl - 1
            go to 30
         endif
         icmb(k(ml-1)) = 0
         k(ml-1) = k(ml-1) + 1
         icmb(k(ml-1)) = ib
         do 40 i = ml, m
            icmb(k(i)) = 0
            k(i) = k(i-1) + 1
            icmb(k(i)) = ib
  40     continue
      else
         icmb(k(ml)) = 0
         k(ml) = k(ml) + 1
         icmb(k(ml)) = ib
      endif

      return
      end

