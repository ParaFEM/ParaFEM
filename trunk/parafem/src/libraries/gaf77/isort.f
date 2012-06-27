c  *********************************************************************
c  *                                                                   *
c  *                       subroutine isort                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Obtained from Numerical Recipes and implemented here by
c  Gordon A. Fenton, TUNS, Jun 11, 1998
c  Latest Update: Sep 21, 2001
c
c  PURPOSE  sorts a sequence of integers into ascending order
c
c  DESCRIPTION
c
c  ARGUMENTS
c      iv	integer vector of length at least n which contains the raw
c		data on input and the sorted (into ascending order) data
c		on output. (input/output)
c
c	n	number of elements in the vector arr. (input)
c
c  REVISION HISTORY:
c  1.1	renamed this routine 'isort', as intended (Sep 21/01)
c-------------------------------------------------------------------------
      subroutine isort(iv,n)
      parameter (m = 7, nstack = 100)
      integer iv(*), a, temp
      integer istack(nstack)

      jstack = 0
      l = 1
      ir = n
   1  if(ir-l.lt.M)then
        do 12 j = l+1,ir
          a = iv(j)
          do 11 i = j-1,1,-1
            if(iv(i).le.a)goto 2
            iv(i+1) = iv(i)
  11      continue
          i = 0
   2      iv(i+1) = a
  12    continue
        if(jstack.eq.0)return
        ir = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack-2
      else
        k = (l+ir)/2
        temp = iv(k)
        iv(k) = iv(l+1)
        iv(l+1) = temp
        if(iv(l+1).gt.iv(ir))then
          temp = iv(l+1)
          iv(l+1) = iv(ir)
          iv(ir) = temp
        endif
        if(iv(l).gt.iv(ir))then
          temp = iv(l)
          iv(l) = iv(ir)
          iv(ir) = temp
        endif
        if(iv(l+1).gt.iv(l))then
          temp = iv(l+1)
          iv(l+1) = iv(l)
          iv(l) = temp
        endif
        i = l+1
        j = ir
        a = iv(l)
   3    continue
          i = i+1
        if(iv(i).lt.a)goto 3
   4    continue
          j = j-1
        if(iv(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp = iv(i)
        iv(i) = iv(j)
        iv(j) = temp
        goto 3
   5    iv(l) = iv(j)
        iv(j) = a
        jstack = jstack+2
c        if(jstack.gt.NSTACK)pause 'NSTACK too small in isort'
        if(ir-i+1.ge.j-l)then
          istack(jstack) = ir
          istack(jstack-1) = i
          ir = j-1
        else
          istack(jstack) = j-1
          istack(jstack-1) = l
          l = i
        endif
      endif
      goto 1
      END
