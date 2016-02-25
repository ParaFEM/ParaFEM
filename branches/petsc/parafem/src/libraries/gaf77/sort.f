c  *********************************************************************
c  *                                                                   *
c  *                        subroutine sort                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Obtained from Numerical Recipes and implemented here by
c  Gordon A. Fenton, TUNS, Jun 11, 1998
c  Latest Update: Jun 11, 1998
c
c  PURPOSE  sorts a sequence of numbers into ascending order
c
c  DESCRIPTION
c
c  ARGUMENTS
c     arr	real vector of length at least n which contains the raw
c		data on input and the sorted (into ascending order) data
c		on output. (input/output)
c
c	n	number of elements in the vector arr. (input)
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      subroutine sort(arr,n)
      parameter (m = 7, nstack = 100)
      dimension arr(*)
      integer istack(nstack)

      jstack = 0
      l = 1
      ir = n
   1  if(ir-l.lt.M)then
        do 12 j = l+1,ir
          a = arr(j)
          do 11 i = j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1) = arr(i)
  11      continue
          i = 0
   2      arr(i+1) = a
  12    continue
        if(jstack.eq.0)return
        ir = istack(jstack)
        l = istack(jstack-1)
        jstack = jstack-2
      else
        k = (l+ir)/2
        temp = arr(k)
        arr(k) = arr(l+1)
        arr(l+1) = temp
        if(arr(l+1).gt.arr(ir))then
          temp = arr(l+1)
          arr(l+1) = arr(ir)
          arr(ir) = temp
        endif
        if(arr(l).gt.arr(ir))then
          temp = arr(l)
          arr(l) = arr(ir)
          arr(ir) = temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp = arr(l+1)
          arr(l+1) = arr(l)
          arr(l) = temp
        endif
        i = l+1
        j = ir
        a = arr(l)
   3    continue
          i = i+1
        if(arr(i).lt.a)goto 3
   4    continue
          j = j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        goto 3
   5    arr(l) = arr(j)
        arr(j) = a
        jstack = jstack+2
c        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
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
