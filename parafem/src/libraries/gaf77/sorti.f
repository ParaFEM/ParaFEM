c  *********************************************************************
c  *                                                                   *
c  *                       subroutine sorti                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Obtained from Numerical Recipes and implemented here by
c  Gordon A. Fenton, DalTech, Mar 3, 2000
c  Latest Update: Mar 3, 2000
c
c  PURPOSE  sorts a sequence of real numbers into ascending order by ordering
c           an index vector
c
c  DESCRIPTION
c
c  ARGUMENTS
c
c     arr	real vector of length at least n which contains the raw
c		data to be sorted. (input)
c
c    indx       integer vector of length at least n which on output will
c		contain the indices into the vector `arr' so that arr(indx(1))
c		is the smallest element of arr, arr(indx(2)) is the second
c		smallest, and so on. Specifically,
c
c		arr(indx(1)) < arr(indx(2)) < ... < arr(indx(n))
c
c               (output)
c
c	n	number of elements in the vector arr. (input)
c
c  REVISION HISTORY:
c
c-------------------------------------------------------------------------
      SUBROUTINE sorti(arr,indx,n)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=200)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sorti'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software !00,]v45.3.
