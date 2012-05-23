c  *********************************************************************
c  *                                                                   *
c  *                         subroutine srtgrp                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Fri May 20 11:30:12 1994
c
c  PURPOSE  sorts a sequence of values into ascending order and determines
c           the number of subgroups and assigns pointers to them.
c
c  The sorting portion of this routine is based on Press etal's, Numerical
c  Recipes, heapsort algorithm `indexx', pg 248. In addition to sorting the
c  values using an index vector, this routine determines the number of
c  subgroups (under the assumption that only a limited number of different
c  values are actually present in the input vector) and assigns pointers
c  into the index array to the beginning of each group.
c  Arguments to this routine are as follows;
c
c	x	real array of length at least max(ind(i); i = i0,...,i1)
c		containing the elements to be sorted. (input)
c
c	ind	integer vector of length at least i1 which on output
c		will contain the indices of successively larger elements
c		of x, ie. such that
c		   x(ind(i0)) <= x(ind(i0+1)) <= ... <= x(ind(i1))
c		On first input, ind is assumed to contain ind(i) = i, for
c		i = i0, i0+1, ..., i1. This must be set by the calling
c		routine. (input/output)
c
c	i0	starting index into the vector ind. (input)
c
c	i1	ending index into the vector ind. (input)
c
c	ng	number of distinct groups found. Each x(i) is assumed to
c		belong to a group having value v_j if
c		    abs((x(i)-v_j)/v_j) < tol
c		where tol = 1.e-3 (see below), for j = 1, 2, ..., ng.
c		(output)
c
c	ig	integer vector of length at most (i1 - i0 + 2) which
c		on output will contain the indices into the vector ind
c		of the beginning of each group. ig(ng+1) is set to equal
c		(i1+1) so that the number of elements in group j can
c		be computed as ig(j+1) - ig(j). (output)
c
c	lg	logical flag which is true if the number of groups and
c		their corresponding indices are to be found (ie. if lg
c		is false, then ng and ig are not computed). (input)
c-------------------------------------------------------------------------
      subroutine srtgrp( x, ind, i0, i1, ng, ig, lg )
      real x(*)
      integer ind(*), ig(*)
      logical lg, debug
      common/dbgsrc/ debug
      data zero/0.0/, tol/1.e-3/
c					dump debug data?
   1  format(a,i4,a)
   2  format(8f8.2)
      if( debug ) then
         np = i1 - i0 + 1
         write(6,1)'Sorting ',np,' values of x:'
         write(6,2) (x(i),i=i0,i1)
      endif
c					first sort according to increasing x(i)
      is = i0 - 1
      n  = i1 - is
      L  = 1 + n/2
      ir = n
  10  if( L .gt. 1 ) then
         L  = L - 1
         it = ind(L+is)
         q  = x(it)
      else
         it = ind(ir+is)
         q  = x(it)
         ind(ir+is) = ind(i0)
         ir = ir - 1
         if( ir .eq. 1 ) then
            ind(i0) = it
            if( lg ) go to 30
            if( debug ) write(6,1)'Done (no group count).'
            return
         endif
      endif
      i = L
      j = L + L
  20  if( j .le. ir ) then
         if( j .lt. ir ) then
            if( x(ind(j+is)) .lt. x(ind(j+i0)) ) j = j + 1
         endif
         if( q .lt. x(ind(j+is)) ) then
            ind(i+is) = ind(j+is)
            i = j
            j = j + j
         else
            j = ir + 1
         endif
         go to 20
      endif
      ind(i+is) = it
      go to 10
c					determine number of groups and pointers
  30  i     = i0
      ig(1) = i
      ng    = 1
      v     = x(ind(i))
      do 40 i = i0+1, i1
         e = abs(x(ind(i)) - v)
         if( v .ne. zero ) e = e/abs(v)
         if( e .gt. tol ) then
            ng = ng + 1
            ig(ng) = i
            v = x(ind(i))
         endif
  40  continue
      ig(ng+1) = i1 + 1
      if( debug ) write(6,1)'Done. Number of groups = ',ng

      return
      end
