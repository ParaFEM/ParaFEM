c  ***********************************************************************
c  *                                                                     *
c  *                          Subroutine vgath                           *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  create a vector by selecting elements from another vector
c
c  This routine selects elements from the vector {x} and puts them in
c  the vector {z}. In particular, beginning with the first element of
c  {x}, every i'th element is put into {z} until n elements have been
c  extracted;
c
c   {z} = { x(1), x(1+i), x(1+2i), ..., x(1+(n-1)i) }
c
c  Arguments to this routine are as follows;
c
c	x	real vector of length at least (1 + (n-1)*i) containing
c		the superset of elements. (input)
c
c	i	integer containing the `stride' at which x is to be sampled.
c		(input)
c
c	n	integer containing the total number of elements to be
c		extracted from x. (input)
c
c	z	real vector of length at least n into which the n elements
c		extracted from x are stored. (output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      subroutine vgath( x, i, n, z )
      dimension x(*), z(*)

      do 10 j = 1, n
         z(j) = x(1 + (j-1)*i)
  10  continue

      return
      end
