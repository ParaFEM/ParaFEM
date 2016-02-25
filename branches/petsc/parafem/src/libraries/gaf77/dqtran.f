c  *********************************************************************
c  *                                                                   *
c  *                         Subroutine dqtran                         *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.11
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to construct a matrix [Q] which rotates a vector {x} into {y}
c
c  This routine accepts two vectors {x} and {y} and constructs a matrix [Q]
c  such that [Q]{x} = {y}. Note that [Q] will be orthogonal if the lengths
c  of {x} and {y} are the same (based on the 2-norm). The matrix is formed
c  by noting that
c
c      [U]{x} = -sign(x_1)*|| x ||*e_1
c
c      [V]{y} = -sign(y-1)*|| y ||*e_1
c
c  where [U] and [V] are orthogonal Householder transformation matrices.
c  From these, [Q] can be determined to be
c
c      [Q] = r*[V-transpose][U]
c
c             sign(y_1)*|| y ||
c  where r = -------------------
c             sign(x_1)*|| x ||
c
c  is essentially the ratio in lengths of {y} and {x}. Arguments to the routine
c  are as follows;
c
c     Q     real array of size at least n x n which on output will contain
c           the elements of the transformation matrix [Q]. (output)
c
c    iq     leading dimension of the array Q exactly as specified in the
c           calling routine. (input)
c
c     x     real vector of length at least n which contains the elements
c           of the vector {x}. (input)
c
c     y     real vector of length at least n which contains the elements
c           of the vector {y}. (input)
c
c     n     length of the vectors {x} and {y}. (input)
c
c  REVISION HISTORY:
c  1.1	eliminated unused local variable `zero' (Dec 5/96)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine dqtran( Q, iq, x, y, n )
      implicit real*8 (a-h,o-z)
      dimension Q(iq,*), x(*), y(*)
      data one/1.d0/
      sign(u,v) = dsign(u,v)
      abs(u)    = dabs(u)
      sqrt(u)   = dsqrt(u)
c					save x(1) and y(1) to restore later
      x1 = x(1)
      y1 = y(1)
c					compute 2-norms
      xn = x1*x1
      yn = y1*y1
      yx = y1*x1
      do 10 i = 2, n
         xn = xn + x(i)*x(i)
         yn = yn + y(i)*y(i)
         yx = yx + y(i)*x(i)
  10  continue
c					compute coefficients
      xn  = sqrt(xn)
      yn  = sqrt(yn)
      sxn = sign( xn, x1 )
      syn = sign( yn, y1 )
      r   = syn/sxn
      a   = r/(xn*(xn + abs(x1)))
      b   = one/(yn*(yn + abs(y1)))
      c   = a*b*(yx + y1*sxn + x1*syn + sxn*syn)
      b   = r*b

c					define Householder vectors
      x(1) = x1 + sxn
      y(1) = y1 + syn
c					set up transformation
      do 20 i = 1, n
         Q(i,i) = r - a*x(i)*x(i) - b*y(i)*y(i) + c*y(i)*x(i)
         do 20 j = i+1, n
            aa = a*x(i)*x(j) + b*y(i)*y(j)
            Q(i,j) = c*y(i)*x(j) - aa
            Q(j,i) = c*y(j)*x(i) - aa
  20  continue
c					restore {x} and {y}
      x(1) = x1
      y(1) = y1
c					all done
      return
      end
