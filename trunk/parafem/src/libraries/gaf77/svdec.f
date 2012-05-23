c  *********************************************************************
c  *                                                                   *
c  *                        Subroutine svdec                           *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.01
c  From Numerical Recipes, Press etal., 1988
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the singular value decomposition [A] = [U][W][V'] of
c           an m x n matrix [A] (m >= n).
c
c  Given an m x n matrix [A], this routine computes its singular value
c  decomposition [A] = [U][W][V'] (where prime denotes transpose). The
c  resulting matrices [U], [W], and [V] are of size m x n, n x n diagonal,
c  and n x n, respectively. [W] is stored as a vector of length n. The
c  dimension m must be greater than or equal to n. If it isn't, then
c  the matrix A should be filled to square with zero rows.
c  Arguments to the routine are as follows;
c
c     A    real array of size at least m x n (m >= n) which on input
c          contains the matrix to be factorized. On output A contains
c          the elements of the matrix [U]. (input/output)
c
c    ia    integer giving the leading dimension of A as specified in the
c          calling routine. (input)
c
c     W    real vector of length at least n which on output will contain
c          the diagonal matrix of singular values. (output)
c
c     V    real array of size at least n x n which on output will contain
c          the elements of the matrix [V] (not the transpose). (output)
c
c    iv    integer giving the leading dimension of V as specified in the
c          calling routine. (input)
c
c   rv1    temporary vector of length at least n.
c
c     m    number of rows in [A]. `m' must be greater than or equal to `n'
c          (which can be accomplished by filling out A with zero rows if
c          necessary). (input)
c
c     n    number of columns in [A]. (input)
c
c  ierr    integer error flag which takes the following values;
c          =  0 if all goes well
c          = -1 if m < n
c          = -2 if convergence not achieved in MAXIT iterations (see param's)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine svdec( A, ia, W, V, iv, rv1, m, n, ierr )
      parameter (MAXIT = 30)
      dimension A(ia,*), W(*), V(iv,*), rv1(*)
      data zero/0.0/, one/1.0/, two/2.0/
c					check m >= n
      if( m .lt. n ) then
         ierr = -1
         return
      endif
c					initialize
      g     = zero
      scale = zero
      anorm = zero
c					Householder reduction to bidiagonal
      do 140 i = 1, n
         l      = i + 1
         rv1(i) = scale*g
         g      = zero
         s      = zero
         scale  = zero
         if( i .le. m ) then
            do 10 k = i, m
               scale = scale + abs(a(k,i))
  10        continue
            if( scale .ne. zero ) then
               do 20 k = i, m
                  a(k,i) = a(k,i)/scale
                  s      = s + a(k,i)*a(k,i)
  20           continue
               f      = a(i,i)
               g      = -sign(sqrt(s),f)
               h      = f*g - s
               a(i,i) = f - g
               if( i .ne. n ) then
                  do 50 j = l, n
                     s = zero
                     do 30 k = i, m
                        s = s + a(k,i)*a(k,j)
  30                 continue
                     f = s/h
                     do 40 k = i, m
                        a(k,j) = a(k,j) + f*a(k,i)
  40                 continue
  50              continue
               endif
               do 60 k = i, m
                  a(k,i) = scale*a(k,i)
  60           continue
            endif
         endif
         w(i)  = scale*g
         g     = zero
         s     = zero
         scale = zero
         if( (i.le.m) .and. (i.ne.n) ) then
            do 70 k = l, n
               scale = scale + abs(a(i,k))
  70        continue
            if( scale .ne. zero ) then
               do 80 k = l, n
                  a(i,k) = a(i,k)/scale
                  s      = s + a(i,k)*a(i,k)
  80           continue
               f      = a(i,l)
               g      = -sign(sqrt(s),f)
               h      = f*g - s
               a(i,l) = f - g
               do 90 k = l, n
                  rv1(k) = a(i,k)/h
  90           continue
               if( i .ne. m ) then
                  do 120 j = l, m
                     s = zero
                     do 100 k = l, n
                        s = s + a(j,k)*a(i,k)
 100                 continue
                     do 110 k = l, n
                        a(j,k) = a(j,k) + s*rv1(k)
 110                 continue
 120              continue
               endif
               do 130 k = l, n
                  a(i,k) = scale*a(i,k)
 130           continue
            endif
         endif
         anorm = amax1( anorm, (abs(w(i))+abs(rv1(i))) )
 140  continue
c					accumulate right-hand [V] transforms
      do 200 i = n, 1, -1
         if( i .lt. n ) then
            if( g .ne. zero ) then
               do 150 j = l, n
                  v(j,i) = (a(i,j)/a(i,l))/g
 150           continue
               do 180 j = l, n
                  s = zero
                  do 160 k = l, n
                     s = s + a(i,k)*v(k,j)
 160              continue
                  do 170 k = l, n
                     v(k,j) = v(k,j) + s*v(k,i)
 170              continue
 180           continue
            endif
            do 190 j = l, n
               v(i,j) = zero
               v(j,i) = zero
 190        continue
         endif
         v(i,i) = one
         g      = rv1(i)
         l      = i
 200  continue
c					accumulate left-hand [U] transforms
      do 270 i = n, 1, -1
         l = i+1
         g = w(i)
         if( i .lt. n ) then
            do 210 j = l, n
               a(i,j) = zero
 210        continue
         endif
         if( g .ne. zero ) then
            g = one/g
            if( i .ne. n ) then
               do 240 j = l, n
                  s = zero
                  do 220 k = l, m
                     s = s + a(k,i)*a(k,j)
 220              continue
                  f = s*g/a(i,i)
                  do 230 k = i, m
                     a(k,j) = a(k,j) + f*a(k,i)
 230              continue
 240           continue
            endif
            do 250 j = i, m
               a(j,i) = a(j,i)*g
 250        continue
         else
            do 260 j = i, m
               a(j,i) = zero
 260        continue
         endif
         a(i,i) = a(i,i) + one
 270  continue
c					diagonalize the bidiagonal form
      do 380 k = n, 1, -1
c						iterate MAXIT times
         do 370 its = 1, MAXIT
c						find `small' values
            do 280 l = k, 1, -1
               nm = l - 1
               if( (abs(rv1(l)) + anorm) .eq. anorm ) go to 320
               if( (abs(w(nm))  + anorm) .eq. anorm ) go to 290
 280        continue
 290        c = zero
            s = one
            do 310 i = l, k
               f = s*rv1(i)
               if( (abs(f)+anorm) .ne. anorm ) then
                  g    =  w(i)
                  h    =  sqrt(f*f+g*g)
                  w(i) =  h
                  h    =  one/h
                  c    =  g*h
                  s    = -f*h
                  do 300 j = 1, m
                     y       =  a(j,nm)
                     a(j,nm) =  c*y + s*a(j,i)
                     a(j,i)  = -s*y + c*a(j,i)
 300              continue
               endif
 310        continue
 320        z = w(k)
            if( l .eq. k ) then
               if( z .lt. zero ) then
                  w(k) = -z
                  do 330 j = 1, n
                     v(j,k) = -v(j,k)
 330              continue
               endif
c						done, find next diagonal term
               go to 380
            endif
c						set up for next iteration
            x  = w(l)
            nm = k - 1
            y  = w(nm)
            g  = rv1(nm)
            h  = rv1(k)
            f  = ((y-z)*(y+z)+(g-h)*(g+h))/(two*h*y)
            g  = sqrt(f*f+one)
            f  = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c  = one
            s  = one
            do 360 j = l, nm
               i      = j + 1
               g      = rv1(i)
               y      = w(i)
               h      = s*g
               g      = c*g
               z      = sqrt(f*f+h*h)
               rv1(j) = z
               c      = f/z
               s      = h/z
               f      = x*c + g*s
               g      =-x*s + g*c
               h      = y*s
               y      = y*c
               do 340 nm = 1, n
                  x       =  v(nm,j)
                  v(nm,j) =  c*x + s*v(nm,i)
                  v(nm,i) = -s*x + c*v(nm,i)
 340           continue
               z    = sqrt(f*f+h*h)
               w(j) = z
               if( z .ne. zero ) then
                  z = one/z
                  c = f*z
                  s = h*z
               endif
               f = c*g + s*y
               x = c*y - s*g
               do 350 nm = 1, m
                  y       =  a(nm,j)
                  a(nm,j) =  c*y + s*a(nm,i)
                  a(nm,i) = -s*y + c*a(nm,i)
 350           continue
 360        continue
            rv1(l) = zero
            rv1(k) = f
            w(k)   = x
 370     continue
c					if we get here, convergence failed
         ierr = -2
         return
 380  continue
c					success!
      ierr = 0
      return
      end
