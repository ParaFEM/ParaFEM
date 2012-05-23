c  ************************************************************************
c  *                                                                      *
c  *                         Function pder2                               *
c  *                                                                      *
c  ************************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, Dec. 1990
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  estimates the second partial derivative of a supplied function
c           with respect to (X_i,X_j).
c
c  Estimates numerically the second partial derivative of the supplied function
c  `G' with respect to its (i,j)'th variables. It is assumed that G returns the
c  function evaluated at the location X = {X_1,X_2,...,X_n} and is called
c  by G(X). If i and j differ, the derivative is evaluated using a centered
c  divided difference formulation that involves 7 calls to G(X) and has an
c  accuracy in the order of h**2. If i = j, then another centered difference
c  formulation is used that involves 5 calls to G(X) and has an accuracy in
c  the order of h**4. In this routine, h is fixed at 0.0001
c  Arguments to the routine are as follows;
c
c    G   externally supplied function.
c
c    X   real vector specifying the point at which the derivative is
c        desired. This vector is modified within this routine, but its
c        elements are restored prior to returning. (input)
c
c  i,j   integer indices: the 2nd partial derivative will be evaluated with
c        respect to X_i and X_j (the i and j'th components of X). (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      real function pder2( G, X, i, j )
      dimension X(*)
      external G
c			twohin = 1./(2*h*h),  twlvh2 = 1./(12*h*h)
      data h/0.001/, twoh/0.002/, threeh/0.003/, hhinv/1.e6/
      data a1/2.722222222222222222222/
      data a2/1.5/
      data a3/0.15/
      data a4/0.011111111111111111111/
      data b1/0.375/
      data b2/0.0375/
      data b3/0.00277777777777777777778/

c				i = j we do by interpolating along a line
      if( i .eq. j ) then
c				first store X_i to restore later
         XI = X(i)
c				obtain 7 function points
         T1   = G(X)
         X(i) = XI + h
         T2   = G(X)
         X(i) = XI + twoh
         T3   = G(X)
         X(i) = XI + threeh
         T4   = G(X)
         X(i) = XI - h
         T5   = G(X)
         X(i) = XI - twoh
         T6   = G(X)
         X(i) = XI - threeh
         T7   = G(X)

         X(i) = XI

         pder2 = hhinv*(-a1*T1 + a2*(T2+T5) - a3*(T3+T6) + a4*(T4+T7))

         return
      endif
c				now if i != j, interpolate in 2-D
c				first store X_i and X_j to restore later
      XI = X(i)
      XJ = X(j)
c				obtain the 12 function points
      X(i) = XI + h
      X(j) = XJ + h
      T1   = G(X)
      X(i) = XI - h
      T2   = G(X)
      X(j) = XJ - h
      T3   = G(X)
      X(i) = XI + h
      T4   = G(X)

      X(i) = XI + twoh
      X(j) = XJ + twoh
      T5   = G(X)
      X(i) = XI - twoh
      T6   = G(X)
      X(j) = XJ - twoh
      T7   = G(X)
      X(i) = XI + twoh
      T8   = G(X)

      X(i) = XI + threeh
      X(j) = XJ + threeh
      T9   = G(X)
      X(i) = XI - threeh
      T10  = G(X)
      X(j) = XJ - threeh
      T11  = G(X)
      X(i) = XI + threeh
      T12  = G(X)

c				reset X
      X(i) = XI
      X(j) = XJ
c				estimate the derivative

      pder2 = hhinv*( b1*(T1-T2+T3-T4) - b2*(T5-T6+T7-T8)
     >               + b3*(T9-T10+T11-T12))

      return
      end
