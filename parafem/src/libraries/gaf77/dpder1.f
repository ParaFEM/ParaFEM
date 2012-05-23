c  ************************************************************************
c  *                                                                      *
c  *                        Function dpder1                               *
c  *                                                                      *
c  ************************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, Dec. 1990
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  estimates the first partial derivative of a supplied function
c
c  Estimates numerically the first partial derivative of the supplied function
c  `G' with respect to its i'th variable. It is assumed that G returns the
c  function evaluated at the location X = {X_1,X_2,...,X_n} and is called
c  by G(X). The derivative is evaluated using a centered divided difference
c  formulation that involves 6 calls to G(X) and has an accuracy in the
c  order of h**6. In this routine, h is fixed at 0.001
c  Arguments to the routine are as follows;
c
c    G   externally supplied function.
c
c    X   real vector specifying the point at which the derivative is
c        desired. This vector is modified within this routine, but its
c        elements are restored prior to returning. (input)
c
c    i   integer index: the partial derivative will be evaluated with
c        respect to X_i (the i'th component of X). (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      real*8 function dpder1( G, X, i )
      implicit real*8 (a-h,o-z)
      dimension X(*)
      external G
      data h/0.001d0/, twoh/0.002d0/, threeh/0.003d0/, hinv/1000.d0/
      data a1/0.01666666666666666666667d0/
      data a2/0.15d0/
      data a3/0.75d0/

c				first store X_i to restore later
      T = X(i)
c				obtain 6 function values
      X(i) = T + threeh
      T1   = G(X)
      X(i) = T + twoh
      T2   = G(X)
      X(i) = T + h
      T3   = G(X)
      X(i) = T - h
      T4   = G(X)
      X(i) = T -twoh
      T5   = G(X)
      X(i) = T - threeh
      T6   = G(X)
c				reset X
      X(i) = T
c				estimate the derivative

      dpder1 = hinv*(a1*(T1-T6) - a2*(T2-T5) + a3*(T3-T4))

      return
      end
