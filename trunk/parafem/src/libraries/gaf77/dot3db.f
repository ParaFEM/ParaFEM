c  *********************************************************************
c  *                                                                   *
c  *                  More Dot Product Functions                       *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.0
c  Written by Gordon A. Fenton, TUNS, Dec. 6, 1993
c
c  PURPOSE  another set of specialized dot product functions designed for
c           use by LAS3G.
c
c  This file includes a set of functions which compute dot products for
c  LAS3G of the form
c
c             dot = {A^T}{Z} + {C^T}{U}
c
c  as well as some simple fixed length dot products.
c  In general arguments are as follows;
c
c    Z    real vector containing the parent cell values. Z(1) is assumed
c         to be the first element included in the dot product. (input)
c
c    A    real vector containing the BLUE coefficients. (input)
c
c    C    real vector containing the Covariance coefficients. (input)
c
c    U    real vector containing the random noise inputs. (input)
c
c    dotn external simple dot product function of fixed length. This is
c         used to compute {C^T}{U}.
c
c    j0-j2
c         pointers to cells in left edge of parent cell neighborhood.
c---------------------------------------------------------------------------
c-------------------------------------------------- 4 element corner -------
      real function dot4c( Z, A, C, U, dotn, j0, j1 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot4c = A(1)*Z(j0) + A(2)*Z(1+j0)
     >      + A(3)*Z(j1) + A(4)*Z(1+j1)
     >      + dotn( C, U )

      return
      end
c-------------------------------------------------- 6 element side (hor) ---
      real function dot6h( Z, A, C, U, dotn, j0, j1 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot6h = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(2+j0)
     >      + A( 4)*Z(j1) + A( 5)*Z(1+j1) + A( 6)*Z(2+j1)
     >      + dotn( C, U )

      return
      end
c-------------------------------------------------- 6 element side (vert) --
      real function dot6v( Z, A, C, U, dotn, j0, j1, j2 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot6v = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(j1) + A( 4)*Z(1+j1)
     >      + A( 5)*Z(j2) + A( 6)*Z(1+j2)
     >      + dotn( C, U )

      return
      end
c-------------------------------------------------- 9 element interior -----
      real function dot9i( Z, A, C, U, dotn, j0, j1, j2 )
      dimension Z(*), A(*), C(*), U(*)
      external dotn

      dot9i = A( 1)*Z(j0) + A( 2)*Z(1+j0) + A( 3)*Z(2+j0)
     >      + A( 4)*Z(j1) + A( 5)*Z(1+j1) + A( 6)*Z(2+j1)
     >      + A( 7)*Z(j2) + A( 8)*Z(1+j2) + A( 9)*Z(2+j2)
     >      + dotn( C, U )

      return
      end
