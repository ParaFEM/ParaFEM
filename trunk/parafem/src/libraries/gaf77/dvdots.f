c  *********************************************************************
c  *                                                                   *
c  *                 Scalar Dot Product Functions                      *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, 1989.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  a set of dot product functions
c
c  This is a set of dot product functions written specifically for LAS3D
c  but with application anywhere. The functions are specialized to fixed
c  vector lengths to make them run as fast as possible, however, a
c  general dot function `DDOTN' has been added to the end to handle
c  vectors of any length N.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-------------------------------------------------------------------------
c                                       vectors of length 2
      real*8 function ddot2(e,f)
      real*8 e(*), f(*)
      ddot2 = e(1)*f(1) + e(2)*f(2)
      return
      end
c                                       vectors of length 3
      real*8 function ddot3(e,f)
      real*8 e(*), f(*)
      ddot3 = e(1)*f(1) + e(2)*f(2) + e(3)*f(3)
      return
      end
c                                       vectors of length 4
      real*8 function ddot4(e,f)
      real*8 e(*), f(*)
      ddot4 = e(1)*f(1) + e(2)*f(2) + e(3)*f(3) + e(4)*f(4)
      return
      end
c                                       vectors of length 5
      real*8 function ddot5(e,f)
      real*8 e(*), f(*)
      ddot5 = e(1)*f(1) + e(2)*f(2) + e(3)*f(3) + e(4)*f(4)
     >      + e(5)*f(5)
      return
      end
c                                       vectors of length 6
      real*8 function ddot6(e,f)
      real*8 e(*), f(*)
      ddot6 = e(1)*f(1) + e(2)*f(2) + e(3)*f(3) + e(4)*f(4)
     >      + e(5)*f(5) + e(6)*f(6)
      return
      end
c                                       vectors of length 7
      real*8 function ddot7(e,f)
      real*8 e(*), f(*)
      ddot7 = e(1)*f(1) + e(2)*f(2) + e(3)*f(3) + e(4)*f(4)
     >      + e(5)*f(5) + e(6)*f(6) + e(7)*f(7)
      return
      end
c                                       vectors of length 8
      real*8 function ddot8(e,f)
      real*8 e(*), f(*)
      ddot8 = e(1)*f(1) + e(2)*f(2) + e(3)*f(3) + e(4)*f(4)
     >      + e(5)*f(5) + e(6)*f(6) + e(7)*f(7) + e(8)*f(8)
      return
      end
c                                       vectors of length 12
      real*8 function ddot12(e,f)
      real*8 e(*), f(*)
      ddot12 = e( 1)*f( 1) + e( 2)*f( 2) + e( 3)*f( 3) + e( 4)*f( 4)
     >       + e( 5)*f( 5) + e( 6)*f( 6) + e( 7)*f( 7) + e( 8)*f( 8)
     >       + e( 9)*f( 9) + e(10)*f(10) + e(11)*f(11) + e(12)*f(12)
      return
      end
c                                       vectors of length 18
      real*8 function ddot18(e,f)
      real*8 e(*), f(*)
      ddot18 = e( 1)*f( 1) + e( 2)*f( 2) + e( 3)*f( 3) + e( 4)*f( 4)
     >       + e( 5)*f( 5) + e( 6)*f( 6) + e( 7)*f( 7) + e( 8)*f( 8)
     >       + e( 9)*f( 9) + e(10)*f(10) + e(11)*f(11) + e(12)*f(12)
     >       + e(13)*f(13) + e(14)*f(14) + e(15)*f(15) + e(16)*f(16)
     >       + e(17)*f(17) + e(18)*f(18)
      return
      end
c                                       vectors of length 27
      real*8 function ddot27(e,f)
      real*8 e(*), f(*)
      ddot27 = e( 1)*f( 1) + e( 2)*f( 2) + e( 3)*f( 3) + e( 4)*f( 4)
     >       + e( 5)*f( 5) + e( 6)*f( 6) + e( 7)*f( 7) + e( 8)*f( 8)
     >       + e( 9)*f( 9) + e(10)*f(10) + e(11)*f(11) + e(12)*f(12)
     >       + e(13)*f(13) + e(14)*f(14) + e(15)*f(15) + e(16)*f(16)
     >       + e(17)*f(17) + e(18)*f(18) + e(19)*f(19) + e(20)*f(20)
     >       + e(21)*f(21) + e(22)*f(22) + e(23)*f(23) + e(24)*f(24)
     >       + e(25)*f(25) + e(26)*f(26) + e(27)*f(27)
      return
      end
c                                       vectors of length N
      real*8 function ddotn(e,f,n)
      real*8 e(*), f(*)
      ddotn = e(1)*f(1)
      do 10 i = 2, n
         ddotn = ddotn + e(i)*f(i)
  10  continue
      return
      end
