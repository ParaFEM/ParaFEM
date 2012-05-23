c  ********************************************************************
c  *                                                                  *
c  *                        subroutine plan3d                         *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, Dec. 6, 1993
c  Latest Update: Dec 5, 1996
c
c  PURPOSE   produces stage 1 of the 3-D LAS subdivision in the event that
c            one or two (but not three) of k1, k2, or k3 is equal to 1.
c
c
c  REVISION HISTORY:
c  1.2	explicitly declared dot function types (Dec 5/96)
c--------------------------------------------------------------------------
      subroutine plan3d(Z,iq,jq,k1,k2,k3,AC,AS,AI,CC,CS,CI,U,iout)
      real Z(*), U(*)
      real AC(4,7,*), AS(6,7,*), AI(9,7)
      real CC(28,*), CS(28,*), CI(*)
      logical lk1, lk2, lk3
      real     dot1, dot2, dot3, dot4, dot5, dot6, dot7
      external dot1, dot2, dot3, dot4, dot5, dot6, dot7
      data eight/8.0/

   1  format(a)

      lk1  = (k1 .eq. 1)
      lk2  = (k2 .eq. 1)
      lk3  = (k3 .eq. 1)

      if( .not. lk1 .and. .not. lk2 .and. .not. lk3 ) return

      ix  = 2*k1
      iy  = 2*k2
      ixy = ix*iy
      j0  = jq + 1
      i0  = iq + 1
      i1  = i0 + ix
      i2  = i0 + ixy
      i3  = i2 + ix
c							1-D case

      if( (lk1.and.lk2) .or. (lk1.and.lk3) .or. (lk2.and.lk3) ) then
         jm = max0(k1,k2,k3)
         if( .not. lk1 ) then
            iz = 2
         elseif( .not. lk2 ) then
            iz = 4
         elseif( .not. lk3 ) then
            iz = 8
         else
            write(iout,1)'Can''t subdivide a single cell in plan3d'
            write(iout,1)'(ie. can''t have k1 = k2 = k3 = 1)'
            write(iout,1)'This should never happen!'
            stop
         endif
         call vnorm( U, 7*jm )

         Z(i0)   = AS(1,1,1)*Z(j0) + AS(2,1,1)*Z(j0+1)
     >             + CS(1,1)*U(1)
         Z(i0+1) = AS(1,2,1)*Z(j0) + AS(2,2,1)*Z(j0+1)
     >             + dot2(CS( 2,1),U(1))
         Z(i1)   = AS(1,3,1)*Z(j0) + AS(2,3,1)*Z(j0+1)
     >             + dot3(CS( 4,1),U(1))
         Z(i1+1) = AS(1,4,1)*Z(j0) + AS(2,4,1)*Z(j0+1)
     >             + dot4(CS( 7,1),U(1))
         Z(i2)   = AS(1,5,1)*Z(j0) + AS(2,5,1)*Z(j0+1)
     >             + dot5(CS(11,1),U(1))
         Z(i2+1) = AS(1,6,1)*Z(j0) + AS(2,6,1)*Z(j0+1)
     >             + dot6(CS(16,1),U(1))
         Z(i3)   = AS(1,7,1)*Z(j0) + AS(2,7,1)*Z(j0+1)
     >             + dot7(CS(22,1),U(1))
         Z(i1+1) = eight*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
     >               - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
         L = 8
         do 10 js = 1, jm-2
            i0 = i0 + iz
            i1 = i1 + iz
            i2 = i2 + iz
            i3 = i3 + iz
            Z(i0)   = AI(1,1)*Z(j0  )+AI(2,1)*Z(j0+1)+AI(3,1)*Z(j0+2)
     >              + CI(1)*U(L)
            Z(i0+1) = AI(1,2)*Z(j0  )+AI(2,2)*Z(j0+1)+AI(3,2)*Z(j0+2)
     >              + dot2(CI(2),U(L))
            Z(i1)   = AI(1,3)*Z(j0  )+AI(2,3)*Z(j0+1)+AI(3,3)*Z(j0+2)
     >              + dot3(CI(4),U(L))
            Z(i1+1) = AI(1,4)*Z(j0  )+AI(2,4)*Z(j0+1)+AI(3,4)*Z(j0+2)
     >              + dot4(CI(7),U(L))
            Z(i2)   = AI(1,5)*Z(j0  )+AI(2,5)*Z(j0+1)+AI(3,5)*Z(j0+2)
     >              + dot5(CI(11),U(L))
            Z(i2+1) = AI(1,6)*Z(j0  )+AI(2,6)*Z(j0+1)+AI(3,6)*Z(j0+2)
     >              + dot6(CI(16),U(L))
            Z(i3)   = AI(1,7)*Z(j0  )+AI(2,7)*Z(j0+1)+AI(3,7)*Z(j0+2)
     >              + dot7(CI(22),U(L))
            Z(i1+1) = eight*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
            j0 = j0 + 1
            L  = L  + 7
  10     continue
         i0 = i0 + iz
         i1 = i1 + iz
         i2 = i2 + iz
         i3 = i3 + iz
         Z(i0)   = AS(2,1,1)*Z(j0) + AS(1,1,1)*Z(j0+1)
     >             + CS(1,1)*U(1)
         Z(i0+1) = AS(2,2,1)*Z(j0) + AS(1,2,1)*Z(j0+1)
     >             + dot2(CS( 2,1),U(1))
         Z(i1)   = AS(2,3,1)*Z(j0) + AS(1,3,1)*Z(j0+1)
     >             + dot3(CS( 4,1),U(1))
         Z(i1+1) = AS(2,4,1)*Z(j0) + AS(1,4,1)*Z(j0+1)
     >             + dot4(CS( 7,1),U(1))
         Z(i2)   = AS(2,5,1)*Z(j0) + AS(1,5,1)*Z(j0+1)
     >             + dot5(CS(11,1),U(1))
         Z(i2+1) = AS(2,6,1)*Z(j0) + AS(1,6,1)*Z(j0+1)
     >             + dot6(CS(16,1),U(1))
         Z(i3)   = AS(2,7,1)*Z(j0) + AS(1,7,1)*Z(j0+1)
     >             + dot7(CS(22,1),U(1))
         Z(i1+1) = eight*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
     >               - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
         return
      endif
c							2-D case
      jx = k1
      iz = 2
      if( lk1 ) then
         jx = k2
         iz = 4
      endif
      jy = k3
      if( lk3 ) jy = k2
      j1 = j0 + jx

      call vnorm( U, 7*jx )
c								corner #1
      Z(i0)   = dot4c(Z,AC(1,1,1),CC( 1,1),U(1),dot1,j0,j1)
      Z(i0+1) = dot4c(Z,AC(1,2,1),CC( 2,1),U(1),dot2,j0,j1)
      Z(i1)   = dot4c(Z,AC(1,3,1),CC( 4,1),U(1),dot3,j0,j1)
      Z(i1+1) = dot4c(Z,AC(1,4,1),CC( 7,1),U(1),dot4,j0,j1)
      Z(i2)   = dot4c(Z,AC(1,5,1),CC(11,1),U(1),dot5,j0,j1)
      Z(i2+1) = dot4c(Z,AC(1,6,1),CC(16,1),U(1),dot6,j0,j1)
      Z(i3)   = dot4c(Z,AC(1,7,1),CC(22,1),U(1),dot7,j0,j1)
      Z(i3+1) = eight*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
     >            - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								side #1
      L = 8
      do 40 js = 1, jx-2
         i0 = i0 + iz
         i1 = i1 + iz
         i2 = i2 + iz
         i3 = i3 + iz
         Z(i0)   = dot6h(Z,AS(1,1,1),CS( 1,1),U(L),dot1,j0,j1)
         Z(i0+1) = dot6h(Z,AS(1,2,1),CS( 2,1),U(L),dot2,j0,j1)
         Z(i1)   = dot6h(Z,AS(1,3,1),CS( 4,1),U(L),dot3,j0,j1)
         Z(i1+1) = dot6h(Z,AS(1,4,1),CS( 7,1),U(L),dot4,j0,j1)
         Z(i2)   = dot6h(Z,AS(1,5,1),CS(11,1),U(L),dot5,j0,j1)
         Z(i2+1) = dot6h(Z,AS(1,6,1),CS(16,1),U(L),dot6,j0,j1)
         Z(i3)   = dot6h(Z,AS(1,7,1),CS(22,1),U(L),dot7,j0,j1)
         Z(i3+1) = eight*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                 - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
         j0 = j0 + 1
         j1 = j1 + 1
         L  = L  + 7
  40  continue
c								corner #2
      i0 = i0 + iz
      i1 = i1 + iz
      i2 = i2 + iz
      i3 = i3 + iz
      Z(i0)   = dot4c(Z,AC(1,1,2),CC( 1,2),U(L),dot1,j0,j1)
      Z(i0+1) = dot4c(Z,AC(1,2,2),CC( 2,2),U(L),dot2,j0,j1)
      Z(i1)   = dot4c(Z,AC(1,3,2),CC( 4,2),U(L),dot3,j0,j1)
      Z(i1+1) = dot4c(Z,AC(1,4,2),CC( 7,2),U(L),dot4,j0,j1)
      Z(i2)   = dot4c(Z,AC(1,5,2),CC(11,2),U(L),dot5,j0,j1)
      Z(i2+1) = dot4c(Z,AC(1,6,2),CC(16,2),U(L),dot6,j0,j1)
      Z(i3)   = dot4c(Z,AC(1,7,2),CC(22,2),U(L),dot7,j0,j1)
      Z(i3+1) = eight*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
     >              - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

      j0 = jq + 1
      do 60 ks = 1, jy-2
         j1 = j0 + jx
         j2 = j1 + jx
         i0 = i1 + 2
         if( .not. lk3 ) i0 = i0 + 2*jx
         i1 = i0 + ix
         i2 = i0 + ixy
         i3 = i2 + ix
         call vnorm( U, 7*jx )
c								side #2
         Z(i0)   = dot6v(Z,AS(1,1,2),CS( 1,2),U(1),dot1,j0,j1,j2)
         Z(i0+1) = dot6v(Z,AS(1,2,2),CS( 2,2),U(1),dot2,j0,j1,j2)
         Z(i1)   = dot6v(Z,AS(1,3,2),CS( 4,2),U(1),dot3,j0,j1,j2)
         Z(i1+1) = dot6v(Z,AS(1,4,2),CS( 7,2),U(1),dot4,j0,j1,j2)
         Z(i2)   = dot6v(Z,AS(1,5,2),CS(11,2),U(1),dot5,j0,j1,j2)
         Z(i2+1) = dot6v(Z,AS(1,6,2),CS(16,2),U(1),dot6,j0,j1,j2)
         Z(i3)   = dot6v(Z,AS(1,7,2),CS(22,2),U(1),dot7,j0,j1,j2)
         Z(i3+1) = eight*Z(j1) - Z(i0) - Z(i0+1) - Z(i1)
     >               - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								interior
         L = 8
         do 50 js = 1, jx-2
            i0 = i0 + iz
            i1 = i1 + iz
            i2 = i2 + iz
            i3 = i3 + iz
            Z(i0)   = dot9i(Z,AI(1,1),CI( 1),U(L),dot1,j0,j1,j2)
            Z(i0+1) = dot9i(Z,AI(1,2),CI( 2),U(L),dot2,j0,j1,j2)
            Z(i1)   = dot9i(Z,AI(1,3),CI( 4),U(L),dot3,j0,j1,j2)
            Z(i1+1) = dot9i(Z,AI(1,4),CI( 7),U(L),dot4,j0,j1,j2)
            Z(i2)   = dot9i(Z,AI(1,5),CI(11),U(L),dot5,j0,j1,j2)
            Z(i2+1) = dot9i(Z,AI(1,6),CI(16),U(L),dot6,j0,j1,j2)
            Z(i3)   = dot9i(Z,AI(1,7),CI(22),U(L),dot7,j0,j1,j2)
            Z(i3+1) = eight*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
            j0 = j0 + 1
            j1 = j1 + 1
            j2 = j2 + 1
            L  = L  + 7
  50     continue
c								side #3
         i0 = i0 + iz
         i1 = i1 + iz
         i2 = i2 + iz
         i3 = i3 + iz
         Z(i0)   = dot6v(Z,AS(1,1,3),CS( 1,3),U(L),dot1,j0,j1,j2)
         Z(i0+1) = dot6v(Z,AS(1,2,3),CS( 2,3),U(L),dot2,j0,j1,j2)
         Z(i1)   = dot6v(Z,AS(1,3,3),CS( 4,3),U(L),dot3,j0,j1,j2)
         Z(i1+1) = dot6v(Z,AS(1,4,3),CS( 7,3),U(L),dot4,j0,j1,j2)
         Z(i2)   = dot6v(Z,AS(1,5,3),CS(11,3),U(L),dot5,j0,j1,j2)
         Z(i2+1) = dot6v(Z,AS(1,6,3),CS(16,3),U(L),dot6,j0,j1,j2)
         Z(i3)   = dot6v(Z,AS(1,7,3),CS(22,3),U(L),dot7,j0,j1,j2)
         Z(i3+1) = eight*Z(j1) - Z(i0) - Z(i0+1) - Z(i1)
     >               - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
         j0 = j0 + 2
  60  continue

      j1 = j0 + jx
      i0 = i1 + 2
      if( .not. lk3 ) i0 =i0 + 2*jx
      i1 = i0 + ix
      i2 = i0 + ixy
      i3 = i2 + ix
      call vnorm( U, 7*jx )
c								corner #3
      Z(i0)   = dot4c(Z,AC(1,1,3),CC( 1,3),U(1),dot1,j0,j1)
      Z(i0+1) = dot4c(Z,AC(1,2,3),CC( 2,3),U(1),dot2,j0,j1)
      Z(i1)   = dot4c(Z,AC(1,3,3),CC( 4,3),U(1),dot3,j0,j1)
      Z(i1+1) = dot4c(Z,AC(1,4,3),CC( 7,3),U(1),dot4,j0,j1)
      Z(i2)   = dot4c(Z,AC(1,5,3),CC(11,3),U(1),dot5,j0,j1)
      Z(i2+1) = dot4c(Z,AC(1,6,3),CC(16,3),U(1),dot6,j0,j1)
      Z(i3)   = dot4c(Z,AC(1,7,3),CC(22,3),U(1),dot7,j0,j1)
      Z(i3+1) = eight*Z(j1) - Z(i0) - Z(i0+1) - Z(i1)
     >            - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								side #4
      L = 8
      do 70 js = 1, jx-2
         i0 = i0 + iz
         i1 = i1 + iz
         i2 = i2 + iz
         i3 = i3 + iz
         Z(i0)   = dot6h(Z,AS(1,1,4),CS( 1,4),U(L),dot1,j0,j1)
         Z(i0+1) = dot6h(Z,AS(1,2,4),CS( 2,4),U(L),dot2,j0,j1)
         Z(i1)   = dot6h(Z,AS(1,3,4),CS( 4,4),U(L),dot3,j0,j1)
         Z(i1+1) = dot6h(Z,AS(1,4,4),CS( 7,4),U(L),dot4,j0,j1)
         Z(i2)   = dot6h(Z,AS(1,5,4),CS(11,4),U(L),dot5,j0,j1)
         Z(i2+1) = dot6h(Z,AS(1,6,4),CS(16,4),U(L),dot6,j0,j1)
         Z(i3)   = dot6h(Z,AS(1,7,4),CS(22,4),U(L),dot7,j0,j1)
         Z(i3+1) = eight*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                 - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
         j0 = j0 + 1
         j1 = j1 + 1
         L  = L  + 7
  70  continue
c								corner #4
      i0 = i0 + iz
      i1 = i1 + iz
      i2 = i2 + iz
      i3 = i3 + iz
      Z(i0)   = dot4c(Z,AC(1,1,4),CC( 1,4),U(L),dot1,j0,j1)
      Z(i0+1) = dot4c(Z,AC(1,2,4),CC( 2,4),U(L),dot2,j0,j1)
      Z(i1)   = dot4c(Z,AC(1,3,4),CC( 4,4),U(L),dot3,j0,j1)
      Z(i1+1) = dot4c(Z,AC(1,4,4),CC( 7,4),U(L),dot4,j0,j1)
      Z(i2)   = dot4c(Z,AC(1,5,4),CC(11,4),U(L),dot5,j0,j1)
      Z(i2+1) = dot4c(Z,AC(1,6,4),CC(16,4),U(L),dot6,j0,j1)
      Z(i3)   = dot4c(Z,AC(1,7,4),CC(22,4),U(L),dot7,j0,j1)
      Z(i3+1) = eight*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
     >              - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

      return
      end
