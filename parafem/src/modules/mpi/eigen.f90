MODULE eigen
USE precision ; USE global_variables; USE mp_interface; USE gather_scatter;
USE maths
CONTAINS
  subroutine lancz1(n,el,er,acc,leig,lx,lalfa,lp,iflag,u,v, &
                     eig,jeig,neig,x,del,nu,alfa,beta,v_store)            
 implicit none ! for use in parallelised programs - start addresses 1
!     taken from Parlett and Reid - see Harwell Library
!     use the Lanczos algorithm to find the spectrum of a symmetric matrix
!     or that part of the spectrum that lies in a given interval. this
!     version returns information for later calculation of eigenvectors
!     by lancz2.
      integer ::   n,leig,lx,lalfa,iflag,lp,neig
      integer ::   nu(lx), jeig(2,leig)
      real(iwp)::   el,er,acc,u(n),v(n), eig(leig),x(lx),del(lx), &
                        alfa(lalfa),beta(lalfa),v_store(n,lalfa)
!
! n must be set to the matrix order. it is not altered.
! el,er must be set to indicate the range within which the spectrum
!     is wanted. if el.ge.er then it is assumed that the whole spectrum
!     is wanted.they are not altered.
! acc must be set to the required precision, relative to the largest
!     eigenvalue of a. if it is very small or negative then as much
!     accuracy as the precision allows will be found. it is not altered.
! leig must be set to length of array eig. it must be as large as the
!     number of distinct eigenvalues in the range (el,er). it is not
!     altered.
! lx must be set to the length of arrays x,del,nu. a value three times
!     number of distinct eigenvalues in the range (el,er) usually
!     suffices. it is not altered.
! lalfa must be set to the length of arrays alfa and beta. it limits
!      the number of lanczos steps possible. it is not altered.
! lp must be set to the unit number for diagnostic messages. if lp.le.0
!     the messages are suppressed. it is not altered.
! iflag must be set prior to the first call for a particular matrix to
!       -1 if the user does not want to specify a start vector, or
!       -2 if the user has placed a required start vector in v.
!     it should not otherwise be changed.
!     on return it has the value
!       0 on successful completion
!       1 on a intermediate return
!       2,3,...,7 on error conditions
! u,v hold the working vectors of the lanczos recurrence. on a return
!     with iflag=1 the user must add to u the vector a*v without
!     altering v.
! eig need not be set by the user. on any return it contains neig
!     computed eigenvalues, stored in increasing order.
! jeig need not be set by the user. on return jeig(1,i) contains the
!     lanczos step at which eig(i) was accepted and jeig(2,i) contains
!     the matching point for the recurrences to get the eigenvector.
! neig need not be set by the user. on any return it contains the number
!     of eigenvalues in eig.
! x is used for a work array.  x(1).lt.x(2).lt.x(3).lt. ... are points
!     approximating eigenvalues of a.
! del is used as a workarray.  del(i),i=1,2,3,... contain the last pivot
!     in the lu factorization of t-lamda*i, where t is the lanczos
!     tridiagonal matrix, i is the identity matrix and lamda=x(i).
! nu is used as a workarray.  iabs(nu(i))-1,i=1,2,3,... contain the
!     number of negative pivots in the lu factorization of t-lamda*i,
!     lamda=x(i). fixed intervals x(i),x(i+1)  about converged
!     eigenvalues are flagged by negative values of nu(i).
! alfa and beta are used as workarrays holding the lanczos tridiagonal
!     matrix t. beta(1)=0.  and a typical row is
!         beta(i)     alfa(i)     beta(i+1)
      integer :: nfix = 0,nxtbnd,i,nlan,jlan1,jlan,mlan,l,kx,maxrz,         &
                 nuk,k=0,idummy,lfp,m=0,m1,j,ne=0,nxtr=0,match,nur,ncand,   &
                 nritz,maxrzo
      real(iwp) ::     drelpr,vv,uu,alf,bet,vi,tol,be2,anorm,oldel,         &
                       older,dr,xmid,xa,xc,fa,fc,tem,den,pi,sig,            &
                       tau,disc,xav,pole,ritz,xl,xr,diff,                   &
                       err,bnd,t,tolc,w,xli,xri,en,g,fa01as
      common/l_block/anorm,oldel,older,tolc,kx,maxrz,maxrzo,nlan,mlan,nxtbnd
! common block l_block contains the following quantities,
!     which are required to be preserved between entries:
!   anorm is an estimate of the norm of the matrix a, based
!      on a gregorshin bound applied to t.
!   oldel,older are the values of el,er on the previous entry.
!   tolc is the tolerance for agreement between sucessive ritz
!      values that decides whether to call the error bounding routine
!      ea15cd (but see also nxtbnd)
!   kx is the number of points currently stored in x, with associated
!      information in del and nu.
!   maxrz is the maximal number of ritz values in a gap between
!      fixed intervals.
!   maxrzo is the previous value of maxrz
!   nlan is the number of lanczos steps performed. alfa(i)<
!      beta(i),i=1,2,...,nlan have been set.
!   mlan is the order of tridiagonal matrix so far used for this
!      interval of the spectrum.
!   nxtbnd is used to delay error bounding when tolc is at roundoff
!      error level. ea15cd is not called until we are looking
!      at the tridiagonal matrix of order nxtbnd.
      real(iwp)  ::     zero,one,two,half,three
! drelpr is the relative precision
      data drelpr/2.2e-16/, zero/0.0e0/, one/1.0e0/, two/2.0e0/
!cray data drelpr/1.4e-14/, zero/0.0e0/, one/1.0e0/, two/2.0e0/
!ibm  data drelpr/2.2e-16/, zero/0.0e0/, one/1.0e0/, two/2.0e0/
      data half/0.5e0/, three/3.0e0/
      en=float(n)
      if(n.le.0 .or. lx.le.5 .or. lalfa.le.0)go to 950
      if(iflag.eq.1)go to 35
      vv=zero
      if(iflag.eq.-2)go to 5; if(iflag.eq.-1)go to 15; if(iflag.eq.0)go to 70
      go to 970
!
! find the norm of the user's start vector.
    5 vv = sum_p(v**2)
   10 continue
      if(vv.gt.zero)go to 25
!
! generate pseudo-random start vector
   15 g=1431655765.e0
      do 20 i=1,n
      g=mod(g*9228907.e0,4294967296.e0)
      if(i.ge.0)fa01as=g/4294967296.e0
      v(i)=fa01as
   20 continue
      vv = sum_p(v**2)
!
! normalize start vector and perform initializations
   25 vv=one/ sqrt(vv)
      v = v*vv  ;  u = zero
      nlan=0    ;  anorm=zero  ;  beta(1)=zero ;   neig=0
      go to 912
!
! perform a lanczos step
   35 nlan=nlan+1
      v_store(:,nlan) = v
      if(nlan.ge.lalfa)go to 940
      alf = DOT_PRODUCT_P(u,v) ;      alfa(nlan)=alf
      u = u - alf * v    ;  uu = norm_p(u)
      bet=beta(nlan)     ;  anorm=max(anorm,bet+ abs(alf)+uu)
      bet=max(uu,anorm*drelpr) ;      uu=one/bet  ;      beta(nlan+1)=bet
      do  i=1,n   ;   vi=u(i)*uu ; u(i)=-bet*v(i); v(i)=vi ; end do
!
! this is the beginning of the loop that analyses t
! normally we advance the analysis of t by one row, but the
!     user may have requested a restart
   70 jlan1=nlan
      if(oldel.eq.el .and. older.eq.er .and. iflag.eq.1)go to 71
      jlan1=1     ;      neig=0
   71 do 910 jlan=jlan1,nlan
      mlan=jlan
      if(jlan-2<0)go to 73 ; if(jlan-2==0)go to 76; if(jlan-2>0)go to 85
!
! jlan=1. check for the trivial case.
   73 eig(1)=alfa(1)
      jeig(1,1)=jlan      ;      jeig(2,1)=jlan
      if(beta(2).gt.drelpr*anorm)go to 910
      if(el.ge.er .or. (el.le.eig(1) .and. eig(1).le.er))neig=1
      go to 960
!
! jlan=2. set up four points that span and separate  the
!     three ritz values at levels 1 and 2
   76 w=half* abs(alfa(1)-alfa(2))
      t= sqrt(w*w+beta(2)**2)
      x(2)=half*(alfa(1)+alfa(2)-t-w)  ;      x(3)=x(2)+t+w
      x(1)=x(2)+w-1.1*t                ;      x(4)=x(3)-w+1.1*t
      do l=1,4;del(l)=alfa(2)-x(l)-beta(2)**2/(alfa(1)-x(l)); nu(l)=l-(l-1)/2
      end do
      kx=4  ;      maxrz=2    ;      maxrzo=1
      if(x(1).lt.x(2) .and. x(3).lt.x(4))go to 910
      eig(1)=(alfa(1)+alfa(2))*half
      jeig(1,1)=jlan    ;      jeig(2,1)=jlan   ;      neig=1
      go to 960
!
! jlan.gt.2
!
! add or remove points to the right, if appropriate
   85 nuk=jlan-1
      if(el.ge.er)go to 100
! find interval (x(k-1),x(k)) containing er
      do 92 i=2,kx
      k=2+kx-i            ;      if(x(k-1).le.er)go to 95
   92 continue
      k=1
   95 nuk=min0(iabs(nu(k)),jlan-1)
! look for point with nu value at least 2 greater than any in interval
!     containing er
      do 97 i=k,kx
      if(iabs(nu(i)).gt.nuk+1)go to 98
   97 continue
      go to 100
!  reduce kx to get rid of unnecessary points
   98 kx=i
      go to 110
! if necessary add extra points to right
  100 do 105 idummy=1,lx
      if(iabs(nu(kx-1)).gt.nuk)go to 110
      if(kx.ge.lx)go to 930
      kx=kx+1
      x(kx)=x(kx-1)*three-x(kx-2)*two
      call lsub_2(jlan-1,alfa,beta,x(kx),del(kx),nu(kx),dr,nur)
  105  continue
!
! copy information to ends of arrays x,del,nu. if any
!     gap between fixed intervals contains no ritz values,
!     remove all points in that gap
  110 lfp=lx
! lfp is the last left-end of a fixed interval
      m=lx
      do 140 idummy=1,lx
      if(nu(kx-1).ge.0)go to 130
      if(nu(kx).eq.iabs(nu(lfp)))m=lfp-1   ;      lfp=m-1
  130 x(m)=x(kx)
      del(m)=del(kx)   ;      nu(m)=nu(kx) ;      kx=kx-1  ;      m=m-1
      if(kx.eq.1)go to 150
  140 continue
  150 x(m)=x(1)
      del(m)=del(1)        ;      nu(m)=nu(1)
!
! add or remove points to the left, if appropriate
      nuk=2
      if(el.ge.er)go to 175
! find interval (x(k-1),x(k)) containing el
      m1=m+1
      do 155 k=m1,lx
      if(x(k).ge.el)go to 160
  155 continue
      k=lx+1
  160 nuk=max0(iabs(nu(k-1)),2)
! look for point with nu value at least 2 less than any in interval
!     containing el
      do 165 j=m1,k
      i=m+k-j
      if(iabs(nu(i)).lt.nuk-1)go to 170
  165 continue
      go to 175
! increase m to get rid of unnecessary points
  170 m=i
      go to 190
! if necessary, add extra points to left
  175 do 180 idummy=1,lx
      if(iabs(nu(m+1)).lt.nuk)go to 190
      if(m.eq.3)go to 930   ;      m=m-1
      x(m)=three*x(m+1)-two*x(m+2)
      call lsub_2(jlan-1,alfa,beta,x(m),del(m),nu(m),dr,nur)
  180 continue
  190 k=3
      if(m.le.5)go to 930
      xri=x(m)    ;  alf=alfa(jlan)  ;  bet=beta(jlan) ;   be2=bet**2
! tol holds the radius of fixed intervals set up around converged
!     eigenvalues.
      tol=max(acc,en*two*drelpr)*anorm
! tolc is the tolerance used for accepting eigenvalues.
      if(jlan.lt.10)tolc=max(tol**2*en/anorm,anorm*drelpr*5.)
! move leading three points to beginning of store and update
!     associated data.
      m=m-1
      do  210 i=1,3
        m=m+1   ;   x(i)=x(m)  ;   nu(i)=nu(m) ;  del(i)=alf-x(m)-be2/del(m)
        if(del(i))195,200,210
   195  nu(i)=nu(i)+isign(1,nu(i))
        go to 210
   200  del(i)=drelpr*anorm
   210 continue
!
! process the points from left to right. the current interval
!     is (x(k-1),x(k))=(x(m-1),x(m)), k.lt.m-2.
! nu(i),del(i),i=1,2,...,k hold new values of nu,delta.
! nu(i),del(i),i=m-2,m-1,...,lx hold old values.
      ne=1
! ne normally holds the number of empty intervals adjacent on the
!     left of the current interval.
      lfp=1
! lfp points to the last fixed point encountered.
!      do 530 until m.eq.lx
  220 if(m.ge.lx)go to 600
      if(nu(k-1).ge.0)go to 240
! we have a fixed interval
      if(iabs(nu(k-1)).eq.iabs(nu(k)))go to 230
! accept fixed interval
      lfp=k
      go to 500
! fixed interval no longer contains ritz value. free it.
  230 nu(k-1)=iabs(nu(k-1))
      do 235 i=1,neig
      if(eig(i).le.x(k))go to 235
      eig(i-1)=eig(i)  ;   jeig(1,i-1)=jeig(1,i) ;  jeig(2,i-1)=jeig(2,i)
  235 continue
      neig=neig-1
  240 if(nu(k-1).lt.iabs(nu(k)))go to 250
! current interval contains no ritz value.
      ne=ne+1
      if(ne.le.3)go to 502
! there are four adjacent empty intervals. remove middle point.
      x(k-2)=x(k-1) ;  nu(k-2)=nu(k-1); del(k-2)=del(k-1); x(k-1)=x(k)
      nu(k-1)=nu(k) ;  del(k-1)=del(k);  ne=3  ;  k=k-1
      go to 502
  250 if(iabs(nu(m-1)).eq.iabs(nu(m)))go to 500
!
! interval is 'interesting'. it contains at least one new and
!     at least one old ritz value.
! jump if interval is a gap between fixed intervals,
!     and contains just one ritz value, unless this is the
!     first time that all such gaps have less than two
!     ritz values.
      if(maxrz.le.1 .and. maxrzo.gt.1)go to 252
      if(nu(k-1)+nu(k).eq.-1)go to 500
  252 xmid=(x(m-1)+x(m))*half
! bisect if interval contains more than one ritz value.
      ritz=xmid
      if(iabs(nu(k)).gt.nu(k-1)+1)go to 275
! check whether progress is so slow that bisection is needed
!     the criterion is that three steps have been taken without
!     reducing the interval length by at least the factor 0.6.
!     (xli,xri) is the original interval of interest or the
!     interval of interest when it was last reduced in length
!      by at least the factor 0.6
      if(x(k-1).ge.xri)go to 255
      if(x(k)-x(k-1).ge.0.6*(xri-xli))go to 257
  255 nxtr=0
      xri=x(k-1)     ;      xli=x(k)
! bisect if progress on this interval is slow
  257 if(nxtr.ge.2)go to 275
      nxtr=nxtr+1
! estimate position (pole) of old ritz value by 2-1 rational
!     interpolation at x(l-1), x(l), x(l+1)
! choose l so that extra point is near
      l=m
      if(xmid-x(m-2).lt.x(m+1)-xmid)l=m-1
      xa=x(l-1)-x(l)  ;   xc=x(l+1)-x(l) ;  fa=xa*(xa+del(l-1))
      fc=xc*(xc+del(l+1)); tem=xa/(xa-xc)
      den=del(l)-del(l-1)-tem*(del(l+1)-del(l-1))
      if(den.eq.zero)go to 275
      pi=(tem*(fa-fc)-fa)/den
      sig=half*((fc-fa)-pi*(del(l+1)-del(l-1)))/(xc-xa)
      tau=fa-two*xa*sig-del(l-1)*pi   ;      disc=sig**2+tau
      if(disc.lt.zero)go to 275
      xav=sig+ sign( sqrt(disc),sig)
      pole=-tau/xav
! jump if pole is in required interval
      if( abs(pole+x(l)-xmid).le.x(m)-xmid)go to 260
! try other root.
      pole=xav
      if( abs(pole+x(l)-xmid).gt.x(m)-xmid)go to 275
! estimate position (ritz) of ritz value by 2-1 rational
!     interpolation at x(k-1),x(k) using (x-pole) for
!     denominator.
  260 if(l.ne.m)go to 265
      tau=-pole*del(k)
      sig=xa+((xa-pole)*del(k-1)-tau)/xa
      go to 270
  265 tau=-pole*del(k-1)  ;  sig=xc+((xc-pole)*del(k)-tau)/xc
  270 sig=sig*half        ;      disc=sig**2+tau
      if(disc.lt.zero)go to 275
      xav=sig+ sign( sqrt(disc),sig)
      tau=-tau/xav
      if( abs(tau+x(l)-xmid).le.x(m)-xmid)go to 280
      tau=xav
      if( abs(tau+x(l)-xmid).le.x(m)-xmid)go to 280
! calculation has failed
! if third point is just outside current interval then take
!     point twice as far inside interval. otherwise bisect.
  275 ritz=xmid
      pole=xmid
      if(x(m)-x(m-1).le.tolc)go to 300
      go to 490
  280 ritz=x(l)+tau   ;   diff= abs(tau-pole) ;      pole=x(l)+pole
      if( abs(tau).ge.50.*max(diff,tolc))go to 480
      if(diff.gt.tolc)go to 500
! we may have a converged eigenvalue
  300 if(jlan.lt.nxtbnd)go to 500
      call lsub_1(jlan,lalfa,alfa,beta,ritz,anorm*drelpr*en,err,bnd,match)
      if(bnd.lt.err*0.1)tolc=diff*(tol/err)**2
      if(tolc.gt.anorm*drelpr*5.0)go to 305
      tolc=anorm*drelpr*5.0
! tolc has become too small to give good criterion for calling ea15cd.
!    do not call lsub_1 again until 1% more steps performed.
      nxtbnd=jlan+jlan/100
  305 if(err.le.tol)go to 310
      if(bnd.gt.err*0.1 .and. x(m-1).lt.ritz .and. ritz.lt.x(m)) go to 480
      go to 500
!
! we have an accepted point
! test whether (ritz-err,ritz+err) overlaps a fixed interval
  310 if(ritz-err.lt.x(lfp) .and. lfp.gt.1)go to 500
      do  i=m,lx
       if(x(i).gt.ritz+err)go to 420   ;      if(nu(i).lt.0)go to 500
      end do
! set up new fixed interval
  420 xl=ritz-tol
      if(xl.le.x(lfp))go to 450
! remove points in interval (xl,ritz)
      do 430 idummy=1,lx
      if(x(k-1).lt.xl)go to 440
      k=k-1
  430 continue
  440 x(k)=xl
      call lsub_2(jlan,alfa,beta,xl,del(k),nu(k),dr,nur)
      go to 455
  450 k=lfp
  455 nu(k)=-nu(k)
      xr=ritz+tol
! remove points in interval (ritz,xr)
      i=m
      do  m=i,lx
       if(x(m).gt.xr)go to 470   ;      if(nu(m).lt.0)go to 473
      end do
      m=lx+1
  470 m=m-1
      x(m)=xr
      call lsub_2(jlan-1,alfa,beta,x(m),del(m),nu(m),dr,nur)
  473 nfix=nfix+1
      if(neig.ge.leig)go to 920    ;      if(el.ge.er)go to 474
      if(ritz.lt.el .or. ritz.gt.er)go to 478
  474 i=neig
      if(neig.le.0)go to 477
      do 475 j=1,neig
      if(eig(i).lt.ritz)go to 477
      eig(i+1)=eig(i); jeig(1,i+1)=jeig(1,i) ; jeig(2,i+1)=jeig(2,i); i=i-1
  475 continue
  477 eig(i+1)=ritz
      jeig(1,i+1)=jlan  ;  jeig(2,i+1)=match;  neig=neig+1
  478 if(m.gt.k+2)go to 505
      go to 930
! extrapolate to estimate position of eigenvalue
  480 ritz=ritz+ritz-pole
      ritz=max(ritz,x(m-1)+(x(m)-x(m-1))*0.01)
      ritz=min(ritz,x(m)-(x(m)-x(m-1))*0.01)
! insert new point
  490 if(k.ge.m-3)go to 930
      m=m-1 ; x(m-2)=x(m-1); del(m-2)=del(m-1);  nu(m-2)=nu(m-1)
      x(m-1)=x(m) ; del(m-1)=del(m); nu(m-1)=nu(m); x(m)=ritz; x(k)=ritz
      call lsub_2(jlan,alfa,beta,ritz,del(k),nu(k),dr,nur)
      del(m)=dr      ;      nu(m)=nur
      go to 530
! interval is acceptable. advance by one point
  500 ne=0
  502 m=m+1
  505 k=k+1
      x(k)=x(m)       ;   del(k)=alf-x(k)-be2/del(m) ;    nu(k)=nu(m)
      if(del(k))510,520,530
  510 nu(k)=nu(k)+isign(1,nu(k))
      go to 530
  520 del(k)=drelpr*anorm
  530 go to 220
!
! scan to find the maximum number of ritz values (maxrz) in a gap
!     between fixed intervals, the number of gaps containing candidates
!     that are being watched (ncand) and the number of fixed intervals
!     (nfix).
  600 kx=k
      ncand=0    ;   maxrzo=maxrz ;   maxrz=0  ;  nfix=0 ; lfp=2 ;  k=k-1
      do 710 i=1,k
      if(nu(i).gt.0)go to 710
      nfix=nfix+1  ; nritz=iabs(nu(i))-iabs(nu(lfp)); maxrz=max0(maxrz,nritz)
      if(nritz.gt.0 .and. i.gt.lfp+1)ncand=ncand+1
      do 680 lfp=i,kx
      if(nu(lfp).gt.0)go to 710
  680 continue
      lfp=kx
  710 continue
      nritz=iabs(nu(k))-iabs(nu(lfp)) ;   maxrz=max0(maxrz,nritz)
      if(nritz.gt.0 .and. k.gt.lfp+1)ncand=ncand+1
      if(bet.le.en*anorm*drelpr)go to 960
      if(maxrz.gt.1 .or. ncand.gt.0)go to 910
      if(nfix.gt.0)go to 915
  910 continue
!
! normal returns
  912 iflag=1
      go to 1000
  915 iflag=0
      go to 1000
!
! error returns
  920 iflag=2
      if(numpe==npes.and.lp.gt.0)write(lp,925)'error return 2 from lancz1.  &
                                      &         leig=',leig
  925 format( a , i6)
      go to 1000
  930 iflag=3
      if(numpe==npes.and.lp.gt.0)write(lp,935)'error return 3 from lancz1.  &
                                      &         lx=',lx
  935 format( a , i6)
      go to 1000
  940 iflag=4
      if(numpe==npes.and.lp.gt.0)write(lp,945)'error return 4 from lancz1.  &
                                      &         lalfa=',lalfa
  945 format( a , i6)
      go to 1000
  950 iflag=5
      if(numpe==npes.and.lp.gt.0)write(lp,955)'error return 5 from lancz1.  &
                                       &        n,lx,lalfa=', n,lx,lalfa
  955 format( a , 3i6)
      go to 1000
  960 iflag=6
      if(numpe==npes.and.lp.gt.0)write(lp,965)'error return 6 from lancz1.  &
                         &   premature termination.',neig,'eigenvalues found'
  965 format( a , i6 , a)
      go to 1000
  970 if(numpe==npes.and.lp.gt.0)write(lp,975)'error return 7 from lancz1.  &
                             & on entry iflag is ',iflag
  975 format( a , i4)
      iflag=7
 1000 oldel=el
      older=er
      return
      end subroutine lancz1
      subroutine lsub_1(j,lalfa,alfa,beta,e,enorm,err,bnd,match)
! given a tridiagonal matrix t generated by the lanczos process applied
!     to a given symmetric matrix a and an approximate eigenvalue of t,
!     find a bound for its error as an eigenvalue of a.
! j is the index of the lanczos step. it is not altered.
! lalfa must be set to the lengths of arrays alfa and beta. it is not
!     altered.
! alfa(i),i=1,j must be set to the diagonal elements of the lanczos
!     tridiagonal matrix t. they are not altered.
! beta(i),i=2,j+1 must be set to the off-diagonal elements of the
!     lanczos tridiagonal matrix t. they are not altered.  note that a
!     typical row of t is
!          beta(i)  alfa(i)  beta(i+1)
! e must be set to the approximate eigenvalue of t. it is not altered.
! enorm must be set to the relative precision times the norm of a times
!     the order of a. it is not altered.
! err is set to a bound for the error of e.
! bnd is set to an estimate of the amount that err might be reduced by
!     using a more accurate e and/or a more accurate eigenvector
!     calculation.
! match is set to the matching point between use of forward and backward
!     recurrences.
   implicit none
      integer      ::       lalfa,j,match
      real(iwp)    ::       alfa(lalfa),beta(lalfa),e,s,psi,bnd,err,enorm
      real(iwp)    ::       w,w0,w1,w2,sw,z1,z2,bndk,errk,wmax10  ,zero,one
      integer      ::       j1,k,ktry,k1,kk,i
      data zero/0.0e0/, one/1.0e0/
      j1=j-1
!
! start the forward recurrence for solving (t-e)*z=psi*ek
      w1=zero  ; w2=one ; sw=zero ;    k=0 ;      w=alfa(1)-e ;  wmax10=zero
!
      do 100 ktry=1,10
      if(k-j1<0)goto 5; if(k-j1==0)goto 15; if(k-j1>0)goto 110
    5 k1=k+1
!     continue the forward recurrence until the first maximum
!     or a maximum at least ten times bigger than the previous
!     maximum is found.
      do 10 k=k1,j1
      w0=w1 ;   w1=w2;   w2=-w/beta(k+1) ;   sw=sw+w1**2
      w=(alfa(k+1)-e)*w2+beta(k+1)*w1
      if( abs(w1).lt.wmax10)go to 10
      if( abs(w1).ge.max( abs(w0), abs(w2)))go to 20
   10 continue
   15 k=j
      w0=w1    ;      w1=w2  ;      sw=sw+w1**2
      if( abs(w1).lt.wmax10)go to 100
      if( abs(w1).lt. abs(w0))go to 100
!
!     forward recurrence shows a maximum at component k. continue
!     solving (t-e)*z=psi*ek by backward recurrence
   20 s=zero
      wmax10= abs(w1)*10.     ;      z1=one ;      psi=alfa(j)-e
      if(k.gt.j1)go to 40
      do  kk=k,j1
      i=k+j1-kk   ;      z2=z1  ;      z1=-psi/beta(i+1)  ;      s=s+z2**2
      psi=(alfa(i)-e)*z1+beta(i+1)*z2
      end do
   40 if(w1.eq.zero)go to 100
!
! form present error bound and store the best so far found
      z1=z1/w1  ;      s=s+sw*z1**2   ;      s= sqrt(s)
      if(k.gt.1)psi=psi+beta(k)*z1*w0
      bndk=1.25* abs(psi)/s  ;      errk=1.25*(enorm+beta(j+1)/s)+bndk
      if(ktry.eq.1)go to 50
      if(err.lt.errk)go to 90
   50 err=errk
      bnd=bndk     ;      match=k
!
! test result for acceptability
   90 if(bndk.lt.errk/5.)go to 110
  100 continue
  110 return
      end subroutine lsub_1
      subroutine lsub_2(n,alfa,beta,x,del,nu,dr,nur)
! computes the triangular factorization of t-x*i, where t is a
!     tridiagonal matrix, x is a scalar and i is the identity matrix. it
!     records the number of negative pivots and the last two pivots.
  implicit none
      integer      :: n,nu,nur
      real(iwp)    :: alfa(n),beta(n),x,del,dr
! n must be set to the order of the matrix. it is not altered.
! alfa must be set ot contain the diagonal elements of the tridiagonal
!     matrix. it is not altered.
! beta(1) must have the value zero and beta(2),...,beta(nlan) must be
!     set to contain the off-diagonal elements. beta is not altered.
! x   must be set to the scalar. it is not altered.
! del is set to the last pivot.
! nu is set to 1+(the number of negative pivots).
! dr is set to the last-but-one pivot.
! nur is set to 1+(the number of negative pivots when the last is
!     excluded).
      real(iwp) ::      zero,one        ; integer :: k
! drelpr is the relative precision.
      real(iwp) ::      drelpr
!ibm  data drelpr/2.2e-16/
!ray  data drelpr/1.4e-14/
      data drelpr/2.2e-16/
      data zero/0.0e0/, one/1.0e0/
      nu=1          ;      del=one
      do 10 k=1,n
      dr=del         ;      del=alfa(k)-x-beta(k)**2/del
      if(del)6,3,10
    3 del=drelpr*( abs(alfa(k))+ abs(x))
      go to 10
    6 nu=nu+1
   10 continue
   20 nur=nu
      if(del.lt.zero)nur=nu-1
      return
      end subroutine lsub_2
      subroutine lancz2(n,lalfa,lp,eig,jeig,neig,  &
                       alfa,beta,lz,jflag,y,w,z,v_store)
! find approximate eigenvectors corresponding to eigenvalues found
!     by lancz1.
! n must be set to the matrix order. it is not altered.
! lalfa must be set to the length of arrays alfa and beta. it
!     is not altered.
! lp must be set to the unit number for diagnostic messages. if lp.le.0
!     the messages are suppressed. it is not altered.
! eig must contain eigenvalues as found by lancz1. it is not altered.
! jeig must hold the information passed down by lancz1, namely jeig(1,i)
!     must contain  the lanczos step at which eig(i) was accepted and
!     jeig(2,i) must contain the matching point for the eigenvector
!     corresponding to eig(i). jeig is not altered.
! neig must hold the number of eigenvectors wanted. it is not altered.
! alfa,beta must be as left by lancz1. they are not altered.
! lz must hold the first dimension of z. it must be at least
!     max(jeig(1,i),i=1,neig). it is not altered.
! iflag need not be set. on return it has one of the values
!      0    successful completion
!      1    n.le.0 .or. lalfa.le.0 .or. neig.le.0
!      2    ly too small
!      3    lz too small
! y is set to the required eigenvectors.
! w is used as a workvector.
! z is used as a workvector.
! lz must be set to the first dimension of z. it must be at least
!     max(jeig(1,i),i=1,2,...,neig).
   implicit none 
      integer   ::       neig,lalfa,n,lz,jflag,lp  
      real(iwp) ::       eig(neig),alfa(lalfa),beta(lalfa),   &
                         y(n,neig),w(n),z(lz,neig),v_store(n,lalfa)
      integer   ::       jeig(2,neig)
      real(iwp) ::       s,zero,one  ; integer :: i,j,l,m
      data zero/0.0e0/, one/1.0e0/
      jflag=0
      if(n.le.0 .or. lalfa.le.0 .or. neig.le.0)go to 100
      if(n.lt.n)go to 120
      m=0
      do 20 l=1,neig
      m=max0(m,jeig(1,l))  ;      if(lz.lt.m)go to 140  ;      y(:,l) = zero
      call lsub_3(jeig(1,l),lalfa,alfa,beta,eig(l),jeig(2,l),z(1,l))
   20 continue
      do 50 j=1,m
      w = v_store(:,j)
      do 40 l=1,neig
      if(j.gt.jeig(1,l))go to 40
      do 30 i=1,n
      y(i,l)=y(i,l)+z(j,l)*w(i)
   30 continue
   40 continue
   50 continue
!
! normalize the vectors
      do  l=1,neig
        s = one/norm_p(y(:,l))   ;     y(:,l) = y(:,l)*s
      end do  
      go to 150
  100 jflag=1
      if(numpe==npes.and.lp.gt.0)write(lp,110)'error return 1 from lancz2. &
          & n=',n,'lalfa=',lalfa,'neig=',neig
  110 format( 3(a,i6))
      go to 150
  120 jflag=2
      if(numpe==npes.and.lp.gt.0)write(lp,130)'error return 2 from lancz2. &
          &   ly =',n, 'and should be at least',n
  130 format( 2(a , i6))
      go to 150
  140 jflag=3
      if(numpe==npes.and.lp.gt.0)write(lp,145)'error return 3 from lancz2. &
         &    lz =',lz, 'and should be at least',m
  145 format( 2(a , i6 ))
  150 return
      end  subroutine lancz2
      subroutine lsub_3(j,lalfa,alfa,beta,e,match,z)
! find an approximate eigenvector of a lanczos tridiagonal matrix by
!     forward and backward recurrence.
! j is the order of the matrix and is not altered.
! alfa,beta must be exactly as left by lancz1. they are not altered.
! e is the eigenvalue corresponding to which an eigenvector is wanted.
!     it is not altered.
! match is the matching point between the recurrences. it is not
!     altered.
! z is set to the required eigenvector.
  implicit none 
      integer      ::   j,lalfa,match
      real(iwp)    ::   alfa(lalfa),beta(lalfa),e,z(j)
      real(iwp)    ::   psi,w,w0,w1,w2,z1,z2  ,zero,one
      integer      ::   j1,k,kk,i 
      data zero/0.0e0/, one/1.0e0/
      j1=j-1
!
! forward recurrence.
      w1=zero  ;      w2=one  ;      w=alfa(1)-e
      if(j1.le.0)go to 15
      do 10 k=1,j1
      w0=w1    ;      w1=w2  ;      w2=-w/beta(k+1)
      w=(alfa(k+1)-e)*w2+beta(k+1)*w1   ;      z(k)=w1
      if(k.eq.match)go to 20
   10 continue
   15 k=j
      w0=w1      ;      w1=w2    ;      z(k)=w1
! backward recurrence
   20 z1=one
      psi=alfa(j)-e
      if(k.gt.j1)go to 40
      do  kk=k,j1
      i=k+j1-kk  ;      z2=z1 ;      z1=-psi/beta(i+1)
      z(i+1)=z2          ;      psi=(alfa(i)-e)*z1+beta(i+1)*z2
      end do
! rescale the last set of components
   40 w1=w1/z1
      if(k.gt.j1)go to 60
      do  i=k,j1
       z(i+1)=z(i+1)*w1
      end do
   60 return
      end subroutine lsub_3 
END MODULE eigen
