C
C    Implicit (semi-implicit with nonlinear hardening) integration for 
C    isotropic and kinematic hardening elasto-viscoplasticity
C      
C       
      IMPLICIT REAL*8(A-H,O-Z)
C
      OPEN (UNIT=11, FILE='stress.dat', STATUS='old')
C
      alpha = 0.2418e-6
      beta = 0.1135
      E = 1000.
      xnu = 0.3
      G = E/(2.*(1+xnu))
      yield = 10.
      h = 100.
C      h = 0.
      b = 0.1
      q = 100.
C      c = 100.
      c = 0.
C
      r = 0.
C
      dt = 1.E-1
      stranrat = 1.e-3
      stran = 0.
C
      do i = 1,10000
      dstran = stranrat*dt
      dpstress = E*dstran
      pstress = stress + E*dstran
C
C   yield function
C
      f = pstress - x - r - yield
      if(f.ge.0.) then
      r0 = r
      x0 = x
      dp = 0.
C
        do j = 1,10
C
C    linear hardening
C
        gamma = h
C
C    nonlinear hardening (semi-implicit integration)
C
C        gamma = b*(q-r)
C
        xphidp = -3.*G*alpha*beta*cosh(beta*
     +  (pstress-x-(3.*G*dp)-r-yield))
        xphi =  alpha*sinh(beta*(pstress-x-(3.*G*dp)-r-yield))
        xphir = -alpha*beta*cosh(beta*(pstress-x-(3.*G*dp)-r-yield))
        xphis = -alpha*beta*cosh(beta*(pstress-x-(3.*G*dp)-r-yield))
        deqpl = (xphi - (dp/dt))/((1/dt) - xphidp-gamma*xphir+c*xphis)
        dp = dp + deqpl
        r = r0 + gamma*dp
        x = x0 + c*dp
        if (deqpl.lt.1.e-12) goto 10
        end do
   10 continue
      iter = j
      end if
C
      stran = stran + dstran
      destran = dstran - dp 
      dstress = E*destran
      stress = stress + dstress
C
      ICOUNT = ICOUNT + 1
      IF(ICOUNT.EQ.100) THEN
      ICOUNT = 0
      WRITE(11,*) STRAN, STRESS, iter
      END IF
C
       end do
C
      CLOSE(11)
      end  