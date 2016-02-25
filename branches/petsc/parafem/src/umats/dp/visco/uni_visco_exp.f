C
C    Explicit integration for isotropic and kinematic hardening
C    elasto-viscoplasticity
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
C      h = 100.
      h = 0.
      b = 0.1
      q = 100.
      c = 100.
C
      r = 0.
      x = 0.
C
      dt = 1.E-1
      stranrat = 1.e-3
      stran = 0.
C
      do i = 1,10000
      dstran = stranrat*dt
      dpstress = E*dstran
      pstress = E*stran
C
C   yield function
C
      f = stress - x - r - yield
      
      if(f.ge.0.) then
      xphi = alpha*sinh(beta*(stress-x-r-yield))
        dp = xphi*dt
      end if
C
      stran = stran + dstran
      destran = dstran - dp 
      dstress = E*destran
      stress = stress + dstress
C
C    linear isotropic hardening
C
      gamma = h
C
C    nonlinear isotropic hardening hardening 
C
C      gamma = b*(q-r)
C
      dr = gamma*dp
      r = r + dr
C
C    linear kinematic hardening
C
      dx = c*dp
      x = x + dx
C
      ICOUNT = ICOUNT + 1
      IF(ICOUNT.EQ.100) THEN
      ICOUNT = 0
      WRITE(11,*) STRAN, STRESS
      END IF
C
       end do
C
      CLOSE(11)
      end  