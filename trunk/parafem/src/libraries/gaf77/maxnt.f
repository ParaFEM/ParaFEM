c  ********************************************************************
c  *                                                                  *
c  *                           Subroutine maxnt                       *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, Princeton, May 4, 1989.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  compute the maximum entropy power spectral estimate
c
c  This routine calculates the power spectral distribution of
c  an input vector X using the maximum entropy method developed
c  by Burg (1970), Akaike (1970), and Kanasewich (1975).
c  Arguments are as follows;
c
c    X   input real vector of length NSTEP which contains the
c        observed series. On output, X will contain the one-sided power
c        spectral estimates in the first NW+1 locations (NW < NSTEP).
c
c  B1,B2 two temporary real vectors of length at least NSTEP
c
c  A,AA  two temporary real vectors of length at least NSTEP/2
c
c  DT    real `time' increment between successive observations
c
c  NSTEP number of time steps in the observed series
c
c  NW    number of frequencies at which to calculate the power spectral
c        estimates. These frequencies lie between 0 and pi/DT inclusive.
c
c  ierr  integer error flag which is set to one of the following values;
c        =  0  if all goes well
c        = -1  if the number of frequencies selected > (nstep - 1)
c              (solution continues with nw = nstep - 1)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c--------------------------------------------------------------------------
      subroutine maxnt( x, b1, b2, a, aa, dt, nstep, nw, ierr )
      dimension x(*), b1(*), b2(*), a(*), aa(*)
      data pi/3.1415926535897932384/
      data zero/0.0/, one/1.0/
c					check data
      nww = nw
      if( nw .ge. nstep ) then
           ierr = -1
           nww  = nstep - 1
      endif
c					initialize variables
      mmax = nstep/2
      pm   = zero
      do 10 i = 1, nstep
         pm = pm + x(i)*x(i)
  10  continue
      pm = pm/float(nstep)
      b1(1) = x(1)
      b2(nstep-1) = x(nstep)
      do 20 i = 2,nstep-1
         b1(i) = x(i)
         b2(i-1) = x(i)
  20  continue
      r = zero
      q = zero
      do 30 i = 1,nstep-1
         r = r + b1(i)*b2(i)
         q = q + b1(i)*b1(i) + b2(i)*b2(i)
  30  continue
      a(1) = 2.*r/q
      pm = pm*(one - a(1)*a(1))
      fpel = (pm*pm)*float(nstep + 2)/float(nstep - 2)

c					calculate digital filter coefficients
      do 70 m = 2,mmax
c
         do 40 i = 1,m-1
            aa(i) = a(i)
  40     continue
         r = zero
         q = zero
         do 50 i = 1,nstep-m
            b1(i) = b1(i) - aa(m-1)*b2(i)
            b2(i) = b2(i+1) - aa(m-1)*b1(i+1)
            r = r + b1(i)*b2(i)
            q = q + b1(i)*b1(i) + b2(i)*b2(i)
  50     continue
         a(m) = 2.*r/q
         pmn = pm*(one - a(m)*a(m))
c
         do 60 i = 1,m-1
            a(i) = aa(i) - a(m)*aa(m-i)
  60     continue
c					minimize final prediction error

         fpe = (pmn*pmn)*float(nstep+m+1)/float(nstep-m-1)
         if(fpe.gt.fpel)go to 80
         fpel = fpe
         pm = pmn
  70  continue

      m = mmax + 1
  80  m = m-1
      dw = pi/(dt*float(nww))
c					calculate power spectral values
      r = zero
      do 90 i = 1,m
         r = r + aa(i)
  90  continue
      q = (one - r)*(one - r)
      x(1) = zero
      if(q.ne.zero) x(1) = pm*dt/(q*pi)
      pm = pm*dt/pi
      w = zero

      do 110 j = 2,nww+1
         w = w + dw
         s = w*dt
         r = zero
         q = zero
         do 100 i = 1,m
            t = s*float(i)
            r = r + aa(i)*cos(t)
            q = q + aa(i)*sin(t)
  100    continue
         qw = (one - r)*(one - r) + q*q
         x(j) = pm/qw
  110 continue

      return
      end
