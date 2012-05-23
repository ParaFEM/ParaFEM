c  *********************************************************************
c  *                                                                   *
c  *                          Function second                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 3.0
c  Written by Gordon A. Fenton, TUNS, 1990
c
c  PURPOSE returns elapsed user execution time in seconds.
c
c  Returns elapsed user execution time of the calling process in seconds.
c  This routine tends to be sysem specific and you may need to customize it
c  for your environment. This routine works for Sun's and VaX's running
c  Ultrix f77.
c---------------------------------------------------------------------------
      real function second()
      real tyme(2)

      second = etime(tyme)
C                               and for the HP-9000
c     second = 1.e-6*float(clock())

      return
      end
