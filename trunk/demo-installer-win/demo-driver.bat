@ECHO OFF

SET PARAFEM_HOME=@INSTDIR@
cd "%PARAFEM_HOME%\bin"

:start
cls
ECHO -------------------------------------------------------------------------------
ECHO Welcome to the ParaFEM Demonstrator
ECHO -------------------------------------------------------------------------------
ECHO This package provides some sample implementations of software utilising ParaFEM
ECHO NB. To build your own software ParaFEM must be installed separately
ECHO -------------------------------------------------------------------------------

:demos
ECHO.
ECHO The following demos are available:
ECHO P121 - Three-dimensional analysis of an elastic solid
ECHO P122 - Three-dimensional analysis of an elasoplastic (Mohr-Colomb) solid
ECHO P126 - Three-dimensional steady-state Navier-Stokes analysis
ECHO.
SET /P DEMO="Which demo would you like to run?: (Enter name, e.g. P121): "
IF /I '%DEMO%'=='p121' goto :p121
IF /I '%DEMO%'=='p122' goto :p122
IF /I '%DEMO%'=='p126' goto :p126
goto demos

:p121
  cls
  call p121.bat
  pause
  goto start

:p122
  cls
  call p122.bat
  pause
  goto start

:p126
  cls
  call p126.bat
  pause
  goto start

