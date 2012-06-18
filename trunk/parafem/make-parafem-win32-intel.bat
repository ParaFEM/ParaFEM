echo
rem  *** EDIT THE NEXT THREE LINES IF NECESSARY ***
rem
rem
set PARAFEM=c:\parafem\parafem
rem !set G95=c:\MinGW\g95\bin\g95
set G95=ifort /integer-size:32 /real-size:64 /Od
set FFLAGS=-c 
set AR=lib /out:
set VER=1
rem
rem 
rem *** CLEAN ***
rem
rem
cd /D %PARAFEM%\bin
del *.a *.obj *.mod *.exe
cd /D %PARAFEM%\lib
del *.a
rem
rem
rem *** BUILD DUMMY_MPI LIBRARY ***
rem
rem
cd /D %PARAFEM%\src\libraries\dummy_mpi
del *.lib *.obj *.mod
%G95% %FFLAGS% mpi_stubs.f90
%AR%libmpi_stubs.lib mpi_stubs.obj 
move libmpi_stubs.lib %PARAFEM%\lib
del *.obj
rem
rem
rem *** BUILD GAF77 LIBRARY ***
rem
rem
cd /D %PARAFEM%\src\libraries\gaf77
del *.lib *.obj *.mod
%G95% %FFLAGS% *.f
%AR%libgaf77.lib *.obj
move libgaf77.lib %PARAFEM%\lib
del *.obj
rem
rem *** BUILD PARAFEM LIBRARY ***
rem
cd /D %PARAFEM%\src\modules\shared
del *.lib *.obj *.mod
cd /D %PARAFEM%\src\modules\mpi
del *.lib *.obj *.mod
cd /D %PARAFEM%\src\modules\shared
%G95% %FFLAGS%  precision.f90
%G95% %FFLAGS%  global_variables.f90
%G95% %FFLAGS%  geometry.f90
%G95% %FFLAGS%  elements.f90
%G95% %FFLAGS%  steering.f90
%G95% %FFLAGS%  timing.f90
%G95% %FFLAGS%  partition.f90
%G95% %FFLAGS%  plasticity.f90
%G95% %FFLAGS%  fluid.f90
move *.obj %PARAFEM%\src\modules\mpi
move *.mod %PARAFEM%\src\modules\mpi
cd /D %PARAFEM%\src\modules\mpi
%G95% %FFLAGS%  mp_interface.f90
%G95% %FFLAGS%  maths.f90
%G95% %FFLAGS%  input.f90
%G95% %FFLAGS%  output.f90
%G95% %FFLAGS%  gather_scatter.f90
%G95% %FFLAGS%  loading.f90
%G95% %FFLAGS%  large_strain.f90
%G95% %FFLAGS%  pcg.f90
%G95% %FFLAGS%  bicg.f90
%AR%libparafem.lib *.obj
move libparafem.lib %PARAFEM%\lib
move *.mod %PARAFEM%\include
del *.obj 
rem
rem
rem *** BUILD RFEMBC ***
rem
rem
cd /D %PARAFEM%\src\tools\preprocessing\rfembc
del rfembc.exe rfembc.o
%G95% rfembc.f90 %PARAFEM%\lib\libparafem.lib -I%PARAFEM%\include -o rfembc.exe
del *.obj *.lib
move rfembc.exe %PARAFEM%\bin
rem
rem 
rem *** BUILD RFEMFIELD ***
rem
rem
cd /D %PARAFEM%\src\tools\preprocessing\rfemfield
del rfemfield.exe rfemfield.obj
%G95% rfemfield.f90 %PARAFEM%\lib\libgaf77.lib %PARAFEM%\lib\libparafem.lib -I%PARAFEM%\include -o rfemfield.exe
del *.obj *.lib
move rfemfield.exe %PARAFEM%\bin
rem
rem *** BUILD RFEMSOLVE ***
rem
cd /D %PARAFEM%\src\programs\rfem
%G95% rfemsolve.f90 %PARAFEM%\lib\libmpi_stubs.lib %PARAFEM%\lib\libparafem.lib -I%PARAFEM%\include -o rfemsolve.exe
del *.obj *.lib
move rfemsolve.exe %PARAFEM%\bin
cd /D %PARAFEM%

