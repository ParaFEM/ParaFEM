rem --------------------------------------------------------------------------------
rem Build ParaFEM and the 5th Edition example driver programs for Microsoft Windows
rem J.Robinson@software.ac.uk
rem Requires:
rem		MS Windows XP or newer
rem 	Visual Fortran Composer XE 2013
rem		MS Visual Studio 2010 Express
rem		Microsoft Visual C++ 2010 SP1 Redistributable Package
rem		OpenMPI 1.6 (Native Windows version NOT Cygwin)
rem --------------------------------------------------------------------------------

rem  *** EDIT THE FOLLOWING LINES IF NECESSARY ***
rem
rem
set PARAFEM=%CD%
set G95=mpif77 /integer-size:32 /real-size:64 /Od
set FFLAGS=-c -I%PARAFEM%\include
set AR=lib /out:
set BUILD_GROUP_5ED=p121 p122 p123 p124 p125 p126 p128 p129 p1210

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
%G95% %FFLAGS%  new_library.f90

move *.obj %PARAFEM%\src\modules\mpi
move *.mod %PARAFEM%\src\modules\mpi
cd /D %PARAFEM%\src\modules\mpi
%G95% %FFLAGS%  mpi_wrapper.f90
rem %G95% %FFLAGS%  mp_interface.f90 -I../../libraries/dummy_mpi/include
%G95% %FFLAGS%  mp_interface.f90
%G95% %FFLAGS%  maths.f90
%G95% %FFLAGS%  input.f90
%G95% %FFLAGS%  output.f90
%G95% %FFLAGS%  gather_scatter.f90
%G95% %FFLAGS%  loading.f90
rem %G95% %FFLAGS%  large_strain.f90
%G95% %FFLAGS%  pcg.f90
%G95% %FFLAGS%  bicg.f90
%G95% %FFLAGS%  eigen.f90

%AR%libparafem.lib *.obj
move libparafem.lib %PARAFEM%\lib
move *.mod %PARAFEM%\include
del *.obj

rem
rem
rem *** Build drivers ***
rem
rem
for %%f in (%BUILD_GROUP_5ED%) do ( cd /D src\programs\5th_ed\%%f
	call make clean
	call make
	cd /D %PARAFEM%
)

:end
