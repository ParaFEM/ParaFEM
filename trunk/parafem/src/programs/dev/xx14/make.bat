@echo off

set DRIVER_NAME=p121
set PARAFEM=c:\parafem\parafem
set G95=mpif77 /integer-size:32 /real-size:64 /fast
set FFLAGS=-c -I%PARAFEM%\include
set LIB=%LIB%;%PARAFEM%\lib

IF "%1"=="clean" (
	del *.obj
	del *.exe
) else (
	echo "Building %DRIVER_NAME% RELEASE"
	%G95% %FFLAGS% *.f90
	%G95% *.obj libparafem.lib -o %DRIVER_NAME%.exe  
	copy *.exe %PARAFEM%\bin
)