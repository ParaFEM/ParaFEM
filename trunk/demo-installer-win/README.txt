===============================================================================
Building ParaFEM and 5th Edition example driver programs binary installer on 
Microsoft Windows
J.Robinson@software.ac.uk
===============================================================================

OS Prerequisites
----------------

Microsoft Windows XP or newer (Tested on Windows 7, 8 and 2008 Server).
	
Software Prerequisites
----------------------

(This assumes you have checked out the ParaFEM Subversion trunk to C:\parafem)

1. Build ParaFEM and the driver executables from source.  For details see:
	C:\parafem\README-parafem-and-drivers-win32-intel.txt

2. Install the Nullsoft Scriptable Install System 3.x [1]

Build Procedure
---------------

1. Open the NSIS script compiler (MakeNSISW) then load and run the build script:
	C:\parafem\demo-installer-win\ParafemWinDemo.nsi

Note that for to add / change the driver programs which are included in the 
demonstrator, the following changes must be made:

1. A driver (pXXX.bat) file should be created in:
	C:\parafem\demo-installer-win\
This should invoke the driver executable with paths the appropriate input files
It should then invoke the bundled paraview with the appropriate saved state.
See p121.bat for an exemplar.

2. The file: 
	C:\parafem\demo-installer-win\demo-driver.bat 
should be updated to describe the new executable, add a menu option to select 
it, and a code block to invoke the pXXX.bat file

3. The NSIS build script:
	C:\parafem\demo-installer-win\ParafemWinDemo.nsi
 should be updated to replace CASE file path in Paraview saved state files for 
 the new driver program.


References
----------

[1] http://nsis.sourceforge.net/
