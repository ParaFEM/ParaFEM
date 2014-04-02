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

Open the NSIS script compiler (MakeNSISW) then load and run the build script:
C:\parafem\demo-installer-win\ParafemWinDemo.nsi

References
----------

[1] http://nsis.sourceforge.net/
