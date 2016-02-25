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
	
	C:\parafem\parafem\README-parafem-and-drivers-win32-intel.txt

2. Install the Nullsoft Scriptable Install System 3.x [1]

Build Procedure
---------------

1. Open the NSIS script compiler (MakeNSISW) then load and run the build script:
	C:\parafem\demo-installer-win\ParafemWinDemo.nsi

Note that for to add / change the driver programs which are included in the 
demonstrator, the following changes must be made:

1. Edit the build script: 

C:\parafem\parafem\make-parafem-and-drivers-win32-intel.bat

Set the variable BUILD_GROUP_5ED to include the drivers you want built

2. Add a directory for input data, e.g:

 C:\parafem\parafem\examples\5th_ed\pXXX\demo
 
Ensure that the input files are named in the format pXXX_demo.x
Edit the Paraview saved state file pXXX_demo.psvm, so that the path to the 
CASE file is replaced by the placeholder %CASE_FILE_PATH% (in order that the
installer can replace this with the real path on the target system). eg.:

<Property name="CaseFileName" id="2148.CaseFileName" number_of_elements="1">
    <Element index="0" value="%CASE_FILE_PATH%"/>
    <Domain name="files" id="2148.CaseFileName.files"/>
</Property>



References
----------

[1] http://nsis.sourceforge.net/
