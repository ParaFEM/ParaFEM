# ParaFEM Demonstrator for Windows installer script
# Requires Nullsoft Scriptable Install System 3.x (http://nsis.sourceforge.net/)
# J.S.Robinson@soton.ac.uk

OutFile "parafem-demo-installer.exe"

!include "MUI2.nsh"
!include "LogicLib.nsh"
!include "include\AdvReplaceInFile.nsh"

Name "ParaFEM Demonstrator"
!define MUI_PRODUCT "ParaFEM_Demo"
!define MUI_ICON "parafem.ico"

#InstallDir "$PROGRAMFILES32\${MUI_PRODUCT}"
InstallDir "C:\${MUI_PRODUCT}"
   
!define MUI_PAGE_HEADER_TEXT "ParaFEM Demonstrator"
!define MUI_WELCOMEPAGE_TITLE "ParaFEM Demonstrator" 
!define MUI_WELCOMEPAGE_TEXT "Welcome to the ParaFEM Demonstrator installer.  This package will install a selection of programs that demonstrate the power of ParaFEM, the freely available portable library for parallel finite element analysis."
!define MUI_WELCOMEFINISHPAGE_BITMAP "parafem_logo_full.bmp"
!insertmacro MUI_PAGE_WELCOME

!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES

!define MUI_FINISHPAGE_LINK "For more information visit http://parafem.org.uk"
!define MUI_FINISHPAGE_LINK_LOCATION "http://parafem.org.uk"
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_LANGUAGE "English"

# default section start
Section "install"

	#OpenMPI
	SetOutPath "$INSTDIR"
	File /r ".\vendor\OpenMPI_v1.6.2-win32"
	#FIXME - This a hack to get OpenMPI binaries to find their resources
	File /r  ".\vendor\OpenMPI_v1.6.2-win32\share"
	
	#ParaView
	SetOutPath "$INSTDIR"
	File /r ".\vendor\ParaView 4.1.0"
	
	#Icons & images
	File /r "parafem.ico"
	File /r "parafem_logo_full.bmp"
	
	#Launcher
	File /r "demo-driver.hta"
	
	#Binaries
	CreateDirectory "$INSTDIR\bin"
	SetOutPath "$INSTDIR\bin"
	File /r "..\parafem\bin\*.exe"
	File "*.bat"
	
	#Examples
	CreateDirectory "$INSTDIR\examples"
	SetOutPath  "$INSTDIR\examples"
	File /r "..\parafem\examples\5th_ed\*"
	
	#Replace @INSTDIR@ in driver batch files
	#demo-driver.bat
	Push @INSTDIR@                 #text to be replaced
	Push $INSTDIR          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\bin\demo-driver.bat    #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
	
	#demo-driver.hta
	Push @INSTDIR@                 #text to be replaced
	Push $INSTDIR          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\demo-driver.hta     #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
	
	
SectionEnd
	
# Start Menu Entries
RequestExecutionLevel user
Section	
	# ParaFEM demo-driver
	createDirectory "$SMPROGRAMS\${MUI_PRODUCT}"
	createShortCut "$SMPROGRAMS\${MUI_PRODUCT}\ParaFEM-Demos.lnk" "$SYSDIR\mshta.exe" "$INSTDIR\demo-driver.hta" "$INSTDIR\parafem.ico"
	
	# ParaView
	createShortCut "$SMPROGRAMS\${MUI_PRODUCT}\ParaView.lnk" "$INSTDIR\ParaView 4.1.0\bin\paraview.exe"
		
	#Uninstaller
	writeUninstaller "$INSTDIR\uninstall.exe"
	createShortCut "$SMPROGRAMS\${MUI_PRODUCT}\Uninstall.lnk" "$INSTDIR\uninstall.exe"	
SectionEnd

Section "uninstall"
	# Clean up Start Menu
	delete "$SMPROGRAMS\${MUI_PRODUCT}\ParaFEM-Demos.lnk"
	delete "$SMPROGRAMS\${MUI_PRODUCT}\ParaView.lnk"
	delete "$SMPROGRAMS\${MUI_PRODUCT}\Uninstall.lnk"
	rmDir "$SMPROGRAMS\${MUI_PRODUCT}"
	
	# Remove files
	rmDir /r "$INSTDIR"
SectionEnd