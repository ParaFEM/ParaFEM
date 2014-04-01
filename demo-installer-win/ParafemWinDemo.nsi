# ParaFEM Demonstrator for Windows installer script
# Requires Nullsoft Scriptable Install System (http://nsis.sourceforge.net/)
# J.S.Robinson@soton.ac.uk

OutFile "parafem-demo-installer.exe"

!include "MUI2.nsh"
!include "LogicLib.nsh"
!include "AdvReplaceInFile.nsh"

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
	
	#Icon
	File /r "parafem.ico"
	
	#Binaries
	CreateDirectory "$INSTDIR\bin"
	SetOutPath "$INSTDIR\bin"
	File /r "..\parafem\bin\*.exe"
	File "*.bat"
	
	#Examples
	CreateDirectory "$INSTDIR\examples"
	SetOutPath  "$INSTDIR\examples"
	File /r "..\parafem\examples\5th_ed\*"
	
	#Replace CASE file path in Paraview saved state files	
	#P121
	Push D:\parafem-data\p121\                 #text to be replaced
	Push $INSTDIR\examples\p121\demo\          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\examples\p121\demo\p121_small.pvsm     #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
	
	#P122
	Push D:\parafem-data\p122\dino_small\      #text to be replaced
	Push $INSTDIR\examples\p122\demo\          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\examples\p122\demo\p122_dino_small.pvsm     #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
	
	#P126
	Push D:\parafem-data\p126\tiny\			   #text to be replaced
	Push $INSTDIR\examples\p126\demo\          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\examples\p126\demo\p126_tiny.pvsm     #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
	
	#Replace @INSTDIR@ in driver batch files
	#demo-driver.bat
	Push @INSTDIR@                 #text to be replaced
	Push $INSTDIR          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\bin\demo-driver.bat     #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
	
	#p121.bat
	Push @INSTDIR@                 #text to be replaced
	Push $INSTDIR          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\bin\p121.bat     #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
	
	#p122.bat
	Push @INSTDIR@                 #text to be replaced
	Push $INSTDIR          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\bin\p122.bat     #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
	
	#p126.bat
	Push @INSTDIR@                 #text to be replaced
	Push $INSTDIR          #replace with
	Push all                                   #replace all occurrences
	Push all                                   #replace all occurrences
	Push $INSTDIR\bin\p126.bat     #file to replace in
	Call AdvReplaceInFile                      #call find and replace function
SectionEnd
	
# Start Menu Entries
RequestExecutionLevel user
Section	
	# ParaFEM demo-driver
	createDirectory "$SMPROGRAMS\${MUI_PRODUCT}"
	createShortCut "$SMPROGRAMS\${MUI_PRODUCT}\ParaFEM-Demos.lnk" "$INSTDIR\bin\demo-driver.bat" "" "$INSTDIR\parafem.ico"
	
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