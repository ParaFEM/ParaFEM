RT_PFEMVIEW=..\rt_pfemview
CPDIR=xcopy /I /S /Y /Q /E
RM=del /Q
PACK_MACHINE=avsenv v bin runtime data\helix data\helix_time pfemview.bat pfemview.ico

# Doesn't do the actual runtime compile but does some post installation
runtime:
	@echo "Copying distributable files to $(RT_PFEMVIEW)"
	@sed -e "s/XP_PATH.*/XP_PATH=./g" avsenv > $(RT_PFEMVIEW)\avsenv
	@sed -e "s/MACHINE=.*/MACHINE=$(MACHINE)/g" rt_pfemview.bat > $(RT_PFEMVIEW)\pfemview.bat
	@copy /Y pfemview.ico $(RT_PFEMVIEW)\pfemview.ico
	@if not exist $(RT_PFEMVIEW)\data $(CPDIR) data $(RT_PFEMVIEW)\data

package: runtime
	@echo "Building package $(RT_PFEMVIEW)\pfemview-$(MACHINE)-$(VER).zip"
        @cd ..\rt_pfemview
	@zip -r pfemview-$(MACHINE)-$(VER).zip $(PACK_MACHINE)

