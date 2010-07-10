RT_PFEMVIEW=../rt_pfemview
SHELL=sh
RM=rm -r
RT_LIBS=libgdprt.so libag.so libagx.so libgdogl.so libgdsw.so libgd.so \
	libvxbx.so libui.so libmods.so libpal.so libdmap.so libdbxx.so \
	libgmod.so libmutil.so libfld.so libgdif.so libwtobj.so libwtool.so \
	libgdvps.so libgdvrm.so libcxxut.so libomx.so libom.so libcomm.so \
	libtool.so libxplic.so
PACK_MACHINE=avsenv v bin runtime lib/$(MACHINE)/*.so* data/helix data/helix_time pfemview

# Why doesn't Express add the .so to the compiled project?
addlibs:
	@echo "Copying shared libraries to $(RT_PFEMVIEW)/lib/$(MACHINE)"
	-@mkdir -p $(RT_PFEMVIEW)/lib/$(MACHINE)
	@for l in ${RT_LIBS} ; do cp $(XP_ROOT)/lib/$(MACHINE)/$$l $(RT_PFEMVIEW)/lib/$(MACHINE); done
	@cd $(RT_PFEMVIEW)/lib/$(MACHINE) ; for l in ${RT_LIBS} ; do ln -fs $$l $$l.7.2.1; done

# Doesn't do the actual runtime compile but does some post installation
runtime: addlibs
	@echo "Copying distributable files to $(RT_PFEMVIEW)"
	@sed -e "s/XP_PATH.*/XP_PATH=./g" avsenv > $(RT_PFEMVIEW)/avsenv
	@sed -e "s/MACHINE=.*/MACHINE=$(MACHINE)/g" rt_pfemview > $(RT_PFEMVIEW)/pfemview
	@chmod a+x $(RT_PFEMVIEW)/pfemview
	@if [ -e $(RT_PFEMVIEW)/data ] ; then $(RM) $(RT_PFEMVIEW)/data; fi
	@cp -r data $(RT_PFEMVIEW)/data

package: runtime
	@echo "Building package $(RT_PFEMVIEW)/pfemview-$(MACHINE)-$(VER).tgz"
	@cd $(RT_PFEMVIEW) ; tar zcf pfemview-$(MACHINE)-$(VER).tgz $(PACK_MACHINE) 
