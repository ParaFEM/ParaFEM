# Makefile for FELIB

INCLUDES=/I.
LIBDIR=..\lib
SHLIB = $(LIBDIR)\fe.dll
LIB = $(LIBDIR)\libfe.lib

CC      = cl /nologo
CXX     = cl /nologo
LINK    = link
AR      = lib /nologo

!ifdef XP_COMPILE_WITH_DEBUG
G = /Od /Zi
!else
G = /DNDEBUG /O2 /fp:fast 
!endif


!if "$(MACHINE)" != "pc_64"
CFLAGS =/MD $(G) /DWIN32 /D_WIN32 /D_CRT_SECURE_NO_DEPRECATE $(INCLUDES) 
!else
CFLAGS =/MD $(G) /DWIN32 /D_WIN32 /DWIN64 /D_WIN64 /D_CRT_SECURE_NO_DEPRECATE $(INCLUDES) 
!endif

CXXFLAGS = $(CFLAGS) /EHsc /DWINVER=0x0500 

# LDFLAGS = /debug
!if "$(MACHINE)" != "pc_64"
LDFLAGS = /machine:X86 /nodefaultlib:libc /nodefaultlib:libcmt /opt:noref
ARFLAGS = /machine:X86
!else
LDFLAGS = /machine:x64 /nodefaultlib:libc /nodefaultlib:libcmt /opt:noref
ARFLAGS = /machine:X64
!endif

OBJS = \
	felib.obj \
	fe_file_d.obj \
	fe_file_dis.obj

# Rule to make .obj from .cxx
.cxx.obj:
        $(CXX) $(CXXFLAGS) /Fo$*.obj /c $<

# Choose either the DLL (SHLIB) or static lib (LIB) library here.
all: $(LIBDIR) $(LIB)

$(SHLIB): $(OBJS)
	@echo Missing $(SHLIB) target in nt.mak. Cannot link. FIXME.

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) /out:$@ $(OBJS)

$(LIBDIR):
	-mkdir $(LIBDIR)

clean:
	del /f $(OBJS) $(SHLIB) $(LIB) *.pdb *.manifest *.ilk

