# Makefile for libfe 

INCLUDES=-I.
LIBDIR=../lib
SHLIB = $(LIBDIR)/libfe.so
LIB   = $(LIBDIR)/libfe.a

OPT=-O
ARCH=-fPIC
SHAREDFLAG=-shared
AR=ar
ARFLAGS=crl
CFLAGS=$(ARCH) $(OPT) $(DEFINES) $(INCLUDES)
LDFLAGS=$(ARCH)

LIBSRC= \
	felib.cxx \
	fe_file_d.cxx \
	fe_file_dis.cxx

LIBHDRS=felib.h fe_file_d.h fe_file_dis.h

LIBOBJ=${LIBSRC:.cxx=.o}

.SUFFIXES: .o .cxx

.cxx.o:
	$(CXX) -c $< $(CFLAGS)

# Choose either the shared (SHLIB) or static (LIB) library here.
all: $(LIBDIR) $(LIB)

$(LIBDIR):
	-mkdir $(LIBDIR)

$(SHLIB): $(LIBOBJ) $(LIBHDRS)
	$(CXX) $(LDFLAGS) $(SHAREDFLAG) -o $@ $(LIBOBJ)

$(LIB): $(LIBOBJ) $(LIBHDRS)
	$(AR) $(ARFLAGS) $@ $(LIBOBJ)

clean::
	-rm -f $(SHLIB) $(LIB) $(LIBOBJ)

