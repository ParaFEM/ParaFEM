#$Id$
# BUILD INCLUDE build/template.inc
#
# SYSTEM: The University of Bristol Bluecrystal cluster
# NOTES: I'm only interested in building XX14 project.
#	There are lots of build errors for other projects,
#	but the ParaFEM library seems to build ok.
#
# --------------------------------------------------------------------------------
# Blank template for specific architectures. Complete the variables below as
# necessary. Use xe6.inc as an example.
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# Generic settings for compilation
# --------------------------------------------------------------------------------
FC=			mpiifort

# --------------------------------------------------------------------------------
# FLAGS for EXTERNAL LIBRARY compliation
# --------------------------------------------------------------------------------
EXT_LIB_ARPACK_FC=	mpiifort
EXT_LIB_ARPACK_FFLAGS=
EXT_LIB_ARPACK_LDFLAGS=

EXT_LIB_METIS_CC=	mpiicc
EXT_LIB_METIS_OPTFLAGS=
EXT_LIB_METIS_COPTIONS=
EXT_LIB_METIS_LDOPTIONS=
EXT_LIB_METIS_AR=	ar
EXT_LIB_METIS_RANLIB=	ranlib
  
# --------------------------------------------------------------------------------
# FLAGS for MODULE/LIBRARY compliation
# --------------------------------------------------------------------------------
# Form is LIB_(MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
# These flags are directly used by the module makefiles
# --------------------------------------------------------------------------------
BC3_FFLAGS=		-O2 -r8 -mt_mpi -warn all -stand f08
LIB_MPI_REL_FFLAGS=			$(BC3_FFLAGS)
LIB_MPI_DBG_FFLAGS=			$(BC3_FFLAGS)
  
LIB_SERIAL_REL_FFLAGS=			$(BC3_FFLAGS)
LIB_SERIAL_DBG_FFLAGS=			$(BC3_FFLAGS)
  
LIB_BEM_REL_FFLAGS=			$(BC3_FFLAGS)
LIB_BEM_DBG_FFLAGS=			$(BC3_FFLAGS)
  
LIB_BEM_P_REL_FFLAGS=			$(BC3_FFLAGS)
LIB_BEM_P_DBG_FFLAGS=			$(BC3_FFLAGS)

# --------------------------------------------------------------------------------
# DEFAULT FLAGS for all program groups
# --------------------------------------------------------------------------------
# Form is DEFAULT_(P12|XX|PROG8|GAITFEM)_(MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
#   where applicable.
# These are NOT directly used in any makefiles; specific program flags are assigned
# these defaults below and may be overridden.
# --------------------------------------------------------------------------------

# generic P12 group FFLAGS and LDFLAGS
DEFAULT_P12_MPI_REL_FFLAGS=		$(BC3_FFLAGS)
DEFAULT_P12_MPI_DBG_FFLAGS=
DEFAULT_P12_MPI_REL_LDFLAGS=
DEFAULT_P12_MPI_DBG_LDFLAGS=

# generic XX group FFLAGS and LDFLAGS
DEFAULT_XX_MPI_REL_FFLAGS=	$(BC3_FFLAGS)
DEFAULT_XX_MPI_DBG_FFLAGS=	$(BC3_FFLAGS)
DEFAULT_XX_MPI_REL_LDFLAGS=
DEFAULT_XX_MPI_DBG_LDFLAGS=

# generic PROG8 group FFLAGS and LDFLAGS
DEFAULT_PROG8_BEM_REL_FFLAGS=		$(BC3_FFLAGS)
DEFAULT_PROG8_BEM_DBG_FFLAGS=		$(BC3_FFLAGS)
DEFAULT_PROG8_BEM_REL_LDFLAGS=
DEFAULT_PROG8_BEM_DBG_LDFLAGS=
DEFAULT_PROG8_BEM_P_REL_FFLAGS=		$(BC3_FFLAGS)
DEFAULT_PROG8_BEM_P_DBG_FFLAGS=		$(BC3_FFLAGS)
DEFAULT_PROG8_BEM_P_REL_LDFLAGS=	
DEFAULT_PROG8_BEM_P_DBG_LDFLAGS=

# generic gaitfem group FFLAGS and LDFLAGS
DEFAULT_GAITFEM_SERIAL_REL_FFLAGS=	$(BC3_FFLAGS)
DEFAULT_GAITFEM_SERIAL_DBG_FFLAGS=	$(BC3_FFLAGS)
DEFAULT_GAITFEM_SERIAL_REL_LDFLAGS=
DEFAULT_GAITFEM_SERIAL_DBG_LDFLAGS=
  
# --------------------------------------------------------------------------------
# DEFAULT FLAGS for all tool groups
# --------------------------------------------------------------------------------
# Form is DEFAULT_(PREPROC|POSTPROC)_(FC|CC|MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
#   where applicable.
# These are NOT directly used in any makefiles; specific tool flags are assigned
# these defaults below and may be overridden.
# --------------------------------------------------------------------------------

# generic PREPROC group FFLAGS and LDFLAGS
DEFAULT_PREPROC_MPI_REL_FFLAGS=		$(BC3_FFLAGS)
DEFAULT_PREPROC_MPI_DBG_FFLAGS=		$(BC3_FFLAGS)
DEFAULT_PREPROC_MPI_REL_LDFLAGS=
DEFAULT_PREPROC_MPI_DBG_LDFLAGS=

DEFAULT_PREPROC_CC_REL_FFLAGS=		-O2
DEFAULT_PREPROC_CC_DBG_FFLAGS=
DEFAULT_PREPROC_FC_REL_FFLAGS=		$(BC3_FFLAGS)
DEFAULT_PREPROC_FC_DBG_FFLAGS=

# generic POSTPROC group FFLAGS and LDFLAGS
DEFAULT_POSTPROC_MPI_REL_FFLAGS=	$(BC3_FFLAGS)
DEFAULT_POSTPROC_MPI_DBG_FFLAGS=	$(BC3_FFLAGS)
DEFAULT_POSTPROC_MPI_REL_LDFLAGS=
DEFAULT_POSTPROC_MPI_DBG_LDFLAGS=

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# SPECIFIC FLAGS for PROGRAMS and TOOLS
# By default these inherit the defaults defined above.
# NOTE: If a specific program needs different flags then modify below.
# NOTE: It may be useful to create a new include file for a modified architecture.
#   e.g., <arch>_test.inc and use MACHINE=<arch>_test in the build script.
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# PROGRAM FLAGS for all programs.
# --------------------------------------------------------------------------------
# All use their program group defaults defined above unless overridden below
# --------------------------------------------------------------------------------
P121_MPI_REL_FFLAGS=$(DEFAULT_P12_MPI_REL_FFLAGS)
P121_MPI_REL_LDFLAGS=$(DEFAULT_P12_MPI_REL_LDFLAGS)
P121_MPI_DBG_FFLAGS=$(DEFAULT_P12_MPI_DBG_FFLAGS)
P121_MPI_DBG_LDFLAGS=$(DEFAULT_P12_MPI_DBG_LDFLAGS)

P126_MPI_REL_FFLAGS=$(DEFAULT_P12_MPI_REL_FFLAGS)
P126_MPI_REL_LDFLAGS=$(DEFAULT_P12_MPI_REL_LDFLAGS)
P126_MPI_DBG_FFLAGS=$(DEFAULT_P12_MPI_DBG_FFLAGS)
P126_MPI_DBG_LDFLAGS=$(DEFAULT_P12_MPI_DBG_LDFLAGS)
  
P129_MPI_REL_FFLAGS=$(DEFAULT_P12_MPI_REL_FFLAGS)
P129_MPI_REL_LDFLAGS=$(DEFAULT_P12_MPI_REL_LDFLAGS)
P129_MPI_DBG_FFLAGS=$(DEFAULT_P12_MPI_DBG_FFLAGS)
P129_MPI_DBG_LDFLAGS=$(DEFAULT_P12_MPI_DBG_LDFLAGS)
  
P1210_MPI_REL_FFLAGS=$(DEFAULT_P12_MPI_REL_FFLAGS)
P1210_MPI_REL_LDFLAGS=$(DEFAULT_P12_MPI_REL_LDFLAGS)
P1210_MPI_DBG_FFLAGS=$(DEFAULT_P12_MPI_DBG_FFLAGS)
P1210_MPI_DBG_LDFLAGS=$(DEFAULT_P12_MPI_DBG_LDFLAGS)
  
PROG81_BEM_REL_FFLAGS=$(DEFAULT_PROG8_BEM_REL_FFLAGS)
PROG81_BEM_REL_LDFLAGS=$(DEFAULT_PROG8_BEM_REL_LDFLAGS)
PROG81_BEM_DBG_FFLAGS=$(DEFAULT_PROG8_BEM_DBG_FFLAGS)
PROG81_BEM_DBG_LDFLAGS=$(DEFAULT_PROG8_BEM_DBG_LDFLAGS)
  
PROG82_BEM_REL_FFLAGS=$(DEFAULT_PROG8_BEM_REL_FFLAGS)
PROG82_BEM_REL_LDFLAGS=$(DEFAULT_PROG8_BEM_REL_LDFLAGS)
PROG82_BEM_DBG_FFLAGS=$(DEFAULT_PROG8_BEM_DBG_FFLAGS)
PROG82_BEM_DBG_LDFLAGS=$(DEFAULT_PROG8_BEM_DBG_LDFLAGS)
  
PROG83_BEM_P_REL_FFLAGS=$(DEFAULT_PROG8_BEM_P_REL_FFLAGS)
PROG83_BEM_P_REL_LDFLAGS=$(DEFAULT_PROG8_BEM_P_REL_LDFLAGS)
PROG83_BEM_P_DBG_FFLAGS=$(DEFAULT_PROG8_BEM_P_DBG_FFLAGS)
PROG83_BEM_P_DBG_LDFLAGS=$(DEFAULT_PROG8_BEM_P_DBG_LDFLAGS)
  
GAITFEM_SERIAL_REL_FFLAGS=$(DEFAULT_GAITFEM_SERIAL_REL_FFLAGS)
GAITFEM_SERIAL_REL_LDFLAGS=$(DEFAULT_GAITFEM_SERIAL_REL_LDFLAGS)
GAITFEM_SERIAL_DBG_FFLAGS=$(DEFAULT_GAITFEM_SERIAL_DBG_FFLAGS)
GAITFEM_SERIAL_DBG_LDFLAGS=$(DEFAULT_GAITFEM_SERIAL_DBG_LDFLAGS)
  
XX1_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX1_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX1_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX1_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

XX2_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX2_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX2_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX2_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

XX3_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX3_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX3_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX3_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

XX4_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX4_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX4_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX4_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

XX5_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX5_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX5_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX5_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

XX6_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX6_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX6_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX6_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

XX7_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX7_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX7_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX7_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

XX8_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX8_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX8_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX8_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

XX9_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX9_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX9_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX9_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

# The use of a coarray config file is required for
# execution across node boundaries, i.e. always!
# The file name should be relative to the directory
# from where xx14 is launched.
# Here I assume it is the same directory where it is built.
# The file name *must* be the same as used for building
# CGPACK. 
XX14_MPI_REL_FFLAGS=	$(DEFAULT_XX_MPI_REL_FFLAGS) -I$(HOME)/mod \
			-coarray=distributed -coarray-config-file=xx14.conf
XX14_MPI_REL_LDFLAGS=	$(DEFAULT_XX_MPI_REL_LDFLAGS) -mt_mpi -I$(HOME)/mod -L$(HOME)/lib -lcgpack \
			-coarray=distributed -coarray-config-file=xx14.conf
#			-coarray=single
XX14_MPI_DBG_FFLAGS=	$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX14_MPI_DBG_LDFLAGS=	$(DEFAULT_XX_MPI_DBG_LDFLAGS)

# --------------------------------------------------------------------------------
# TOOL FLAGS for selected tools.
# --------------------------------------------------------------------------------
# All use their tool group defaults defined above unless overridden below
# --------------------------------------------------------------------------------

# Pre Processing tools
# Note: There are no flags for INP2D
MG2D_MPI_REL_FFLAGS=$(DEFAULT_PREPROC_MPI_REL_FFLAGS)
MG2D_MPI_REL_LDFLAGS=$(DEFAULT_PREPROC_MPI_REL_LDFLAGS)
MG2D_MPI_DBG_FFLAGS=$(DEFAULT_PREPROC_MPI_DBG_FFLAGS)
MG2D_MPI_DBG_LDFLAGS=$(DEFAULT_PREPROC_MPI_DBG_LDFLAGS)

PARTITIONER_REL_FFLAGS=$(DEFAULT_PREPROC_CC_REL_FFLAGS)
PARTITIONER_DBG_FFLAGS=$(DEFAULT_PREPROC_CC_DBG_FFLAGS)

DINOMESH_REL_FFLAGS=$(DEFAULT_PREPROC_FC_REL_FFLAGS)
DINOMESH_DBG_FFLAGS=$(DEFAULT_PREPROC_FC_DBG_FFLAGS)
  
FEARFY_REL_FFLAGS=$(DEFAULT_PREPROC_FC_REL_FFLAGS)
FEARFY_DBG_FFLAGS=$(DEFAULT_PREPROC_FC_DBG_FFLAGS)
  
# Post Processing tools
BND2DIS_MPI_REL_FFLAGS=$(DEFAULT_POSTPROC_MPI_REL_FFLAGS)
BND2DIS_MPI_REL_LDFLAGS=$(DEFAULT_POSTPROC_MPI_REL_LDFLAGS)
BND2DIS_MPI_DBG_FFLAGS=$(DEFAULT_POSTPROC_MPI_DBG_FFLAGS)
BND2DIS_MPI_DBG_LDFLAGS=$(DEFAULT_POSTPROC_MPI_DBG_LDFLAGS)

D2OFF_MPI_REL_FFLAGS=$(DEFAULT_POSTPROC_MPI_REL_FFLAGS)
D2OFF_MPI_REL_LDFLAGS=$(DEFAULT_POSTPROC_MPI_REL_LDFLAGS)
D2OFF_MPI_DBG_FFLAGS=$(DEFAULT_POSTPROC_MPI_DBG_FFLAGS)
D2OFF_MPI_DBG_LDFLAGS=$(DEFAULT_POSTPROC_MPI_DBG_LDFLAGS)

MAT2NUMPE_MPI_REL_FFLAGS=$(DEFAULT_POSTPROC_MPI_REL_FFLAGS)
MAT2NUMPE_MPI_REL_LDFLAGS=$(DEFAULT_POSTPROC_MPI_REL_LDFLAGS)
MAT2NUMPE_MPI_DBG_FFLAGS=$(DEFAULT_POSTPROC_MPI_DBG_FFLAGS)
MAT2NUMPE_MPI_DBG_LDFLAGS=$(DEFAULT_POSTPROC_MPI_DBG_LDFLAGS)

UPVW2DIS_MPI_REL_FFLAGS=$(DEFAULT_POSTPROC_MPI_REL_FFLAGS)
UPVW2DIS_MPI_REL_LDFLAGS=$(DEFAULT_POSTPROC_MPI_REL_LDFLAGS)
UPVW2DIS_MPI_DBG_FFLAGS=$(DEFAULT_POSTPROC_MPI_DBG_FFLAGS)
UPVW2DIS_MPI_DBG_LDFLAGS=$(DEFAULT_POSTPROC_MPI_DBG_LDFLAGS)
