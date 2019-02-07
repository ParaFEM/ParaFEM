# BUILD INCLUDE build/csf-gpu.inc
# Author: George Leaver, Louise M. Lever
# --------------------------------------------------------------------------------
# Include file for the CSF GPU Nodes (Nvidia) architecture.
# --------------------------------------------------------------------------------
# Source changes:
# src/tools/postprocessing/endian_convert/endian_convert.f90 Line 50.
# The FORTRAN OPEN() function uses CONVERT='SWAP' but SWAP is not valid with ifort.
# Instead we use NATIVE and the user can set the env var FORT_CONVERT10=BIG_ENDIAN
# or LITTLE_ENDIAN so select the correct format.
#
# This MACHINE file:
# The CSF2 has Nvidia K20 GPUs hence use nvcc -arch=sm_35 
# The CSF3 has Nvidia v100 GPUs hence use nvcc -arch=sm_70
#
# I have remove -fast from the ifort command. It switches on static linking and we
# do not have static OpenMPI libraries. -fast also trades floating point accuracy
# for speed so the other settings enabled by -fast:
# -ipo, -O3, -no-prec-div, -static, -fp-model fast=2
# will not be used. That's OK - I prefer more accurate results.
#
# The compiler will compile routines for automatic CPU dispatch - several versions
# of routines will be compiled, optimized for different Intel CPUs. At runtime the
# executable will automatically run the version suitable for the CPU on which the
# job is running. You still only have one executable.
#
# GPU demo programs:
# The following src/programs/dev/ programs are GPU enabled:
# CUDA:    xx3 xx9
# OpenCL:  xx4 xx10
#
# To build on CSF2 (build CPU version then add in the GPU demos, for example):
# 
# module load mpi/intel-17.0/openmpi/1.8.3m
# module load libs/cuda/9.0.176
# ./make-parafem MACHINE=csf-gpu 2>&1 | tee make.csf-gpu.log
# ./make-parafem MACHINE=csf-gpu -gpu >&1 | tee make.csf-gpu.gpu.log
#
# To build on CSF3:
# module load mpi/intel-17.0/openmpi/3.1.1
# module load libs/cuda/9.2.148
# ./make-parafem MACHINE=csf-gpu 2>&1 | tee make.csf-gpu.log
# ./make-parafem MACHINE=csf-gpu -gpu 2>&1 | tee make.csf-gpu.gpu.log
#
# george.leaver@manchester.ac.uk - Research Infrastructure
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# Generic settings for compilation
# --------------------------------------------------------------------------------
FC=mpif90
CC=mpicc
INSTALL_CMD=mv
# ********* CHANGE TO CORRECT SETTING FOR GPU BUILD *********
# CSF2
NVSM=sm_35
# CSF3
# NVSM=sm_70

# --------------------------------------------------------------------------------
# FLAGS for EXTERNAL LIBRARY compliation
# --------------------------------------------------------------------------------
# CSF only has MPI libs for Intel and GNU compiler. Ignore PGI here.
PGI_EXT_FLAGS=-V -fastsse -r8
GNU_EXT_FLAGS=-fdefault-real-8
# GWL: Remove -fast for CSF compilation (we have no static OpenMPI libs)
INTEL_EXT_FLAGS=-r8 -msse4.2 -axCORE-AVX512,CORE-AVX2,AVX

EXT_FLAGS=$(INTEL_EXT_FLAGS)

EXT_LIB_ARPACK_FC=mpif90
EXT_LIB_ARPACK_FFLAGS=$(EXT_FLAGS)
EXT_LIB_ARPACK_LDFLAGS=

EXT_LIB_GAF77_FC=ifort
EXT_LIB_GAF77_FFLAGS=$(EXT_FLAGS)
EXT_LIB_GAF77_LDFLAGS=

EXT_LIB_METIS_CC=icc
EXT_LIB_METIS_OPTFLAGS=-O2 -msse4.2 -axCORE-AVX512,CORE-AVX2,AVX
EXT_LIB_METIS_COPTIONS=
EXT_LIB_METIS_LDOPTIONS=
EXT_LIB_METIS_AR=ar rv
EXT_LIB_METIS_RANLIB=ranlib

# --------------------------------------------------------------------------------
# FLAGS for GPU compliation
# --------------------------------------------------------------------------------

# Nvidia OpenCL - On CSF the CUDA_HOME var is set when the Nvidia module is loaded
GPU_OCL_NV_ROOT  = $(CUDA_HOME)
GPU_OCL_NV_INC   = -I$(GPU_OCL_NV_ROOT)/include 
GPU_OCL_NV_LIB   = -lOpenCL -lm
GPU_OCL_NV_FLAGS = 

# AMD OpenCL - No AMD h/w in CSF. Ignore.
GPU_OCL_AMD_ROOT     = /export/ati-app-sdk/lnx64/2.5
GPU_OCL_AMD_FLAGS    = -DAMDBLAS
GPU_OCL_AMD_BLASROOT = $(HOME)/clAmdBlas-1.4.182
GPU_OCL_AMD_BLASLIB  = -L$(GPU_OCL_AMD_BLASROOT)/lib64 -lclAmdBlas
GPU_OCL_AMD_INC      = -I$(GPU_OCL_AMD_ROOT)/include -I$(GPU_OCL_AMD_BLASROOT)/include
GPU_OCL_AMD_LIB      = -L$(GPU_OCL_AMD_ROOT)/lib/x86_64/ -lOpenCL $(GPU_OCL_AMD_BLASLIB) -lm

# Select OpenCL vendor to be used by program makefiles
GPU_OCL_ROOT  = $(GPU_OCL_NV_ROOT)
GPU_OCL_FLAGS = $(GPU_OCL_NV_FLAGS)
GPU_OCL_INC   = $(GPU_OCL_NV_INC)
GPU_OCL_LIB   = $(GPU_OCL_NV_LIB)

# Nvidia CUDA
GPU_CU_ROOT  = $(CUDA_HOME)
GPU_CU       = nvcc
GPU_CU_FLAGS = -arch=$(NVSM) --ptxas-options=-v
GPU_CU_INC   = -I$(GPU_CU_ROOT)/include
GPU_CU_LIB   = -L$(GPU_CU_ROOT)/lib64 -lcudart

# --------------------------------------------------------------------------------
# FLAGS for MODULE/LIBRARY compliation
# --------------------------------------------------------------------------------
# Form is LIB_(MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
# These flags are directly used by the module makefiles
# --------------------------------------------------------------------------------

# CSF only has MPI libs for Intel and GNU compiler. Ignore PGI here.
PGI_REL=-fastsse -Mlarge_arrays
PGI_DBG=-g -Mbounds -traceback
GNU_REL=
GNU_DBG=-g
# GWL: Remove -fast for CSF compilation (we have no static OpenMPI libs)
INTEL_REL=-msse4.2 -axCORE-AVX512,CORE-AVX2,AVX
INTEL_DBG=-g -traceback -msse4.2 -axCORE-AVX512,CORE-AVX2,AVX

REL_FLAGS=$(INTEL_REL)
DBG_FLAGS=$(INTEL_DBG)

LIB_MPI_REL_FFLAGS=$(REL_FLAGS)
LIB_MPI_DBG_FFLAGS=$(DBG_FLAGS)
 
LIB_SERIAL_REL_FFLAGS=$(REL_FLAGS)
LIB_SERIAL_DBG_FFLAGS=$(DBG_FLAGS)
  
LIB_BEM_REL_FFLAGS=$(REL_FLAGS)
LIB_BEM_DBG_FFLAGS=$(DBG_FLAGS)
  
LIB_BEM_P_REL_FFLAGS=$(REL_FLAGS)
LIB_BEM_P_DBG_FFLAGS=$(DBG_FLAGS)

# --------------------------------------------------------------------------------
# DEFAULT FLAGS for all program groups
# --------------------------------------------------------------------------------
# Form is DEFAULT_(P12|XX|PROG8|GAITFEM)_(MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
#   where applicable.
# These are NOT directly used in any makefiles; specific program flags are assigned
# these defaults below and may be overridden.
# --------------------------------------------------------------------------------

# generic P12 group FFLAGS and LDFLAGS
DEFAULT_P12_MPI_REL_LDFLAGS=$(REL_FLAGS)
DEFAULT_P12_MPI_DBG_LDFLAGS=$(DBG_FLAGS)

# generic XX group FFLAGS and LDFLAGS
DEFAULT_XX_MPI_REL_LDFLAGS=$(REL_FLAGS)
DEFAULT_XX_MPI_DBG_LDFLAGS=$(DBG_FLAGS)

# generic PROG8 group FFLAGS and LDFLAGS
DEFAULT_PROG8_BEM_REL_FFLAGS=
DEFAULT_PROG8_BEM_DBG_FFLAGS=
DEFAULT_PROG8_BEM_REL_LDFLAGS=
DEFAULT_PROG8_BEM_DBG_LDFLAGS=
DEFAULT_PROG8_BEM_P_REL_FFLAGS=
DEFAULT_PROG8_BEM_P_DBG_FFLAGS=
DEFAULT_PROG8_BEM_P_REL_LDFLAGS=
DEFAULT_PROG8_BEM_P_DBG_LDFLAGS=

# generic gaitfem group FFLAGS and LDFLAGS
DEFAULT_GAITFEM_SERIAL_REL_LDFLAGS=$(REL_FLAGS)
DEFAULT_GAITFEM_SERIAL_DBG_LDFLAGS=$(DBG_FLAGS)
  
# --------------------------------------------------------------------------------
# DEFAULT FLAGS for all tool groups
# --------------------------------------------------------------------------------
# Form is DEFAULT_(PREPROC|POSTPROC)_(FC|CC|MPI|SERIAL|OPENMP|GPU)_(REL|DBG)_(FFLAGS|LDFLAGS)
#   where applicable.
# These are NOT directly used in any makefiles; specific tool flags are assigned
# these defaults below and may be overridden.
# --------------------------------------------------------------------------------

# generic PREPROC group FFLAGS and LDFLAGS
DEFAULT_PREPROC_MPI_REL_LDFLAGS=$(REL_FLAGS)
DEFAULT_PREPROC_MPI_DBG_LDFLAGS=$(DBG_FLAGS)

DEFAULT_PREPROC_CC_REL_FFLAGS=-c -O2
DEFAULT_PREPROC_CC_DBG_FFLAGS=-c -g
#DEFAULT_PREPROC_FC_REL_FFLAGS=$(REL_FLAGS)
#DEFAULT_PREPROC_FC_DBG_FFLAGS=$(DBG_FLAGS)

# generic POSTPROC group FFLAGS and LDFLAGS
DEFAULT_POSTPROC_MPI_REL_LDFLAGS=$(REL_FLAGS)
DEFAULT_POSTPROC_MPI_DBG_LDFLAGS=$(DBG_FLAGS)

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# SPECIFIC FLAGS for PROGRAMS and TOOLS
# By default these inherit the defaults defined above.
# NOTE: If a specific program needs different flags then modify below.
# NOTE: It may be useful to create a new include file for a modified architecture.
#   e.g., xe6_test.inc and use MACHINE=xe6_test in the build script.
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

XX10_MPI_REL_FFLAGS=$(DEFAULT_XX_MPI_REL_FFLAGS)
XX10_MPI_REL_LDFLAGS=$(DEFAULT_XX_MPI_REL_LDFLAGS)
XX10_MPI_DBG_FFLAGS=$(DEFAULT_XX_MPI_DBG_FFLAGS)
XX10_MPI_DBG_LDFLAGS=$(DEFAULT_XX_MPI_DBG_LDFLAGS)

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

