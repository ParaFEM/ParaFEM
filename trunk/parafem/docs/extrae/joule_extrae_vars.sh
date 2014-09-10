#!/bin/bash

export EXTRAE_HOME='/gpfs/home/HCP010/lxm03/rxm73-lxm03/extrae-2.5.1-linux-bgq-papi-5.0.1'
export PAPI_HOME='/gpfs/packages/ibm/papi/5.1.0'
export PAPI_LIBS='-L$(PAPI_HOME)/lib -lpapi -L/bgsys/drivers/ppcfloor/spi/lib -lSPI -lSPI_cnk -lstdc++ -lrt'
export XML2_HOME='/gpfs/home/HCP010/lxm03/rxm73-lxm03/libxml2-2.9.1'
export XML2_LIBS='-lxml2'
export BFD_LIBS='-L/gpfs/home/HCP010/lxm03/rxm73-lxm03/binutils/binutils-2.24/bfd -lbfd -L/gpfs/home/HCP010/lxm03/rxm73-lxm03/binutils/binutils-2.24/libiberty -liberty'
export ZLIB_HOME='/gpfs/home/HCP010/lxm03/rxm73-lxm03/zlib-1.2.8'
export ZLIB_LIBS='-L/gpfs/home/HCP010/lxm03/rxm73-lxm03/zlib-1.2.8/lib -lz'
