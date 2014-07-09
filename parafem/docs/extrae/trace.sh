#!/bin/sh

export EXTRAE_HOME=/gpfs/apps/NVIDIA/CEPBATOOLS/extrae/2.3.4/default/64
export EXTRAE_CONFIG_FILE=extrae.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so

$*
