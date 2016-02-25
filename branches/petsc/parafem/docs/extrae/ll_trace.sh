# @ job_name = trace
# @ executable = /bgsys/drivers/ppcfloor/hlcs/bin/runjob
# @ arguments =--envs "EXTRAE_CONFIG_FILE=extrae.xml" --exe ./xx1 --np 16 -p 8 --args "p121_medium"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL;
# @ wall_clock_limit = 00:40:00
# @ job_type = bluegene
# @ bg_size = 128
# @ class = prod
# @ queue

# Enable these if you don't have libxml2 support in Extrae
# export EXTRAE_ON=1
# export EXTRAE_MPI_COUNTERS_ON=1
# export EXTRAE_COUNTERS=PAPI_TOT_INS,PAPI_TOT_CYC,PAPI_L1_LDM,PAPI_BR_MSP,PAPI_FP_INS,PAPI_TLB_DM

# Enable this if you have libxml2 support in Extrae
#export EXTRAE_CONFIG_FILE=extrae.xml
export EXTRAE_CONFIG_FILE=extrae.xml
