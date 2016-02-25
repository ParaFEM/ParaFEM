#!/bin/bash
# @ partition        = debug
# @ class            = debug
# @ job_name         = parafem_p121
# @ initialdir       = .
# @ output           = parafem_%j.out
# @ error            = parafem_%j.err
# @ total_tasks      = 1
## @ gpus_per_node   = 2
# @ tracing          = 1
# @ cpus_per_task    = 1
# @ wall_clock_limit = 00:02:00

# Set EXE to your application binary
export EXE=$HOME/parafem/parafem/bin/p121_5

export TRACENAME=parafem_p121_trace.prv

# Run the application
time srun ./trace.sh $EXE p121_tiny

