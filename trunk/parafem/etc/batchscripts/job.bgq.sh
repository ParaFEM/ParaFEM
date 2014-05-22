#@bg_size=1
##@executable=/bgsys/drivers/ppcfloor/hlcs/bin/runjob
#@executable=job.sh
#@job_type=bluegene
##@arguments= --exe /gpfs/home/HCP010/lxm03/lxm77-lxm03/parafem/parafem/bin/xx7 --np 4 --args "xx7-tiny-f"
##@class=prod
#@class=interactive
#@environment=COPY_ALL
#@output=stdout.$(jobid).txt
#@error=stderr.$(jobid).txt
#@wall_clock_limit=00:20:00
#@notification=complete
#@queue

source /etc/profile.d/modules.sh

TVSCRIPT=gpfs/packages/ibm/toolworks/totalview.8.12.0-0/bin/tvscript
RUNJOB=/bgsys/drivers/ppcfloor/hlcs/bin/runjob

/gpfs/packages/ibm/toolworks/totalview.8.12.0-0/linux-power/bin/totalview $RUNJOB -a  --env-all --exe /gpfs/home/HCP010/lxm03/lxm77-lxm03/parafem/parafem/bin/xx7 --np 4 --args "xx7-tiny-f"

