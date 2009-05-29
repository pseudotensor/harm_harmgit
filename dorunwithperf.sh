#!/bin/sh
# e.g.
#  nohup sh dorun.sh test1103.fluxreconfull.fluxcthll.4x1 &
#
# below prefix to avoid mistakes
MAKEDIR=makedir.$1
RUNDIR=run
RUNPROG=grmhd.$MAKEDIR
#
mkdir -p $MAKEDIR
rm -rf $MAKEDIR/*
cp * $MAKEDIR
cd $MAKEDIR
make superclean
make prep
make
mv grmhd $RUNPROG
mkdir -p $RUNDIR
rm -rf $RUNDIR/*
cp *.dat $RUNDIR
cp $RUNPROG $RUNDIR
cd $RUNDIR
#time ./$RUNPROG 1 1 1


# Note that using pfmon will slow code down, just like -pg and gprof.  About 2X slow for below settings, but depends upon event looking at.

#################################
# OF INTEREST:
#
# L1I_MISSES
# BR_MISSP_EXEC
# BR_RET_MISSP_EXEC

# LAST_LEVEL_CACHE_REFERENCES,LAST_LEVEL_CACHE_MISSES

# BRANCH_INSTRUCTIONS_RETIRED,MISPREDICTED_BRANCH_RETIRED

# DTLB_MISSES,MEMORY_DISAMBIGUATION,PAGE_WALKS

# L2 cache reads, misses, evicted,all requests
# L2_LD:any,L2_LINES_IN:any,L2_LINES_OUT:any,L2_RQSTS:any

# total cycles
# UNHALTED_CORE_CYCLES

##################
# DOESN'T WORK:
# L2_LD:any : doesn't seem to work, give no result, or causes all events to give no result:
# L2_RQSTS:any : no result

###################
# WORKS: (do one at a time):
# UNHALTED_CORE_CYCLES
# DTLB_MISSES:any
# L2_LINES_IN:any
# LAST_LEVEL_CACHE_MISSES
# L1D_PEND_MISS
# BR_INST_EXEC : total branches executed
# BR_MISSP_EXEC : total mispredicted branches executed
# BR_RET_MISSP_EXEC
# X87_OPS_RETIRED : for flops [then don't do sample related things)

###############################
# Some attempted commands
#
#pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=10000 -e UNHALTED_CORE_CYCLES,L2_LD:any,L2_LINES_IN:any,L2_LINES_OUT:any,L2_RQSTS:any ./$RUNPROG 1 1 1


# GOOD ONE:
# pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=10000 -e L2_LD:any,L2_LINES_IN:any

# pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time -e L2_LD:any,L2_LINES_IN:any


#pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=10000 -eUNHALTED_CORE_CYCLES,BR_RET_MISSP_EXEC,BR_MISSP_EXEC,LAST_LEVEL_CACHE_MISSES,L1D_PEND_MISS,DTLB_MISSES:any ./$RUNPROG 1 1 1

#pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=100 -e L2_LD:any,L2_LINES_IN:any ./$RUNPROG 1 1 1


# below works and counts L2 cache misses per function call
#pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=100 -e L2_LINES_IN:any ./$RUNPROG 1 1 1


make
cp grmhd grmhd.makedir.testpfmon

rm -rf run/* ; cp grmhd.makedir.testpfmon run ; cd run
export RUNPROG=./grmhd.makedir.testpfmon
pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=10000 --short-smpl-period=2000 -e UNHALTED_CORE_CYCLES  ./$RUNPROG 1 1 1




cd ..
cd ..



# CPU and Cache info:
# dmesg | grep CPU 
# cat /proc/cpuinfo

# ki-rh42:
# Intel(R) Core(TM)2 Extreme CPU Q6800  @ 2.93GHz stepping 0b
# Core2 CPU: CPU: L1 I cache: 32K, L1 D cache: 32K
# Core2 CPU: CPU: L2 cache: 4096K


