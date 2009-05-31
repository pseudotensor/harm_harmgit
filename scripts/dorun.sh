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
#pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=10000 --short-smpl-period=2000 -e UNHALTED_CORE_CYCLES ./$RUNPROG 1 1 1
#pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=10000 --short-smpl-period=2000 -e LAST_LEVEL_CACHE_MISSES ./$RUNPROG 1 1 1
#pfmon --resolv --smpl-per-function --smpl-outfile smploutput --show-time --long-smpl-period=10000 --short-smpl-period=2000 -e DTLB_MISSES:any ./$RUNPROG 1 1 1
#pfmon -e FP_COMP_OPS_EXE ./$RUNPROG 1 1 1
#source /opt/intel/itt/tcheck/bin/32e/tcheckvars.sh
#tcheck_cl ./$RUNPROG 1 1 1
time ./$RUNPROG 1 1 1
cd ..
cd ..

