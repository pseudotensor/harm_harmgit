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
time ./$RUNPROG
cd ..
cd ..

