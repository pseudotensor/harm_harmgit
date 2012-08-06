#!/bin/bash

# scp jon@ki-rh42.slac.stanford.edu:/data/jon/latestcode/harmgit/scripts/krakenrestartsustain_thickdisk.sh .

EXPECTED_ARGS=5


if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` {whichsystem} {basejob} {basesim} {old ext} {new ext}"
  echo "E.g. `basename $0` 1 td7 thickdisk7 nextnexte nextnextf"
  echo "E.g. `basename $0` 2 t3 thickdisk3 r s"
  exit
fi

whichsystem=$1

# 1 : lonestar
# 2 : kraken
# 3 : NAS PFE

if [ $whichsystem -eq 1 ]
then
    WORKDIR=$SCRATCH
    batchpart=tacclonestar4
fi

if [ $whichsystem -eq 2 ]
then
    WORKDIR=/lustre/scratch/$USER/
    batchpart=kraken
fi

if [ $whichsystem -eq 3 ]
then
    WORKDIR=/nobackup/$USER/
    batchpart=pfe
fi


basejob=$2
basesim=$3
####basedir="sasha"
basedir=$3

oldpwd=`pwd`
cd ~/
if [ ! -d "$basedir" ]; then
    echo "Directory does not exist, reassign basedir"
    exit
fi
cd $oldpwd

oldname=$3$4
newname=$3$5
oldjob=$2$4
newjob=$2$5

echo Using old $oldname to setup $newname



####################
# copy over stderr/stdout file to run directory
mkdir -p $WORKDIR/$oldname/stderrstdout/
mv ~/$basedir/$oldname.o* ~/$basedir/$oldname.err* $WORKDIR/$oldname/stderrstdout/


####################
# copy over last rdump (assumes reasonable size)
if [ $whichsystem -eq 1 ] ||
    [ $whichsystem -eq 3 ]
then
    oldrdump=`ls -lart $WORKDIR/$oldname/dumps/rdump-?.bin | tail -1 | awk '{print $9}'`
    oldrdumpsize=`ls -lart $WORKDIR/$oldname/dumps/rdump-?.bin | tail -1 | awk '{print $5}'`
    oldrdump2=`ls -lart $WORKDIR/$oldname/dumps/rdump-?.bin | tail -2 | head -1 | awk '{print $9}'`
    oldrdumpsize2=`ls -lart $WORKDIR/$oldname/dumps/rdump-?.bin | tail -2 | head -1| awk '{print $5}'`

    oldrdumpupperpole=`ls -lart $WORKDIR/$oldname/dumps/rdumpupperpole-?.bin | tail -1 | awk '{print $9}'`
    oldrdumpupperpolesize=`ls -lart $WORKDIR/$oldname/dumps/rdumpupperpole-?.bin | tail -1 | awk '{print $5}'`
    oldrdumpupperpole2=`ls -lart $WORKDIR/$oldname/dumps/rdumpupperpole-?.bin | tail -2 | head -1 | awk '{print $9}'`
    oldrdumpupperpolesize2=`ls -lart $WORKDIR/$oldname/dumps/rdumpupperpole-?.bin | tail -2 | head -1| awk '{print $5}'`
else
    oldrdump=`ls -lart $WORKDIR/$oldname/dumps/rdump-?.bin | tail -1 | awk '{print $8}'`
    oldrdumpsize=`ls -lart $WORKDIR/$oldname/dumps/rdump-?.bin | tail -1 | awk '{print $5}'`
    oldrdump2=`ls -lart $WORKDIR/$oldname/dumps/rdump-?.bin | tail -2 | head -1 | awk '{print $8}'`
    oldrdumpsize2=`ls -lart $WORKDIR/$oldname/dumps/rdump-?.bin | tail -2 | head -1| awk '{print $5}'`

    oldrdumpupperpole=`ls -lart $WORKDIR/$oldname/dumps/rdumpupperpole-?.bin | tail -1 | awk '{print $8}'`
    oldrdumpupperpolesize=`ls -lart $WORKDIR/$oldname/dumps/rdumpupperpole-?.bin | tail -1 | awk '{print $5}'`
    oldrdumpupperpole2=`ls -lart $WORKDIR/$oldname/dumps/rdumpupperpole-?.bin | tail -2 | head -1 | awk '{print $8}'`
    oldrdumpupperpolesize2=`ls -lart $WORKDIR/$oldname/dumps/rdumpupperpole-?.bin | tail -2 | head -1| awk '{print $5}'`
fi


percdiff=$(( (100*($oldrdumpsize-$oldrdumpsize2))/$oldrdumpsize2 ))


mkdir -p $WORKDIR/$newname/dumps/
cd $WORKDIR/$newname/dumps/

uselast=1
if [ $percdiff -gt 1 ]
then
    uselast=0
fi
if [ $percdiff -lt -1 ]
then
    uselast=0
fi

if [ $uselast -eq 1 ]
then
    echo "Using last rdump: $oldrdump"
    cp $oldrdump rdump-0.orig.bin
    cp $oldrdumpupperpole rdumpupperpole-0.bin
    cp rdump-0.orig.bin rdump-0.bin
else
    echo "Using 2nd to last rdump: $oldrdump2"
    cp $oldrdump2 rdump-0.orig.bin
    cp $oldrdumpupperpole2 rdumpupperpole-0.bin
    cp rdump-0.orig.bin rdump-0.bin
fi

# need to copy over coordparms every time for sashatilt models
cp $WORKDIR/$oldname/coordparms.dat $WORKDIR/$newname/


####################
# Setup new batch for new number using old restart file
cd ~/$basedir/
alias mv='mv'
alias cp='cp'
mv batch.qsub.$batchpart.$basesim batch.qsub.$batchpart.$basesim.orig
sed 's/'$oldname'/'$newname'/g' batch.qsub.$batchpart.$basesim.orig > batch.qsub.$batchpart.$basesim.1
sed 's/'$oldjob'/'$newjob'/g' batch.qsub.$batchpart.$basesim.1 > batch.qsub.$batchpart.$basesim.2

if [ $whichsystem -eq 1 ]
then
    sed 's/$WORK/$WORKDIR/g' batch.qsub.$batchpart.$basesim.2 > batch.qsub.$batchpart.$basesim.3
    sed 's/ibrun .\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 > $FILENAME.so/ibrun .\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 1 0 > $FILENAME.so/g' batch.qsub.$batchpart.$basesim.3 > batch.qsub.$batchpart.$basesim.4
fi

if [ $whichsystem -eq 2 ]
then
    sed 's/$SCRATCH/$WORKDIR/g' batch.qsub.$batchpart.$basesim.2 > batch.qsub.$batchpart.$basesim.3
    sed 's/aprun -n $NTOT $RUNDIR\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 *$/aprun -n $NTOT $RUNDIR\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 1 0/g' batch.qsub.$batchpart.$basesim.3 > batch.qsub.$batchpart.$basesim.4
fi

if [ $whichsystem -eq 3 ]
then
    sed 's/$SCRATCH/$WORKDIR/g' batch.qsub.$batchpart.$basesim.2 > batch.qsub.$batchpart.$basesim.3
    sed 's/mpiexec -np $NTOT $RUNDIR\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 *$/mpiexec -np $NTOT $RUNDIR\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 1 0 > output.txt/g' batch.qsub.$batchpart.$basesim.3 > batch.qsub.$batchpart.$basesim.4
fi

mv batch.qsub.$batchpart.$basesim.4 batch.qsub.$batchpart.$basesim
rm -rf batch.qsub.$batchpart.$basesim.1 batch.qsub.$batchpart.$basesim.2  batch.qsub.$batchpart.$basesim.3  batch.qsub.$batchpart.$basesim.4


# start new job
echo "Starting new job $newname"
qsub < batch.qsub.$batchpart.$basesim



