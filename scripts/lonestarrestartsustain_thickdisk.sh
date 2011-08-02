#!/bin/bash

# scp jon@ki-rh42.slac.stanford.edu:/data/jon/latestcode/harmgit/scripts/lonestarrestartsustain_thickdisk.sh .

EXPECTED_ARGS=4


if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` {basejob} {basesim} {old ext} {new ext}"
  echo "E.g. `basename` td7 thickdisk7 nextnexte nextnextf"
  exit
fi


basejob=$1
basesim=$2
oldname=$2$3
newname=$2$4
oldjob=$1$3
newjob=$1$4

echo Using old $oldname to setup $newname

####################
# copy over stderr/stdout file to run directory
mkdir -p $SCRATCH/$oldname/stderrstdout/
mv ~/$basesim/$oldname.o* $SCRATCH/$oldname/stderrstdout/

####################
# copy over last rdump (assumes reasonable size)
oldrdump=`ls -lart $SCRATCH/$oldname/dumps/rdump* | tail -1 | awk '{print $9}'`
oldrdumpsize=`ls -lart $SCRATCH/$oldname/dumps/rdump* | tail -1 | awk '{print $5}'`

oldrdump2=`ls -lart $SCRATCH/$oldname/dumps/rdump* | tail -2 | head -1 | awk '{print $9}'`
oldrdumpsize2=`ls -lart $SCRATCH/$oldname/dumps/rdump* | tail -2 | head -1| awk '{print $5}'`

percdiff=$(( (100*($oldrdumpsize-$oldrdumpsize2))/$oldrdumpsize2 ))


mkdir -p $SCRATCH/$newname/dumps/
cd $SCRATCH/$newname/dumps/

uselast=1
if [ $percdiff -gt 10 ]
then
    uselast=0
fi
if [ $percdiff -lt -10 ]
then
    uselast=0
fi

if [ $uselast -eq 1 ]
then
    echo "Using last rdump: $oldrdump"
    cp $oldrdump rdump-0.orig.bin
    cp rdump-0.orig.bin rdump-0.bin
else
    echo "Using 2nd to last rdump: $oldrdump2"
    cp $oldrdump2 rdump-0.orig.bin
    cp rdump-0.orig.bin rdump-0.bin
fi


####################
# Setup new batch for new number using old restart file
cd ~/$basesim/
alias mv='mv'
alias cp='cp'
mv batch.qsub.tacclonestar4.$basesim batch.qsub.tacclonestar4.$basesim.orig
sed 's/'$oldname'/'$newname'/g' batch.qsub.tacclonestar4.$basesim.orig > batch.qsub.tacclonestar4.$basesim.1
sed 's/'$oldjob'/'$newjob'/g' batch.qsub.tacclonestar4.$basesim.1 > batch.qsub.tacclonestar4.$basesim.2
sed 's/$WORK/$SCRATCH/g' batch.qsub.tacclonestar4.$basesim.2 > batch.qsub.tacclonestar4.$basesim.3
sed 's/ibrun .\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 > $FILENAME.so/ibrun .\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 1 0 > $FILENAME.so/g' batch.qsub.tacclonestar4.$basesim.3 > batch.qsub.tacclonestar4.$basesim.4

mv batch.qsub.tacclonestar4.$basesim.4 batch.qsub.tacclonestar4.$basesim
rm -rf batch.qsub.tacclonestar4.$basesim.1 batch.qsub.tacclonestar4.$basesim.2  batch.qsub.tacclonestar4.$basesim.3  batch.qsub.tacclonestar4.$basesim.4


# start new job
echo "Starting new job $newname"
qsub < batch.qsub.tacclonestar4.$basesim



