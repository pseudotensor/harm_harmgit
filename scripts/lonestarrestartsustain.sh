#!/bin/bash

# scp jon@ki-rh42.slac.stanford.edu:/data/jon/latestcode/harmgit/scripts/lonestarrestartsustain.sh .

EXPECTED_ARGS=2


if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` {new job number} {batch run dir and stderr/stdout location"
  exit
fi

# orig location
#/work/01014/tg802609/zb1/

mynum=$1
myoldnum=$(($mynum-1))
batchrundir=$2

echo Using old $myoldnum to setup $mynum

####################
# copy over stderr/stdout file to run directory
mkdir -p $SCRATCH/zb$myoldnum/stderrstdout/
mv $batchrundir/zb$myoldnum.o* $SCRATCH/zb$myoldnum/stderrstdout/

####################
# copy over last rdump (assumes reasonable size)
oldrdump=`ls -lart $SCRATCH/zb$myoldnum/dumps/rdump* | tail -1 | awk '{print $9}'`
oldrdumpsize=`ls -lart $SCRATCH/zb$myoldnum/dumps/rdump* | tail -1 | awk '{print $5}'`

oldrdump2=`ls -lart $SCRATCH/zb$myoldnum/dumps/rdump* | tail -2 | head -1 | awk '{print $9}'`
oldrdumpsize2=`ls -lart $SCRATCH/zb$myoldnum/dumps/rdump* | tail -2 | head -1| awk '{print $5}'`

percdiff=$(( (100*($oldrdumpsize-$oldrdumpsize2))/$oldrdumpsize2 ))


mkdir -p $SCRATCH/zb$mynum/dumps/
cd $SCRATCH/zb$mynum/dumps/

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
cd $batchrundir
alias mv='mv'
alias cp='cp'
mv batch.qsub.tacclonestar4.zakamskabig batch.qsub.tacclonestar4.zakamskabig.orig
sed 's/zb'$myoldnum'/zb'$mynum'/g' batch.qsub.tacclonestar4.zakamskabig.orig > batch.qsub.tacclonestar4.zakamskabig.1
sed 's/$WORK/$SCRATCH/g' batch.qsub.tacclonestar4.zakamskabig.1 > batch.qsub.tacclonestar4.zakamskabig.2
sed 's/ibrun .\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 > $FILENAME.so/ibrun .\/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 1 0 > $FILENAME.so/g' batch.qsub.tacclonestar4.zakamskabig.2> batch.qsub.tacclonestar4.zakamskabig.3

mv batch.qsub.tacclonestar4.zakamskabig.3 batch.qsub.tacclonestar4.zakamskabig
rm -rf batch.qsub.tacclonestar4.zakamskabig.1 batch.qsub.tacclonestar4.zakamskabig.2  batch.qsub.tacclonestar4.zakamskabig.3


# start new job
echo "Starting new job zb$mynum"
qsub < batch.qsub.tacclonestar4.zakamskabig


