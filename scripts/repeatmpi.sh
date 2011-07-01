#!/bin/bash

EXPECTED_ARGS=1

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` {0 = new job 1 = restart job}"
  exit
fi


a=$1
while [ ! -e dumps/dump0040 ]
  do


    if [ $a -gt 0 ]
    then
        ####################
        # copy over last rdump (assumes reasonable size)
	# lonestar default output of date is different, so:
	#oldrdump=`ls -lart dumps/rdump* | tail -1 | awk '{print $9}'`
        # ki-jmck format:
	oldrdump=`ls -lart dumps/rdump* | tail -1 | awk '{print $8}'`
	oldrdumpsize=`ls -lart dumps/rdump* | tail -1 | awk '{print $5}'`
	
	#oldrdump2=`ls -lart dumps/rdump* | tail -2 | head -1 | awk '{print $9}'`
	oldrdump2=`ls -lart dumps/rdump* | tail -2 | head -1 | awk '{print $8}'`
	oldrdumpsize2=`ls -lart dumps/rdump* | tail -2 | head -1| awk '{print $5}'`
	
	percdiff=$(( (100*($oldrdumpsize-$oldrdumpsize2))/$oldrdumpsize2 ))

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
	    cp $oldrdump dumps/rdump-0.orig
	    cp dumps/rdump-0.orig dumps/rdump-0
	else
	    echo "Using 2nd to last rdump: $oldrdump2"
	    cp $oldrdump2 dumps/rdump-0.orig
	    cp dumps/rdump-0.orig dumps/rdump-0
	fi

	echo "Starting restart job: $a"
	mpirun -np 4 ./grmhd 1 2 2 1 1 0

    else
	echo "Starting new job: $a"
	mpirun -np 4 ./grmhd 1 2 2 1
    fi


    a=$(($a+1))

done







