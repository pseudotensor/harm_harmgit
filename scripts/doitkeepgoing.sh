#!/bin/bash

# run as (e.g.) : (nohup sh doitfull3d7tilt1.5708.sh &> doitfull3d7tilt1.5708.out &)
system=3
user=jmckinn2
jobprefix=tt1.5
jobdirprefix=thickdiskfull3d7tilt1.5708
cd ~/thickdiskfull3d7tilt1.5708/


#########################
listlet1=`eval echo {a..z}`
listlet2=`echo $listlet1 | sed 's/\([a-z]\)/next\1/g'`
listlet3=`echo $listlet2 | sed 's/next/nextnext/g'`
listlet="$listlet1 $listlet2 $listlet3"
arraylet=($listlet)

#################
numjobs=`echo $listlet | wc | awk '{print $2}'`

seqstart=1 # if set this to restart the job, need to ensure batch file is chosen with same letter.

seqend=$(($numjobs-1))


for job in `seq $seqstart $seqend`
do

    oldlet=`echo ${arraylet[$job]}`
    newlet=`echo ${arraylet[$(($job+1))]}`
    sleep 0
    echo "sh krakenrestartsustain_thickdisk.sh $system $jobprefix $jobdirprefix $oldlet $newlet"
    sh krakenrestartsustain_thickdisk.sh $system $jobprefix $jobdirprefix $oldlet $newlet

    # wait for job to get submitted after copying files
    sleep 60
    didqueue=0
    inqueue=1
    while [ $inqueue -eq 1 ]
    do
	myjobstatus=`qstat | grep $user | grep $jobprefix | awk '{print $8}' | grep "Q" | wc | awk '{print $2}'`

	if [ $myjobstatus -eq 1 ]
	then
	    echo "job=$job in queue"
	    inqueue=1
	    didqueue=1
	else
	    inqueue=0
	fi
	sleep 60
    done
    
    # wait for job to show up as running
    sleep 20
    didrun=0
    inrun=1
    while [ $inrun -eq 1 ]
    do
	myjobstatus=`qstat | grep $user | grep $jobprefix | awk '{print $8}' | grep "R" | wc | awk '{print $2}'`

	if [ $myjobstatus -eq 1 ]
	then
	    echo "job=$job running"
	    inrun=1
	    didrun=1
	else
	    inrun=0
	fi
	sleep 60
    done
    
    if [ $didrun -eq 0 ]
    then
	echo "Problem with job=$job.  Did not ever run apparently."
	exit
    fi
    
done

echo "Done with all jobs"

cd ~/
