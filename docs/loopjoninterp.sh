#!/bin/bash
#################### Start LOOP


# takes about 82 minutes for interp+bin2txt+v5dppm
# seems to take 7 minutes per file just for v5d->ppm part
# need about 4GB per core
# need 64 cores

# qsub -I -A TG-PHY120005  -q analysis -l ncpus=64,mem=256GB,walltime=2:30:00

# qsub -I -A TG-PHY120005  -q analysis -l ncpus=256,mem=1024GB,walltime=2:00:00


# http://stackoverflow.com/questions/10987246/xvfb-multiple-displays-for-parallel-processing

(nohup Xvfb -noreset -ac :2 &)

cd /lustre/medusa/jmckinne/data3/jmckinne/jmckinne/sashaa99t1.5708/




#startdump=0
#enddump=5737
startdump=$1
enddump=$2

#skipdump=100

skipdump=$3

dumpi=$startdump

while [ $dumpi -le $enddump ]
do

    tempdumpi=`echo $(( 10#$dumpi ))`

    god=`printf "%04d" "$dumpi"`

    echo "Doing file $god"
    
    sed 's/dumpnum=5736/dumpnum='$god'/g' howtouse_joninterp.sh > howtouse_joninterp.$god.sh

    # skip over, assume on Nautilus and have many cores.
    #(nohup sh howtouse_joninterp.$god.sh &> out.$god.txt &)
    sh howtouse_joninterp.$god.sh &> out.$god.txt &

    dumpi=$(($dumpi+$skipdump))
done

wait
