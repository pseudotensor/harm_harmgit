#!/bin/bash

# call using, e.g.:
# (nohup sh ./loopjoninterp_simple.sh 0 1664 24 &> interps.txt &)

mkdir outs/
startdump=$1
enddump=$2
skipdump=$3

dumpi=$startdump
while [ $dumpi -le $enddump ]
do

    nextdumpi=$(($dumpi+$skipdump))

    echo "Doing file $dumpi through $nextdumpi"
    inneri=0
    while [ $inneri -le $skipdump ]
    do
        tempdumpi=`echo $(( 10#$dumpi ))`
        textdumpi=`printf "%04d" "$dumpi"`

        # assume howtouse_joninterp.sh can be called with number as first argument and whichres as second argument
        #(nohup sh howtouse_joninterp.biggerbox.sh $textdumpi 1 &> outs/out.$dumpi.txt &)
        sh howtouse_joninterp.biggerbox.sh $textdumpi 1 &> outs/out.$dumpi.txt &

        dumpi=$(($dumpi+1))
        inneri=$(($inneri+1))
    done

    # wait for all $skipdump dumps to be done to limit use of limited cores
    wait

done

echo "All Done."
