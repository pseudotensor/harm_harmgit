#!/bin/bash
alias mydiff='diff -bBdp'
alias lssdir2='ls -ap| grep / | tail -n +3 | sed "s/\///"'

mydirs=`lssdir2`
mydirs="./ $mydirs"

localdir=`pwd`

for mydir in $mydirs
do
    echo $mydir
    cd $mydir
    #
    if [ "$mydir" -eq "./" ]
    then
	thatdir="../$1/$mydir/"
    else
	thatdir="$1/$mydir/"
    fi
    #
    thisdir="$localdir/$mydir/"
    #
    echo "$localdir/$mydir/$fil ../$1/$fil"
    #
    files=`ls`
    #
    numfiles=`echo $files |wc | awk '{print $2}'`
    #
    if [ $numfiles -ne 0 ]
    then
	for fil in $files; do echo $fil; mydiff $thisdir/$fil $thatdir/$fil; done
    fi
    #
    cd ..
    #
    break
    #
done
