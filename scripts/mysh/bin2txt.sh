#!/bin/bash
#
echo $1
echo $2
echo $3
for fil in `/bin/ls dump????.dat.bin`
do 
	prefa=`echo $fil | sed "s/\./ /"`
	pref=`echo $prefa | awk '{print $1}'`
	echo $pref
	bin2txt 1 2 0 0 2 $1 $2 1 1 $pref.dat.bin $pref.dat.txt d 34
	cat $pref.dat.bin.head $pref.dat.txt > $pref.dat
        rm $pref.dat.txt
	echo $pref done
done
