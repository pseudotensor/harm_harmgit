#!/bin/bash
#
nx=$1
ny=$2
# check for existence of arguments
if ([ -z $1 ] || [ -z $2 ]) ; then {
	echo "usage: sh mkfli.sh nx ny"
	exit
} ; fi

#
# make the fli
# add -fgunzip for gzipped items
for (( i=0 ; i < 10 ; i++ ))
do
for (( j=0 ; j < 2 ; j++ ))
do
for (( k=0 ; k < 2 ; k++ ))
do
# make a list file
echo "doing $i $j $k"
#ppm2fli -p/home/jon/research/current/bin/i/john.pal -N -g $nx'x'$ny -z $nx'x'$ny -s 0 tmp.lis im${i}p.fli
#((ls images/im${i}p*.r8 > tmp.lis.$i ; ppm2fli -p/home/jon/research/current/bin/i/john.pal -N -g $nx'x'$ny  -s 0 tmp.lis.$i im${i}p.fli ; rm tmp.lis.$i) > output.txt.$i 2>&1 &)
ls images/im${i}p${j}s${k}l????.r8 > tmp.lis ; ppm2fli -p/home/jon/research/current/bin/i/john.pal -N -g $nx'x'$ny  -s 0 tmp.lis im${i}p${j}s${k}l.fli ; rm tmp.lis
ls images/im${i}c${j}s${k}l????.r8 > tmp.lis ; ppm2fli -p/home/jon/research/current/bin/i/john.pal -N -g $nx'x'$ny  -s 0 tmp.lis im${i}c${j}s${k}l.fli ; rm tmp.lis
echo "done $i $j $k"
done
done
done
#ppm2fli -g$1x$2 -N tmp.lis rho.fli
#

#
