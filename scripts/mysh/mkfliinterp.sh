# e.g. sh mkfliinterp.sh 0 0 1 1 256 256 1.0 0 0 256 512 1.96 20 0 20 20 .05 0 0 0 0 i1
# e.g. sh mkfliinterp.sh 0 0 1 1 456 456 1.0 0 0 256 512 1.321 40 0 40 40 0.3 0 0 0 0 i3
# sh mkfliinterp.sh 0 2 1 1 512 256 1  1.0 0.0 0.5 0  512 1024 1  0 500 -500 500 0 0  1.1 10000 -3 0.3 9 0 0 0 ii
#!/bin/bash
#
datatype=$1
interptype=$2
readheader=$3
writeheader=$4
nx=$5
ny=$6
nz=$7
refinefactor=$8
filter=$9
sigma=${10}
gridtype=${11}
nxnew=${12}
nynew=${13}
nznew=${14}
xin=${15}
xout=${16}
yin=${17}
yout=${18}
zin=${19}
zout=${20}
rin=${21}
rout=${22}
R0=${23}
hslope=${24}
defcoord=${25}

whichp=${26}
whichs=${27}
whichk=${28}
newdir=${29}
# check for existence of arguments
if ([ -z $1 ] || [ -z $2 ]) ; then {
	echo "usage: sh mkfli.sh nx ny"
	exit
} ; fi

#
# make the fli
# add -fgunzip for gzipped items
#ppm2fli -p/home/jon/research/current/bin/i/john.pal -N -g $nx'x'$ny -z $nx'x'$ny -s 0 tmp.lis im${i}p.fli
#((ls images/im${i}p*.r8 > tmp.lis.$i ; ppm2fli -p/home/jon/research/current/bin/i/john.pal -N -g $nx'x'$ny  -s 0 tmp.lis.$i im${i}p.fli ; rm tmp.lis.$i) > output.txt.$i 2>&1 &)
i=$whichp
j=$whichs
k=$whichk

mkdir ${newdir}images
ls images/im${i}p${j}s${k}l????.r8 > tmp.lis
for fil in `cat tmp.lis`
  do
  echo "doing $fil"
  ~/sm/iinterp $datatype $interptype $readheader $writeheader $nx $ny $nz $refinefactor $filter $sigma $gridtype $nxnew $nynew $nznew $xmin $xmax $ymin $ymax $zmin $zmax $rin $rout $R0 $hslope $defcoord < $fil > ${newdir}${fil}
  echo "done $fil"
done
ls ${newdir}images/im${i}p${j}s${k}l????.r8 > ${newdir}tmp.lis
#
#ppm2fli -p/home/jon/research/current/bin/i/john.pal -N -g $nxnew'x'$nynew  -s 0 tmp.lis ${newdir}im${i}p${j}s${k}l.fli ; rm tmp.lis
# for no header
ppm2fli -p/home/jon/research/current/bin/i/john.pal -N -g $nxnew'x'$nynew -z $nxnew'x'$nynew  -s 0 ${newdir}tmp.lis ${newdir}im${i}p${j}s${k}l.fli ; rm tmp.lis

#

#
