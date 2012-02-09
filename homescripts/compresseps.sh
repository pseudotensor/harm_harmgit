#!/bin/bash
#Usage:
#
# compresseps.sh input.eps output.eps
#
fnamein=$1
fnameout=$1

mv $fnamein $fnamein.orig
fnamein=$1.orig


fnamenoext=${1%\.*} #input file without extension
#choose max resolution:
maxdimenx=2048
maxdimeny=2048
#choose snug bounding box
ps2epsi $fnamein ${fnamenoext}.epsi
#convert the EPS/PS file to a PPM file with dimensions given by $maxdimenx and $maxdimeny:
pstopnm -xborder=0 -yborder=0 ${fnamenoext}.epsi -xmax=$maxdimenx -ymax=$maxdimeny -portrait  -stdout -ppm >${fnamein}.ppm
#first, convert to JPG to compress
echo "1"
echo ${fnamein}.ppm ${fnamein}.jpg
sam2p ${fnamein}.ppm ${fnamein}.jpg
#then, convert to desired format (presumably, EPS):
echo "2"
echo ${fnamein}.jpg ${fnameout}
sam2p ${fnamein}.jpg ${fnameout}


