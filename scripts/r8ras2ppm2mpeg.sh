#~/bin/bash

# check for existence of arguments
if ([ -z $1 ] || [ -z $2 ] || [ -z $3 ]) ; then {
        echo "usage: sh r8ras2ppm2mpeg.sh which doimages nx ny"
	echo "which: 0=all 1=fli only 2=mp4 only"
	echo "doimages: 0=no 1=yes.  If 0, assume they were already created and just redoing movie part of script."
        exit
} ; fi

which=$1
doimages=$2
nx=$3
ny=$4

# start in run/images/

if [ $doimages -eq 1 ]
then
if [ $which -eq 0 ] ||
    [ $which -eq 2 ]
then
    # for avconv:
    
    # r8 2 ras
    rm -rf *.ras
    for fil in `ls *.r8` ; do echo $fil ; r8torasjon 0 ~/bin/john.pal $fil ; done
    
    # ras to lossless png
    rm -rf *.png
    for fil in `ls *.ras` ; do echo $fil ; convert $fil $fil.png ; done
fi
fi


# png 2 mpg
# use avconv instead ubuntu says, but still works.

# Default ubuntu's avconv or ffmpeg is very limited due to license restrictions, but see:
# http://www.medibuntu.org/

#options1="-vcodec libx264 -threads 8 -acodec aac -ab 128k -r 30 -b 65536k"
#options1="-vcodec mpeg4 -threads 8 -acodec aac -ab 128k -r 30 -b 65536k"
options1="-vcodec mpeg4 -threads 8 -r 30 -b 65536k"
#options1="-vcodec ljpeg -threads 8 -r 30 -b 65536k"

if [ $nx -lt 50 ] ||
    [ $ny -lt 50 ]
then
    # avconv screws up if too small in resolution for unknown reason.  Generates boundary artifacts, for example.
    options2="-s 50x50"
else
    options2=""
fi


for (( i=0 ; i < 12 ; i++ ))
do
for (( j=0 ; j < 2 ; j++ ))
do
for (( k=0 ; k < 1 ; k++ ))
do
    # make a list file
    echo "doing $i $j $k"

    if [ $which -eq 0 ] ||
	[ $which -eq 2 ]
    then
        # avconv
	rm -rf im${i}p${j}s${k}l.mp4
	avconv -i im${i}p${j}s${k}l%04d.ras.png $options1 $options2 im${i}p${j}s${k}l.mp4
	rm -rf im${i}c${j}s${k}l.mp4
	avconv -i im${i}c${j}s${k}l%04d.ras.png $options1 $options2 im${i}c${j}s${k}l.mp4

        # view with smplayer:
        # smplayer <moviename>
        # mplayer seems to have problems on my latest ubuntu system.
        # also, avconv very lossy and even introduces artifacts at low resolutions (e.g. 20x20 (28x28 with FULLOUTPUT) with RADBEAMFLAT)

    fi



    if [ $which -eq 0 ] ||
	[ $which -eq 1 ]
    then
        #ppm2fli
	ls im${i}p${j}s${k}l????.r8 > tmp.lis
	rm -rf im${i}p${j}s${k}l.fli
	ppm2fli -p/home/jon/bin/john.pal -N -g $nx'x'$ny  -s 0 tmp.lis im${i}p${j}s${k}l.fli
	rm -rf tmp.lis
	#
	ls im${i}c${j}s${k}l????.r8 > tmp.lis
	rm -rf im${i}c${j}s${k}l.fli
	ppm2fli -p/home/jon/bin/john.pal -N -g $nx'x'$ny  -s 0 tmp.lis im${i}c${j}s${k}l.fli
	rm -rf tmp.lis

	# use xanim:
	# xanim -WrT2 <moviename>
	# the fli file will be lossless:
	# http://www.aos.wisc.edu/fli.html
	# mplayer should run .fli as well, but doesn't work on my new ubuntu system.

    fi

    echo "done $i $j $k"

done
done
done


