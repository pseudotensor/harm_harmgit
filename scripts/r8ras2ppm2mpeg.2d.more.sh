#~/bin/bash

# check for existence of arguments
if ([ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ] || [ -z $5 ] || [ -z $6 ]) ; then {
        echo "usage: sh r8ras2ppm2mpeg.sh doimages dogif dofli domp4 nx ny"
	echo "doimages: 0=no 1=yes.  If 0, assume they were already created and just redoing movie part of script."
	echo "dogif: 0=no 1=yes"
	echo "dofli: 0=no 1=yes"
	echo "domp4: 0=no 1=yes"
        exit
} ; fi

doimages=$1
dogif=$2
dofli=$3
domp4=$4
nx=$5
ny=$6

echo "nx=$nx ny=$ny"

# force 2 lines of data so can see image of results even if 1D
if [ $ny -eq 1 ]
then
    force2=1
    ny=10
else
    force2=0
fi

echo "new nx=$nx ny=$ny"

# start in run/images/

if [ $doimages -eq 1 ]
then

    if [ $dogif -eq 1 ]||
	[ $domp4 -eq 1 ]
    then

	rm -rf *.new.r8 *.ras
	
        # r8 2 ras
# BELOW OPTIMIZED FOR WHAT RECENTLY ONLY NEED
#	for fil in `ls *.r8`
	for fil in `ls im8p*.r8 im9p*.r8 im10p*.r8 im8c*.r8 im9c*.r8 im10c*.r8`
	do
	    echo $fil

	    if [ $force2 -eq 1 ]
	    then
		echo "force2=$force2"
		rm -rf $fil.new.r8
		#head -4 $fil >> $fil.new.r8
		# STACK 1D images to get at least $ny in size
		for ii in `seq 0 $ny`
		do
		    tail -n +5 $fil >> $fil.new.r8
		done
		r8torasjon 0 ~/bin/john.pal $fil.new.r8 $nx $ny
	    else
		cp $fil $fil.new.r8
		r8torasjon 0 ~/bin/john.pal $fil.new.r8
	    fi

	    mv $fil.new.ras $fil.ras
	   
	done
    fi

    if [ $domp4 -eq 1 ]
    then
        # for avconv:
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

if [ $nx -lt 256 ] ||
    [ $ny -lt 256 ]
then
    # avconv screws up if too small in resolution for unknown reason.  Generates boundary artifacts, for example.
    options2="-s 256x256"
#    optionsgif="-loop 0 -resize 256x256" # interpolates non-constantly
    optionsgif="-loop 0 -sample 256x256" # just fast constant interpolation so just bigger version of same data.  Still keeps gif small and quick to make.
else
    options2=""
    optionsgif="-loop 0"
fi



for (( i=0 ; i < 12 ; i++ ))
do
for (( j=0 ; j < 2 ; j++ ))
do
for (( k=0 ; k < 1 ; k++ ))
do
    # make a list file
    echo "doing $i $j $k"


    if [ $dogif -eq 1 ]
    then
	rm -rf  im${i}p${j}s${k}l.gif
	convert $optionsgif im${i}p${j}s${k}l*.ras im${i}p${j}s${k}l.gif
	rm -rf  im${i}c${j}s${k}l.gif
	convert $optionsgif im${i}c${j}s${k}l*.ras im${i}c${j}s${k}l.gif
    fi



    if [ $domp4 -eq 1 ]
    then
        # avconv
	rm -rf im${i}p${j}s${k}l.mp4
	avconv -i im${i}p${j}s${k}l%04d.r8.ras.png $options1 $options2 im${i}p${j}s${k}l.mp4
	rm -rf im${i}c${j}s${k}l.mp4
	avconv -i im${i}c${j}s${k}l%04d.r8.ras.png $options1 $options2 im${i}c${j}s${k}l.mp4

        # view with smplayer:
        # smplayer <moviename>
        # mplayer seems to have problems on my latest ubuntu system.
        # also, avconv very lossy and even introduces artifacts at low resolutions (e.g. 20x20 (28x28 with FULLOUTPUT) with RADBEAMFLAT)

    fi



    if [ $dofli -eq 1 ]
    then
        #ppm2fli
	# $fil.new.r8
	ls im${i}p${j}s${k}l????.r8.new.r8 > tmp.lis
	rm -rf im${i}p${j}s${k}l.fli
	ppm2fli -p/home/jon/bin/john.pal -N -g $nx'x'$ny  -s 0 tmp.lis im${i}p${j}s${k}l.fli
	rm -rf tmp.lis
	#
	ls im${i}c${j}s${k}l????.r8.new.r8 > tmp.lis
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


#
if [ $dogif -eq 1 ]
then
    #try1
    #convert im8p1s0l.gif im8c1s0l.gif im9p1s0l.gif im9c1s0l.gif -append immerge1.gif

    #try2
    # http://www.multipole.org/discourse-server/viewtopic.php?f=1&t=22117
    # http://www.imagemagick.org/Usage/anim_basics/#montage
    # http://www.imagemagick.org/Usage/montage/
    # http://www.imagemagick.org/Usage/resize/
    #
    #rm -rf  immerge1.gif immerge1_sidebyside1.gif immerge1_sidebyside2.gif
    #convert im8p1s0l.gif im8c1s0l.gif +append immerge1_sidebyside1.gif
    #convert im9p1s0l.gif im9c1s0l.gif +append immerge1_sidebyside2.gif
    #convert immerge1_sidebyside1.gif immerge1_sidebyside2.gif -append immerge1.gif

    
    #try3
    # http://www.imagemagick.org/Usage/anim_mods/
    FINALTEMP=immerge1.temp.gif
    FINAL=immerge1.gif
    #
    LEFT=im8p1s0l.gif
    RIGHT=im8c1s0l.gif
    TOP=immerge1_sidebyside1.gif
    convert $LEFT'[0]' -coalesce \( $RIGHT'[0]' -coalesce \) \
          +append -channel A -evaluate set 0 +channel \
          $LEFT -coalesce -delete 0 \
          null: \( $RIGHT -coalesce \) \
          -gravity East  -layers Composite  $TOP
    #
    LEFT=im9p1s0l.gif
    RIGHT=im9c1s0l.gif
    BOTTOM=immerge1_sidebyside2.gif
    convert $LEFT'[0]' -coalesce \( $RIGHT'[0]' -coalesce \) \
          +append -channel A -evaluate set 0 +channel \
          $LEFT -coalesce -delete 0 \
          null: \( $RIGHT -coalesce \) \
          -gravity East  -layers Composite  $BOTTOM
    #
    #
    convert $TOP'[0]' -coalesce \( $BOTTOM'[0]' -coalesce \) \
          -append -channel A -evaluate set 0 +channel \
          $TOP -coalesce -delete 0 \
          null: \( $BOTTOM -coalesce \) \
          -gravity South  -layers Composite  $FINALTEMP

    #
    LEFT=im10p1s0l.gif
    RIGHT=im10c1s0l.gif
    BOTTOM2=immerge1_sidebyside3.gif
    convert $LEFT'[0]' -coalesce \( $RIGHT'[0]' -coalesce \) \
          +append -channel A -evaluate set 0 +channel \
          $LEFT -coalesce -delete 0 \
          null: \( $RIGHT -coalesce \) \
          -gravity East  -layers Composite  $BOTTOM2
    #
	###############
    convert $FINALTEMP'[0]' -coalesce \( $BOTTOM2'[0]' -coalesce \) \
          -append -channel A -evaluate set 0 +channel \
          $FINALTEMP -coalesce -delete 0 \
          null: \( $BOTTOM2 -coalesce \) \
          -gravity South  -layers Composite  $FINAL

fi

