#~/bin/bash

# start in run/images/


#then:

# r8 2 ras
for fil in `ls *.r8` ; do echo $fil ; r8torasjon 0 ~/bin/john.pal $fil ; done

# ras to jpg
for fil in `ls *.ras` ; do echo $fil ; convert $fil $fil.jpg ; done

# jpg 2 mpg
# use avconv instead ubuntu says, but still works.
list='0 1 2 3 4 5 6 7 8 9 10 11'
for it in $list
do
    echo "Doing: $it"
    ffmpeg -y -fflags +genpts -r 30 -i im${it}p1s0l%04d.ras.jpg -sameq -qmax 5 im${it}p1s0l.mov
done

# view with smplayer:
# smplayer <moviename>

# mplayer seems to have problems on my system.