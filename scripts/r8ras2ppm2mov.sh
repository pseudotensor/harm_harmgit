for fil in `ls *.r8` ; do echo $fil ; r8torasjon 0 ~/bin/genimages.pal $fil ; done

for fil in `ls *.ras` ; do echo $fil ; rm -rf $fil.png ; convert -quality 100 -crop 256x256+0+0   $fil $fil.png ; done

rm -rf imagesmall*.png
list=`ls *.png`

ii=0
for fil in $list
do
    echo $fil
    ii=$(($ii+1))
    ln -s $fil imagesmall$ii.png
done

rm -rf imagesmall2.mp4
#avconv -r 10 -b 65536k -qscale 0.1   -i imagesmall%d.png imagesmall2.mp4
#avconv  -vcodec qtrle   -i imagesmall%d.png imagesmall2.mov

avconv -f image2 -r 30 -i imagesmall%d.png -vcodec qtrle -pix_fmt rgb24 -t 15 god.mov

#-qscale 1
#exit

