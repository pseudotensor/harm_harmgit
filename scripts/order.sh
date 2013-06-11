iter=0
list1=`find -maxdepth 1 -name "lrho[0-9][0-9][0-9][0-9]_Rzxym1.png" | sort -g | sed 's/.\///g'`
list2=`find -maxdepth 1 -name "lrho[0-9][0-9][0-9][0-9][0-9]_Rzxym1.png" | sort -g | sed 's/.\///g'`
listfull="$list1 $list2"
for fil in $listfull
do
    myfil=`echo $fil | sed 's/lrho[0-9]\+/lrho/g'`
    myfil=$myfil.1.jon$iter.png
    echo $fil $myfil
    ln -s $fil $myfil
    iter=$(($iter+1))
done

iter=0
list1=`find -maxdepth 1 -name "lrhosmall[0-9][0-9][0-9][0-9]_Rzxym1.png" | sort -g | sed 's/.\///g'`
list2=`find -maxdepth 1 -name "lrhosmall[0-9][0-9][0-9][0-9][0-9]_Rzxym1.png" | sort -g | sed 's/.\///g'`
listfull="$list1 $list2"
for fil in $listfull
do
    myfil=`echo $fil | sed 's/lrhosmall[0-9]\+/lrhosmall/g'`
    myfil=$myfil.2.jon$iter.png
    echo $fil $myfil
    ln -s $fil $myfil
    iter=$(($iter+1))
done

modelname="runnorad1"

#NEW: 
# avconv -i lrho_Rzxym1.png.1.jon%d.png lrho.$modelname.mov
# avconv -i lrhosmall_Rzxym1.png.2.jon%d.png lrhosmall.$modelname.mov

#NEWER:
options1="-vcodec mpeg4 -threads 8 -r 30 -b 65536k"
avconv -i lrho_Rzxym1.png.1.jon%d.png $options1 lrho.$modelname.mp4
avconv -i lrhosmall_Rzxym1.png.2.jon%d.png $options1 lrhosmall.$modelname.mp4

# OLD:
# ffmpeg -y -fflags +genpts -r $fps -i lrho%04d_Rzxym1.png -sameq -qmax 5 lrho.$modelname.mov

# ffmpeg -y -fflags +genpts -r $fps -i lrho_Rzxym1.png.1.jon%d.png -sameq -qmax 5 lrho.$modelname.mov
