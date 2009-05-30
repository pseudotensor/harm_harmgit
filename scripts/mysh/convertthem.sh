echo $1
echo $2
echo $3
for fil in `/bin/ls image????.eps`
do
        prefa=`echo $fil | sed "s/\./ /"`
        pref=`echo $prefa | awk '{print $1}'`
        echo $pref
        convert $fil $pref.ppm
        echo $pref done
done
