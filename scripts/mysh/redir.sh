cd $1
newdir=`echo $1 | sed 's/-newfl//'`
for fil in `/bin/cat pg.lst | awk '{print $1}'`
do
 echo $fil doing
 ssh $fil "mv $1 $newdir"
 echo $fil done
done
cd .. ; mv $1 $newdir


