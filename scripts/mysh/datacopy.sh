for fil in `/bin/cat pg.lst`
do
 echo $fil doing
 #ssh $fil "mkdir -p `pwd`"
 rsh $fil "mkdir -p `pwd`"
 #scp vel*.in* $fil:`pwd`
 rcp vel*.in* $fil:`pwd`
 echo $fil done
done
# (done!)


