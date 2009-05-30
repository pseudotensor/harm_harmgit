for fil in `/bin/cat pg.lst | awk '{print $1}'`
do
 echo $fil doing
 rsh $fil "rm -rf `pwd`"
 echo $fil done
done



