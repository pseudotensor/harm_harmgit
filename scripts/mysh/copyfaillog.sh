cdir=`pwd`
result=`echo $HOSTNAME | sed  's/\..*//'`

for fil in `/bin/cat pg.lst | awk '{print $1}'`
do
  if [ $result != $fil ]
      then
      echo $fil doing
      rcp $fil:$cdir/0_*.out.* .
      rcp $fil:$cdir/0_*.out .
      rcp $fil:$cdir/probe* .
 #ssh $fil "cd /$cdir/ ; rm -rf 0_*.out.*"
      echo $fil done
  fi
done



