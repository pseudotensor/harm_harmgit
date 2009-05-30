cdir=`pwd`
result=`echo $HOSTNAME | sed  's/\..*//'`

echo This host: $result
echo

for fil in `/bin/cat pg.lst | awk '{print $1}'`
do
  if [ $result != $fil ]
      then
      echo $fil to do
  fi
done
read cancontinue

#if [ "$cancontinue" == "y" ]
#then
#    echo yes
#fi
#if [ "$cancontinue" == "n" ]
#then
#    echo no
#fi

if [ "$cancontinue" == "y" ]
then
    for fil in `/bin/cat pg.lst | awk '{print $1}'`
      do
      if [ $result != $fil ]
	  then
	  echo $fil doing
	  #ssh $fil "rm -rf `pwd`"
	  echo $fil done
      fi
    done
fi



