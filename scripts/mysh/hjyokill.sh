for fil in `/bin/cat hosts.lst`
do
 echo $fil doing
 ssh $fil "killall $1"
 echo $fil done
done



