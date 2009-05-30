for fil in `cat /etc/hostsalone`
do
 echo $fil doing
 ssh $fil "cd .. ; chmod a+rx jon"
 echo $fil done
done



