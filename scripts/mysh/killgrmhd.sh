for fil in `/bin/cat /etc/hostsallandbh | awk '{print $1}'`
do
 echo $fil doing
 ssh $fil 'for ofil in `ps -A | grep grmhd | awk '{print $1}'` ; do echo $ofil; done'
 echo $fil done
done



