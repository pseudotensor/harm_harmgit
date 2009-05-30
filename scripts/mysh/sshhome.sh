# Goto: http://www-unix.mcs.anl.gov/mpi/mpich/
cd /home/jon
tar cvzpf sshhome.tgz .ssh
for fil in `/bin/cat /etc/hostsalone`
do
 echo $fil doing
 scp /home/jon/sshhome.tgz $fil:/home/jon/
 ssh $fil "cd /home/jon ; tar xvzpf sshhome.tgz"
 echo $fil done
done
# (done!)


