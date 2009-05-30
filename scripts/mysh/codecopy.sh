# Goto: http://www-unix.mcs.anl.gov/mpi/mpich/
for fil in `/bin/cat /etc/hostsalone`
do
 echo $fil doing
 scp .bash_profile $fil:
 echo $fil done
done
# (done!)


