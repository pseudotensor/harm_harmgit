
# from kraken (jmckinne) to ranch:
alias ls='ls'
cd /lustre/scratch/jmckinne/
mydirlist=`ls -aprt | grep / | sed "s/\///" | egrep 'thickdisk*'`
mydirlist2=`ls -aprt | grep / | sed "s/\///" | egrep 'runlocaldipole3dfiducial*'`
finaldirlist="$mydirlist $mydirlist2"
for thedir in $finaldirlist
  do
  


badend=1
counterbad=0
while [ $badend -eq 1 ]
  do


  rm -rf stderr_ranch.$thedir

  ~/bin/bbcp -z -s 20 -b +500 -a -k -f -r -P 5 -V -l stderr_ranch.$thedir -T 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' -S 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' $thedir tg802609@ranch.tacc.utexas.edu:

#  ~/bin/bbcp -z -a -k -f -r -P 5 -V -l stderr_ranch.$thedir -T 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' -S 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' $thedir tg802609@ranch.tacc.utexas.edu:


  badend=`grep "Connection closed" stderr_ranch.$thedir | wc -l | awk '{print $1}'`
  echo "counterbad=$counterbad and badend=$badend"
  counterbad=$(($counterbad+1))
  sleep 2
done




done