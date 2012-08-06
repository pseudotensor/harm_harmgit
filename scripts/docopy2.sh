

# from kraken (jmckinne) to ranch:
alias ls='ls'
cd /lustre/scratch/jmckinne/
#mydirlist=`ls -aprt | grep / | sed "s/\///" | egrep 'thickdisk*'`
#mydirlist2=`ls -aprt | grep / | sed "s/\///" | egrep 'runlocaldipole3dfiducial*'`
mydirlist=`ls -aprt | grep / | sed "s/\///" | egrep 'thickdiskhr3'`
mydirlist2=`ls -aprt | grep / | sed "s/\///" | egrep 'thickdisk3'`
mydirlist3=`ls -aprt | grep / | sed "s/\///" | egrep 'thickdisk10'`
mydirlist4=`ls -aprt | grep / | sed "s/\///" | egrep 'thickdiskr3'`
mydirlist5=`ls -aprt | grep / | sed "s/\///" | egrep 'thickdiskr15'`
#

# only those at end of last copy or not yet copied at all according to ranch dirs
#mydirlist="thickdiskhr3nextf thickdiskhr3nextg thickdiskhr3nexth thickdiskhr3nexti thickdiskhr3nextj thickdiskhr3nextk"
#mydirlist2="thickdisk3nextf"
#mydirlist3="thickdisk10nextf thickdisk10nextg thickdisk10nexth thickdisk10nexti"
#mydirlist4="thickdiskr3nexti thickdiskr3nextc thickdiskr3nextb thickdiskr3nextg thickdiskr3nextj thickdiskr3nextk thickdiskr3nextl thickdiskr3nextm thickdiskr3nextn"
#mydirlist5="thickdiskr15v thickdiskr15w thickdiskr15x thickdiskr15y thickdiskr15z thickdiskr15nexta"

#finaldirlist="$mydirlist"
#finaldirlist="$mydirlist $mydirlist2"
finaldirlist="$mydirlist $mydirlist2 $mydirlist3 $mydirlist4 $mydirlist5"


streams=20
#streams=4

buffers=200

# test
#scp -rp docopy2.sh tg802609@ranch.tacc.utexas.edu: &> stderr.test

# for globus-url-copy
module load globus
echo "Did you run: myproxy-logon -t 300 -l jmckinne -s myproxy.teragrid.org"

if [ 1 -eq 1 ]
then

#
for thedir in $finaldirlist
  do

  ~/bin/bbcp -z -s $streams -b +$buffers -F -a -k -f -r -P 5 -V -l stderr.$thedir -T 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' -S 'ssh -x -a -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' $thedir tg802609@ranch.tacc.utexas.edu:

#  scp -rp $thedir tg802609@ranch.tacc.utexas.edu: &> stderr.$thedir

## gridftp
#  streams=2
#  # or use -rp then relative to starting directory
#  globus-url-copy -vb -p $streams -stripe -tcp-bs 11M -r \
#      file:///lustre/scratch/jmckinne/$thedir/ \
#      gsiftp://gridftp.ranch.tacc.xsede.org/home/01014/tg802609/$thedir/ \
#      &> stderr.$thedir

done

fi


