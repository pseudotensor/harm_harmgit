alias ls='ls'
cd /lustre/scratch/jmckinne/
mydirlist=`ls -aprt | grep / | sed "s/\///" | egrep 'thickdisk*'`
mydirlist2=`ls -aprt | grep / | sed "s/\///" | egrep 'runlocaldipole3dfiducial*'`
mydirlist="$mydirlist $mydirlist2"

dirorig=`pwd`


rm -rf stdfull_ranch* stderr_ranch*


#streams=8
streams=20
#buffers=2000
for thedir in $mydirlist
do

mysubdirlist="/dumps/ /images/ /"
# try just getting dumps and images over for now
#mysubdirlist="/dumps/ /images/"
mysubdirlistname="dumps images base"
isubdir=0
for thesubdir in $mysubdirlist
do

if [ $isubdir -eq 2 ]
then
    thesubdirname="base"
fi
if [ $isubdir -eq 0 ]
then
    thesubdirname="dumps"
fi
if [ $isubdir -eq 1 ]
then
    thesubdirname="images"
fi
isubdir=$(($isubdir+1))


badmkdir=1
while [ $badmkdir -eq 1 ]
  do
  ssh tg802609@ranch.tacc.utexas.edu "mkdir -p $thedir/$thesubdir/" &> $dirorig/sshmkdir.$thedir.$thesubdirname
  badmkdir=`grep "Connection closed" $dirorig/sshmkdir.$thedir.$thesubdirname | wc -l | awk '{print $1}'`
  echo "badmkdir=$badmkdir"
  sleep 5
done


rm -rf stderr_ranch.$thedir.$thesubdirname*
rm -rf stdfull_ranch.$thedir.$thesubdirname*
rm -rf sshmkdir.$thedir.$thesubdirname*


badend=1
badend1=1
badend2=1
counterbad=0
counterbad1=0
counterbad2=0
while [ $badend -gt 0 ]
  do

  rm -rf stderr_ranch.$thedir.$thesubdirname.$counterbad
  rm -rf stdfull_ranch.$thedir.$thesubdirname.$counterbad
  rm -rf sshmkdir.$thedir.$thesubdirname.$counterbad
  
  cd $thedir/$thesubdir/
  #filestocopy=`find . -maxdepth 1 -type f`
  #numfiles=`echo $filestocopy | wc -w | awk '{print $1}'`
  find . -maxdepth 1 -type f > $dirorig/filestocopy.txt.$thedir.$thesubdirname
  numfiles=`cat $dirorig/filestocopy.txt.$thedir.$thesubdirname | wc -w | awk '{print $1}'`
  echo "numfiles=$numfiles"
  #
  # -I $dirorig/filestocopy.txt.$thedir.$thesubdirname
  filelist=`cat $dirorig/filestocopy.txt.$thedir.$thesubdirname`
  echo $filelist
  #
  #-b +$buffers
  ~/bin/bbcp -z -s $streams  -a -k -f -P 5 -V -l $dirorig/stderr_ranch.$thedir.$thesubdirname.$counterbad -T 'ssh -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' -S 'ssh -oConnectTimeout=0 -oFallBackToRsh=no %I -l %U %H ~/bin/bbcp' -I $dirorig/filestocopy.txt.$thedir.$thesubdirname  tg802609@ranch.tacc.utexas.edu:$thedir/$thesubdir/ &> $dirorig/stdfull_ranch.$thedir.$thesubdirname.$counterbad

  #~/bin/bbcp -z -a -k -f -P 5 -V -l $dirorig/stderr_ranch.$thedir.$thesubdirname.$counterbad  -I $dirorig/filestocopy.txt.$thedir.$thesubdirname tg802609@ranch.tacc.utexas.edu:$thedir/$thesubdir/ &> $dirorig/stdfull_ranch.$thedir.$thesubdirname.$counterbad


  cd $dirorig

  badend1=`grep "Connection closed" $dirorig/stderr_ranch.$thedir.$thesubdirname.$counterbad | wc -l | awk '{print $1}'`
  goodend2=`grep "files copied at effectively" $dirorig/stdfull_ranch.$thedir.$thesubdirname.$counterbad | wc -l | awk '{print $1}'`
  numskip=`grep "copy skipped" $dirorig/stderr_ranch.$thedir.$thesubdirname.$counterbad |wc -l | awk '{print $1}'`
  #numfiles=`find | wc -l | awk '{print $1}'`
  goodend3=0
  if [ $numfiles -le $numskip ]
      then
      goodend3=1
  fi

  badend2=0
  if [ $goodend2 -eq 0 ] &&
      [ $goodend3 -eq 0 ]
      then
      badend2=1
  fi

  #badend=$(($badend1+$badend2))
  badend=$(($badend1))

  #
  echo "---------------------------------------------"
  echo "thedir=$thedir"
  echo "thesubdirname=$thesubdirname"
  echo "numfiles=$numfiles and numskip=$numskip"
  echo "goodend2=$goodend2 and goodend3=$goodend3"
  echo "counterbad=$counterbad and badend=$badend"
  counterbad=$(($counterbad+1))
  echo "counterbad1=$counterbad1 and badend1=$badend1"
  counterbad1=$(($counterbad1+1))
  echo "counterbad2=$counterbad2 and badend2=$badend2"
  counterbad2=$(($counterbad2+1))

  sleep 10
# end over bad attempts
done

# end over subdirs
done


# scp (no progress shown)
#  scp -rp $thedir tg802609@ranch.tacc.utexas.edu:  &>  stderr.$thedir

# end over thedirs
done
