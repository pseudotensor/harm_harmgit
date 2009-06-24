#!/bin/bash

# need to set:
# MPI_MAX_CLUSTER_SIZE=16

# Create pg.lst file that contains:
# name numprocsthisnode
# name2 numprocsthisnode
# etc.
# E.g.:
#
# ki-rh42.slac.stanford.edu 1
# ki-rh42.slac.stanford.edu 1
#
# then run (for 1 node, many procs):
# mympirun.sh 0 pg.lst ./cpi 100
#




#defaults from command line

# number of arguments related to mpi (non-program args) + the program name
numargs=3

if [ ${#} -lt $numargs ]
    then
    echo "Need at least $numargs arguments, you gave ${#}"
    echo "e.g. mympirun 1 pg.lst grmhd"
    exit 1
fi
# default number of procs per machine
defnp=1
# local user
defuser=`whoami`
defpath=`pwd`
# whether to copy local file to new path/file
docopy=$1
pgfile=$2
newpgfile=$pgfile.new
rm $newpgfile
# can only support different users if given full path in proc file
#
defprogold=$3
# local file (full path)
localfullpath=$defpath/`basename $defprogold`
defprog=$localfullpath
#
# defaults for script
p4np=$defnp
p4user=$defuser
p4prog=$localfullpath


# first extract the desired arguments
start=$(($numargs+1))
end=$((${#}+0))

a=""
counter=0
for i in "$@"
  do
  #echo \"$i\"
  counter=$(($counter+1))
  #echo $counter
  if [[ $counter -ge $start && $counter -le $end ]]
      then
      #a="$a \"$i\""
      a="$a $i"
      fi
  #echo $a
done
echo "program arguments:"
echo $a

# now copy and run the code
nodenum=1
numprocs=0

#totalnodes=`wc -l $pgfile | awk '{print $1}'`
#while [ $nodenum -le $totalnodes ]
dolr='$'
for node in `cat $pgfile | awk '{print $1}'`
do
  #node=cat $pgfile | sed -n '$nodenum p' | awk '{print $1}'
  echo "doing machine #$nodenum named $node"
  # determine np
  p4nptemp=`cat $pgfile | sed -n "$nodenum p" | awk '{print $2}'`
  echo $p4nptemp
  if [ $p4nptemp ]
      then
      p4np=$p4nptemp
      if [ $nodenum -eq 1 ]
	  then
	  p4np=$(($p4np-1))
      fi
  else
      if [ $nodenum -eq 1 ]
	  then
	  p4np=$(($defnp-1))
      else
	  p4np=$defnp
      fi
  fi
  # determine program name
  p4progtemp=`cat $pgfile  | grep $node | awk '{print $3}'`
  if [ $p4progtemp ]
      then
      p4prog=$p4progtemp
  else
      p4prog=$defprog
  fi
	#	p4path=${p4prog%\/*}
  # determine path
  p4path=`dirname $p4prog`
  # determine user
  p4usertemp=`cat $pgfile  | grep $node | awk '{print $4}'`
  if [ $p4usertemp ]
      then
      p4user=$p4usertemp
  else
      p4user=$defuser
  fi
  if [ $nodenum -eq 1 ]
      then
      masternode=$node
      masterpath=$p4path
      masterprog=$p4prog
      masteruser=$p4user
  fi
  #echo $node $p4np $p4prog $p4user
  #echo $node $p4path $pgfile $localfullpath $node:$p4path
  if [ $docopy -eq 1 ]
      then
      #( nohup ssh $p4user@$node "mkdir -p $p4path" ; nohup scp $localfullpath $p4user@$node:$p4path  ) &
      #ssh $p4user@$node "mkdir -p $p4path" ; scp $localfullpath $p4user@$node:$p4path
      rsh -l $p4user $node "mkdir -p $p4path" ; rcp $localfullpath $node:$p4path
      scpproclist[$nodenum]=$!
  fi
  if [ $nodenum -eq 1 ]
      then
      numprocs=$(($numprocs + $p4np + 1))
  else
      numprocs=$(($numprocs + $p4np))
  fi
  echo $node $p4np $p4prog $p4user >> $newpgfile
  echo "done with machine #$nodenum named $node : $p4np added to get $numprocs total procs"
  nodenum=$(($nodenum+1))
done
cat $newpgfile
# execute one call

#echo $$
#echo $!
#echo $PPID
#sleep 132352135135

#( nohup scp bigone jon@bh1: ) &
#scpproclist[4]=$!

#dolr='$'
#for (( i=1 ; i < $numprocs ; i++ ))
#do
#  a=0
#  while [ -n "$a" ]
#  do
#    a=`ps -co pid,ppid | awk '($1=='${scpproclist[$i]}')&&($2=='$$') {print}'`
#    sleep .01
#    echo $i ${scpproclist[$i]} $a
#    done
#done
#exit
#while `ps -co ppid,cmd | grep scp | awk '{print $1}' | grep "$$"`

#do
#  sleep .01
#  echo "a"
#done


##################################
# CHOOSE METHOD OF STARTING BELOW:

# FAST MANY NODES
#((rsh -l $masteruser $masternode "mkdir -p $p4path" ; rcp $newpgfile $masteruser@$masternode:$masterpath ; rsh -l $masteruser $masternode "cd $masterpath; $p4prog $a -p4pg $newpgfile -p4wd $masterpath" 2>&1 ) &)

# SLOW MANY NODES
#((ssh -l $masteruser $masternode "mkdir -p $p4path" ; scp $newpgfile $masteruser@$masternode:$masterpath ; ssh -l $masteruser $masternode "cd $masterpath; $p4prog $a -p4pg $newpgfile -p4wd $masterpath" 2>&1 ) &)

# FAST ONE NODE MANY PROCS
echo "$p4prog $a -p4pg $newpgfile -p4wd $masterpath"
$p4prog $a -p4pg $newpgfile -p4wd $masterpath

 
