#!/bin/bash

# put script.sh in home directory or somewhere not in run directory
# then run in the home directory:
# nohup sh script.sh <run directory> &
# e.g. nohup sh script.sh myrundir &


# set where storing on relativity
RELSTORAGE=/mnt/removable2/dummy/
RELUSER=dummy
RELNAME=relativity.cfa.harvard.edu

# first argument of script is run directory where done.txt will be generated.
RUNDIR=$1

# sleep interval must be longer than check interval for first step to make sense
# normal
CHECKINTERVALMINS=5
SLEEPINTERVALMINS=10m

# test
#CHECKINTERVALMINS=1
#SLEEPINTERVALMINS=1s

# make run dir on cluster if not there
mkdir -p $RUNDIR

# set directory on relativity
RELDIR=$RELSTORAGE/$RUNDIR/

# make directory on relativity if not there
ssh -x ${RELUSER}@${RELNAME} "mkdir -p $RELDIR"

# go to run directory on cluster
cd $RUNDIR

# set initial value of indicator of whether on very last pass
VERYLASTPASS="0"

# add some regular files to exclusion list.  These will be copied but NOT removed.
# directories themselves (not files in directories) are always excluded.
# multi-CPU files are automatically excluded from deletion in below procedure
EXLIST="./grmhd ./coordparms.dat ./0_logdtfull.out ./probe.dat ./lumvsr.out ./gener.out ./generjet4.out ./generjet3.out ./generjet2.out ./generjet1.out ./generjet0.out ./flener.out ./flenerother5.out ./flenerother4.out ./flenerother3.out ./flenerother2.out ./flenerother1.out ./ener.out ./enerother5.out ./enerother4.out ./enerother3.out ./enerother2.out ./enerother1.out ./debug.out ./nohup.out ./done.txt ./0_log.out ./0_logfull.out ./0_logdt.out ./0_fail.out ./numcpus.txt"

# make remote directories needed (needed for scp to work)
ssh -x ${RELUSER}@${RELNAME} "cd $RELDIR ; mkdir -p images dumps"


# LOOP OVER SLEEP/COPY INTERVALS
while [ 1 ]
  do

# sleep for a while before checking if any new files
  date
  echo "Sleeping"
  sleep ${SLEEPINTERVALMINS}

#######################
# get current file list on cluster in a sorted way so diff doesn't get confused
# clean files before making list
  rm -rf cksumcluster.txt filelistcluster.txt difffilelist.txt diffcksum.txt cksumrel.txt filelistrel.txt
  find -type f -mmin +$CHECKINTERVALMINS -printf "%C@ && %p\n" | sort -n +1 | sed -e 's/.*&& //'  > filelistcluster.txt

######################
# get full sorted list on relativity (clean files before making list)
  ssh -x ${RELUSER}@${RELNAME} "cd $RELDIR ; rm -rf cksumrel.txt filelistrel.txt ; find -type f -mmin +0 -printf '%C@ && %p\n' | sort -n +1 | sed -e 's/.*&& //'  > filelistrel.txt"
  scp -rp ${RELUSER}@${RELNAME}:$RELDIR/filelistrel.txt .

###############################
# check if cluster has new files compared to relativity
  diff -Bbd filelistcluster.txt filelistrel.txt | grep "<" | sed 's/< //' > difffilelist.txt

# GODMARK: Could also add check to see if files already exist and after checksumming them can delete on cluster.  For now assume any file on relativity got deleted from cluster.

  if [ -s difffilelist.txt ]
      then
    # form exclusion list
      DIFFLIST=`cat difffilelist.txt`
      RMLIST=""
      COPYLIST=""
      for fil in $DIFFLIST
	do
	if [ $VERYLASTPASS -eq "0" ]
	    then
	    rm -rf testex.txt
	    echo $EXLIST > exlist.txt
	    for filex in $EXLIST
	      do
            # account for those files with CPU extensions foo.####
            # This doesn't affect other files
	      prefa=`echo $fil | sed "s/\.[0-9][0-9][0-9][0-9]/ /"`
	      pref=`echo $prefa | awk '{print $1}'`
	      grep "$pref" exlist.txt >> testex.txt
	    #echo $fil $filex $pref
	    done
	    rm -rf exlist.txt
	    if [ -s testex.txt ]
		then
	  # then if in exclusion list don't add to removal list
		touch testex.txt
	    else
	  # add to removal list since not on exclusion list
		RMLIST="$RMLIST $fil"
		COPYLIST="$COPYLIST $fil"
	    fi
      # finally remove temporary file
	    rm -rf testex.txt
	else
	  # add all files to remove list when VERYLASTPASS==1
	    RMLIST="$RMLIST $fil"
	    COPYLIST="$COPYLIST $fil"
	fi
      done
  fi



##################################
# copy over new files
##################################
  DIDCOPY=0
  for fil in $COPYLIST
    do
    echo "Copying new file $fil"
# scp sucks at copying files inside directories since it won't create or put files in a remote directory if providing full path locally.  Have to add directory extension
    REMOTEDIR=`dirname $fil`
    FILETOCOPY=`basename $fil`
    scp -rp $fil ${RELUSER}@${RELNAME}:$RELDIR/$REMOTEDIR/$FILETOCOPY
    DIDCOPY=1
  done


##################################
# If copied files, check that copied files are the same
##################################
  if [ $DIDCOPY -eq "1" ]
      then
######################
# rather than getting full file list, only checksum copied files
# cksum only operates on files, not directories (warns about them)
      ssh -x ${RELUSER}@${RELNAME} "cd $RELDIR ; rm -rf cksumrel.txt ; for fil in "$COPYLIST" ; do cksum \$fil >> cksumrel.txt ; done"
# get remote checksum
      scp -rp ${RELUSER}@${RELNAME}:$RELDIR/cksumrel.txt .
# get cluster checksum
      rm -rf cksumcluster.txt ; for fil in $COPYLIST ; do cksum $fil >> cksumcluster.txt ; done
# difference list of checksums
      diff -Bbd cksumcluster.txt cksumrel.txt > diffcksum.txt
# check if any differences
      if [ -s diffcksum.txt ]
	  then
	  echo "Checksums were different! Don't remove any new cluster files, remove relativity files that tried to update, and try again in next iteration"
	  CHECKSUMWAS="0"
	  ssh -x ${RELUSER}@${RELNAME} "cd $RELDIR ; for fil in "$COPYLIST" ; do rm -rf \$fil ; done"
      else
	  
# Then good comparison, so remove local files
#Only remove things NOT in exclusion list
	  for fil in $RMLIST
	    do
	    echo "Removing new file $fil"
	    rm -rf $fil
	  done
	  CHECKSUMWAS="1"
      fi # end if good checksum
  else # end if files to copy
      date
      echo "NO New files to copy"
    #cat difffilelist.txt
      CHECKSUMWAS="1" # indicates ok to be done if no new files
  fi




##################################
# Check if final pass
##################################
  if [ $VERYLASTPASS -eq "1" ]
      then
      date
      echo "Done copying files"
      break
  fi



#########################################
# check if done.txt has been generated by HARM, which says done with simulation
# only do this check if cksum succeeded
#########################################
  if [ $CHECKSUMWAS -eq "1" ]
      then
      if [ -e done.txt ]
	  then
	  date
	  echo "End of computation so copy script ending after one more pass to copy final files"
	  echo "Also will now copy and remove EXLIST files"
	  CHECKINTERVALMINS=0
	  SLEEPINTERVALMINS=1s
	  VERYLASTPASS="1"
      fi
  fi
# end if checksum was good


done
