#!/bin/bash
#  Sample Batch Script for a MVAPICH-Intel job
#
# http://www.loni.org/teragrid/users_guide.php
#
# REMEMBER for HARM code:
# 1) PRODUCTION 1 in init.h
# 2) Set N1,N2,N3, and MAXBND in init.h
# 3) set lim[] in init.c
# 4) mpd &
#
#
# qsub <this filename>
#
#
# $HOME/.soft contains:
#
#  @teragrid-basic
#  @globus-4.0
#  @teragrid-dev
#
# $HOME/.mpd.conf contains:
#
#  MPD_SECRETWORD=XXXXXXX     # random alphanumeric chars
#                             # (MUST contain at least one alphabetic char)
#
# (make sure the file .mpd.conf has permissions 700)
#
#  Submit this script using the command: qsub <script_name>
#
#  Use the "qstat" command to check the status of a job.
#
# The following are embedded QSUB options. The syntax is #PBS (the # does
# _not_  denote that the lines are commented out so do not remove).
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=24:00:00
#
# workq, checkpt, preempt, priority (up to 128*8 procs)
#PBS -q medium
#
#PBS -A TG-PHY120005
#
# nodes: number of 8-core nodes
#   ppn: how many cores per node to use (1 through 8)
#       (you are always charged for the entire node)
#PBS -l size=4608
#
# export all my environment variables to the job
#####PBS -V
#
# job name (default = name of script file)
#PBS -N sa99t1.5708b
#
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
#PBS -o sashaa99t1.5708b.out
#
# filename for standard error (default = <job_name>.e<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
#PBS -e sashaa99t1.5708b.err
#
#PBS -m be
#
#PBS -M jmckinne@stanford.edu
#
#
# End of embedded QSUB options
#
# set echo               # echo commands before execution; use for debugging
#

#cd $PBS_O_WORKDIR
cd /lustre/scratch/rblandfo/
          
###################
#
# setup run
#
##################
#
export NCPUX1=18
export NCPUX2=16
export NCPUX3=16
export NTOT=4608
export FILENAME="sashanoreset"
export DIRFILE="/nics/d/home/rblandfo/sasha"
export RUNDIR=/lustre/scratch/rblandfo/sashaa99t1.5708b/


#############
echo "ncpux1 $NCPUX1"
echo "ncpux2 $NCPUX2"
echo "ncpux3 $NCPUX3"
echo "ntot $NTOT"
echo "filename $FILENAME"
echo "dirfile $DIRFILE"
echo "rundir $RUNDIR"
############################
#
# rest shouldn't change
#
###############################
export BEFOREDIR=`pwd`
mkdir -p $RUNDIR
cd $RUNDIR

mkdir -p dumps
mkdir -p images

# get executable and input files from mass storage
cp $DIRFILE/$FILENAME .
cp $DIRFILE/*.dat .
cp $DIRFILE/*.txt .

wait
chmod u+x $FILENAME
wait

##################################################
#
#    start job
#
##################################################


date

aprun -n $NTOT $RUNDIR/$FILENAME $NCPUX1 $NCPUX2 $NCPUX3 1 0

date
