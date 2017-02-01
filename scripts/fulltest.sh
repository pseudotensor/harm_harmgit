#!/bin/bash

########
### Code to run through complete end-to-end test of harm or harmrad
########



########
# assume already done this outside this script:
# git clone git@github.com:pseudotensor/harm_harmgit.git
# cd harm_harmgit/

########
# now can script actions
# get version setup for master
git checkout master

# default for master branch:
# for non-radiation problem
ln -s initboundcode/init.fishmon.c init.c ; ln -s initboundcode/init.fishmon.h init.h ; ln -s initboundcode/bounds.fishmon.c bounds.c

# for radiation problem
#ln -s initboundcode/init.koral.c init.c ; ln -s initboundcode/init.koral.h init.h ; ln -s initboundcode/bounds.koral.c bounds.c

# compile
make superclean ; make prep ; make -j 16

# default for master branch:
#mkdir run ; cd run ; cp ../grmhd . ; ./grmhd 1 1 1 1 # mpi and openmp both on
mkdir run ; cd run ; cp ../grmhd . ; nohup mpirun -np 16 ./grmhd 1 4 4 1 # mpi and openmp both on

# default for koralinsert or auto branches
# run with 32 openmp threads (efficient use of all 16 cores on 1 node)
#mkdir run ; cd run ; cp ../grmhd . ; nohup ./grmhd 32 1 1 1 &
#mkdir run ; cd run ; cp ../grmhd . ; nohup mpirun -np 1 ./grmhd 1 1 1 1 &

# run for about an hour to get at least a few files.
# this doesn't end the job, just pauses for 1 hour to get to analysis
sleep 1h

##################################################
##################################################


#========================
#Lately-used Branch Details:
#========================

#koralinsert: Latest working (non-devel) for harmrad: b89a412aa9988dde020ebcba011d0b50af7548b2 .

#koralinsert: Latest devel branch for harmrad: HEAD

#master: Latest working (non-devel) for harm3d: HEAD

#subedd: From non-devel harmrad, a setup for sub-Eddington type simulation.

#superedd: From non-devel harmrad, a setup for super-Eddington type simulation.

#========================
#Common problems for master branch:
#========================

# 1) By default, I've assumed the "mpicc" installation handles the
# argument -cc, but some installations don't.  If your compiler doesn't
# take -cc arguments, then your mpicc compiler probably only compiles
# with "gcc" effectively.  Then in the file "makefile", you can set
# USEMCCSWITCH=0 and USEMCCSWITCHFORGCC=0 (first instances of those two
# variables just inside the USEMPI section).

# 2) If you have issues with OpenMP, set USEOPENMP=0 in makehead.inc .  Then when running master branch, do for 1 MPI task:

# mkdir run ; cd run ; cp ../grmhd . ; nohup mpirun -np 1 ./grmhd 1 1 1 &

# 3) If you have issues with MPI but OpenMP is ok, then set USEMPI=0 in makehead.inc .  Then to run do for 16 cores on a node:

# mkdir run ; cd run ; cp ../grmhd . ; ./grmhd 16 &

# 4) If you have issues with MPI and OpenMP, then set USEOPENMP=0 and USEMPI=0 in makehead.inc .  Then to run do:

# mkdir run ; cd run ; cp ../grmhd . ; ./grmhd &


# basic analysis, follows:
# See koralinsert branch ./docs/quick_start_guide/* for tutorials on various aspects and utilities
# See koralinsert branch ./docs/general_plotting_guide.txt for tutorial on data analysis
# See koralinsert branch ./docs/stampede.txt
# See koralinsert branch ./docs/usingpfe.txt

##################################################
##################################################


# basic analysis steps as example:

# move "run" to specific name
mv run mhd1

# ensure you have python utilities
mypath=`pwd`
DIRECTORY=~/py/
if [ -d "$DIRECTORY" ]; then
    echo "directory $DIRECTORY already exists, not overwritting"
else
    echo "directory $DIRECTORY does not exist, creating"
    cd ~ ; git clone git@github.com:pseudotensor/harm_pythontools.git ; mv harm_pythontools py ; cd py ; git branch jon
    cd $mypath
fi

cp ~/py/scripts/makeallmovie.sh .
# edit makeallmovie.sh and change line with dircollect="..." and dirruns="..." so that "..." is "mhd1"
# edit makeallmovie.sh and set numkeep to 50 for test
cat makeallmovie.sh | sed 's/^dirruns=.*/dirruns=\"mhd1\"/g' | sed 's/^numkeep=.*/numkeep=50/g' > makeallmovielocal.sh

# assuming on single node type job (default if system not detected), then do:
# setup links, copy files, and makeavg step that makes avg2d_*.npy.
# then merges to avg2d.npy (average used for rest of analysis)
export moviename="movie1"
bash ./makeallmovielocal.sh ${moviename} 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0
# check if mhd1/$moviename/avg2d.npy created

# do rest of analysis
bash ./makeallmovielocal.sh ${moviename} 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1
# should generate analysis results in mhd1/$moviename/*.png *.txt etc.

# can make movie by doing (new2 for avconv and new3 for ffmpeg)
cd mhd1/$moviename/
cp ~/py/scripts/makelinkimagenew2.sh .
bash ./makelinkimagenew2.sh
# should generate *.mp4 movies from frames
# can use smplayer or other things to view mp4's


# see "system" in makeallmovielocal.sh and how used in ~/py/scripts/makemovie.sh to see what kinds of systems already setup and how you can setup your own system by example.

# common issues on single node systems:
# 1) dash incomplete, need to use bash to run scripts.  Best to remove dash in favor of bash
