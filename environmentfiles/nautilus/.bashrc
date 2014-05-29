# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions

source ~/.bashrc.jon


module load globus


module load cue-login-env
if [ $CUE_HOST_PROMPT = nautilus ]; then
    module load PE-intel

    module load texlive/2010

    module unload python/2.6
    module load python/2.7.1
    export MKL_DYNAMIC=FALSE
    export MKL_NUM_THREADS=8
    export MPLCONFIGDIR=/lustre/medusa/jmckinne/matplotlibdir/
    unset MPLCONFIGDIR

    # Sasha determined required for many cores to run
    export MPI_TYPE_DEPTH=20 #for using ROMIO
    export MPI_TYPE_MAX=65536

    

fi


if [ $CUE_HOST_PROMPT = kraken ]; then
#module unload PrgEnv-pgi/2.2.74
#module load PrgEnv-intel/1.0.0
#module unload pgi/11.9.0
#module load intel/12.1.2.273
    module unload PrgEnv-gnu/3.1.72
    module unload PrgEnv-pgi/3.1.72
    module load PrgEnv-intel/3.1.72

    # sasha uses below
    #module swap xt-mpich2 xt-mpt
    
    module unload python/2.6
    module load python/2.7.1-cnl # has to be specifically this 2.7 version or else compute node complains that not compiled for compute nodes (because this library is  on lustre as required)
    export MKL_DYNAMIC=FALSE
    export MKL_NUM_THREADS=8
    export MPLCONFIGDIR=/lustre/scratch/jmckinne/matplotlibdir/

#see https://mailman.ucsd.edu/pipermail/enzo-users-l/Week-of-Mon-20110103/000539.html
#Try increasing the value of env var MPICH_UNEX_BUFFER_SIZE (curvalue is 62914560), and/or reducing the size of MPICH_MAX_SHORT_MSG_SIZE (cur value is 128000).
    
    export MPICH_UNEX_BUFFER_SIZE=251658240
    export MPICH_MAX_SHORT_MSG_SIZE=64000
#[356] MPICH PtlEQPoll error (PTL_EQ_DROPPED): An event was dropped on the UNEX EQ handle.  #Try increasing the value of env var MPICH_PTL_UNEX_EVENTS (cur size is 20480).                                    
    export MPICH_PTL_UNEX_EVENTS=100000
    export MPICH_PTL_OTHER_EVENTS=400000
    
    
fi

export PYTHONPATH=$HOME/py:$PYTHONPATH

# just ask NICS to change primary group, else globus doesn't work.
#newgrp tug1111



