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

    module unload python/2.6
    module load python/2.7.1
    export MKL_DYNAMIC=FALSE
    export MKL_NUM_THREADS=8
    export MPLCONFIGDIR=/lustre/medusa/$USER/matplotlibdir/
    unset MPLCONFIGDIR


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

# imagemagick
    module load imagemagick/6.6.9-1
# GSL
    module load cue-gsl
    module load gsl
# below not created for some reason on kraken
    export GSL_LIB=${GSL_DIR}/lib/
    export GSLLIB=${GSL_DIR}/lib/
    # below is so harm will use gsl-config to create GSLCFLAGS
    export PATH=$PATH:${GSL_DIR}/bin/

    module unload python/2.6
    module load python/2.7.1-cnl # has to be specifically this 2.7 version or else compute node complains that not compiled for compute nodes (because this library is  on lustre as required)
    export MKL_DYNAMIC=FALSE
    export MKL_NUM_THREADS=8
    export MPLCONFIGDIR=/lustre/scratch/$USER/matplotlibdir/

    #export MPICH_UNEX_BUFFER_SIZE=1073741824
    #export MPICH_UNEX_BUFFER_SIZE=1073741824
    #export MPICH_MAX_SHORT_MSG_SIZE=32000
    #export MPICH_PTL_UNEX_EVENTS=100000
    #export MPICH_PTL_OTHER_EVENTS=400000
    # prevent unexpected event queue from being exhausted in any situation, but may hurt performance.
    #export MPICH_PTL_SEND_CREDITS=-1
    
    #export MPICH_PTL_SEND_CREDITS=-1

    export MPI_TYPE_DEPTH=20
    export MPICH_MAX_SHORT_MSG_SIZE=16000
    export MPICH_PTL_UNEX_EVENTS=80000
    export MPICH_UNEX_BUFFER_SIZE=768M
fi

export PYTHONPATH=$HOME/py:$PYTHONPATH

# causes infinite loop, so just had them change my primary group
#newgrp tug1111
