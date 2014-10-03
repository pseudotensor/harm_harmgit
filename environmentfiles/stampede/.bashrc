# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions

source ~/.bashrc.jon


module load globus
#module load intel
module load python
export MKL_DYNAMIC=FALSE
export MKL_NUM_THREADS=16
export MPLCONFIGDIR=$SCRATCH/matplotlibdir/
unset MPLCONFIGDIR
export MPI_TYPE_DEPTH=20 #for using ROMIO
export MPI_TYPE_MAX=65536

export MKL_MIC_ENABLE=1
# -1 means auto-load balancing
export MKL_MIC_WORKDIVISION0=-1 # HOST
export MKL_MIC_WORKDIVISION1=-1 # PHI COPROC1
export MKL_MIC_WORKDIVISION2=-1 # PHI COPROC2
export MKL_MIC_WORKDIVISION=1 # auto-load
    #module load imagemagick/6.6.9-1
# GSL
module load gsl
# below not created for some reason on kraken
export GSL_LIB=${GSL_DIR}/lib/
export GSLLIB=${GSL_DIR}/lib/
    # below is so harm will use gsl-config to create GSLCFLAGS
export PATH=$PATH:${GSL_DIR}/bin/


export MPI_TYPE_DEPTH=20
export MPICH_MAX_SHORT_MSG_SIZE=16000
export MPICH_PTL_UNEX_EVENTS=80000
export MPICH_UNEX_BUFFER_SIZE=768M


export PYTHONPATH=$HOME/py:$PYTHONPATH

# causes infinite loop, so just had them change my primary group
#newgrp tug1111
