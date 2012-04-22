# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions

source ~/.bashrc.jon

module unload python/2.6
module load python/2.7.1
export MKL_DYNAMIC=FALSE
export MKL_NUM_THREADS=8

module load globus


module load cue-login-env
if [ $CUE_HOST_PROMPT = nautilus ]; then
    module load PE-intel
    export PYTHONPATH=$HOME/py:$PYTHONPATH
fi

if [ $CUE_HOST_PROMPT = kraken ]; then
#module unload PrgEnv-pgi/2.2.74
#module load PrgEnv-intel/1.0.0
#module unload pgi/11.9.0
#module load intel/12.1.2.273
    module unload PrgEnv-pgi/3.1.72
    module load PrgEnv-intel/3.1.72
    export PYTHONPATH=$HOME/py:$PYTHONPATH
fi


