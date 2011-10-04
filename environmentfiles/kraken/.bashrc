# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions

source ~/.bashrc.jon


module unload PrgEnv-pgi/2.2.74
module load PrgEnv-intel/1.0.0
