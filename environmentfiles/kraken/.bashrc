# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions

source ~/.bashrc.jon

# Jon's choice to use ICC:
module unload PrgEnv-pgi/2.2.74
module load PrgEnv-intel/1.0.0

# Sasha's choice to use PathScale:
# module swap PrgEnv-pgi/2.2.48B PrgEnv-pathscale/2.2.48B 2>/dev/null >/dev/null
# To use the latest MPI libraries from Cray:
# module swap xt-mpt/5.0.0 xt-mpt/5.1.0 2>/dev/null >/dev/null
# RECOMPILE!
