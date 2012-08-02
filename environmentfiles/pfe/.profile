# $Header: /cvsroot/bcfg2/bcfg2/Cfg/etc/skel_NAS/.profile/.profile,v 1.1 2009/12/11 16:05:13 dtalcott Exp $
# These commands are executed on a login or start of a PBS job.

# First, run the NAS standard setup.

if [ -e /usr/local/lib/init/global.profile ]; then
	. /usr/local/lib/init/global.profile
fi

# Add your commands here to extend your PATH, etc.

PATH=$PATH:$HOME/bin			# Add private commands to PATH

source ~/.bashrc.jon

# PFE modules
module load comp-intel/2012.0.032
#module load comp-intel/2011.2
#module load mpi-sgi/mpt.2.04.10789
module load mpi-sgi/mpt.2.06r6
#module load mpi-sgi/mpt.2.06a67
#module load mpi-intel/4.0.2.003
export MPICC_CC=icc
export MPIF90_F90=ifort
export MPICXX_CXX=icpc

# other modules
module load petsc/3.1-p7/intel/mpt
module load gnuplot/4.4.0
module load math/intel_mkl_64_10.0.011
module load mathematica/7.0.1
module load matlab/2010b


export MPI_TYPE_DEPTH=20
export MPICH_MAX_SHORT_MSG_SIZE=16000
export MPICH_PTL_UNEX_EVENTS=80000
export MPICH_UNEX_BUFFER_SIZE=768M


