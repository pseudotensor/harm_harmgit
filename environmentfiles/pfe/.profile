# $Header: /cvsroot/bcfg2/bcfg2/Cfg/etc/skel_NAS/.profile/.profile,v 1.1 2009/12/11 16:05:13 dtalcott Exp $
# These commands are executed on a login or start of a PBS job.

# First, run the NAS standard setup.

if [ -e /usr/local/lib/init/global.profile ]; then
	. /usr/local/lib/init/global.profile
fi

# Add your commands here to extend your PATH, etc.

PATH=$HOME/bin/$PATH			# Add private commands to PATH

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
module load math/intel_cmkl_64_9.1.023
module load mathematica/7.0.1
module load matlab/2010b
module load ffmpeg/1.2
module load gnuplot/4.6.3
#module load python/2.7.3
module load tcl-tk/8.5.11
module load hdf5/1.8.3/intel/mpt
module load imagemagick/6.4.0-3

export MPI_TYPE_DEPTH=20
export MPICH_MAX_SHORT_MSG_SIZE=16000
export MPICH_PTL_UNEX_EVENTS=80000
export MPICH_UNEX_BUFFER_SIZE=768M



module load git/1.7.7.4


# if using user python, then do below:

SRCDIR=/nobackup/jmckinn2/tarballs/
BASE=/nobackup/jmckinn2/
export PYTHONPATH=$BASE/lib/python/:$BASE/py/:$HOME/py/
export PATH=$BASE/bin:$PATH
export PYTHON_LIB=$BASE/lib/
export PYTHON_INC=$BASE/include/python2.7/
export LD_LIBRARY_PATH=$BASE/lib/:$LD_LIBRARY_PATH
export LIBRARY_PATH=$BASE/lib/:$LIBRARY_PATH


module load texlive/2008
