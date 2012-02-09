7 Additional optimization on Kraken (few percent speedup, but maybe could be sped up more)
I found a way to optimally place MPI tasks on nodes so that closest MPI tasks are grouped together on the same node.  This, unfortunately, speeds up the code only by a few percent, so see if its worth the effort of figuring out how to use it.

NOTE: While I split up cores between nodes such that each node gets a group of close neighbors, I do not optimize the mapping between groups and nodes.  This might be important for performance and is worth looking into in the future.  

NOTE 2: There are many specialized MPI options that in principle can speed some codes up (but slow other codes down) that we can try.  I have tried some of them, but without luck.  See http://www.nics.tennessee.edu/user-support/mpi-tips-for-cray-xt5

In order to generate optimal MPI rank placement, I wrote a program, which I called "mpiplace", that accepts the (i) CPU geometry, (ii) the geometry of CPUs on a node and produces a file used by Cray's MPI.  

Here is an example of Kraken submit script that makes use this technique (in this example I use 288 cores in a 18 x 8 x 2 geometry and "ask" mpiplace to arrange MPI tasks in tiles of 3 x 2 x 2 -- this groups MPI tasks into groups of 12 so that each group fits on its own node):

#!/bin/bash
#PBS -A TG-ASTxxxxxx
#PBS -l size=288,walltime=24:00:00
#PBS -q small
#PBS -N h2fullpd4nosrcdtuthcour08_pathscale_sse3_mpt51_mpiplace_pd2_0_0_0

cd $PBS_O_WORKDIR

#Generate optimal rank placement (syntax: ./mpiplace Ncpu1 Ncpu2 Ncpu3 Ncorespercpu1 Ncorespercpu2 Ncorespercpu3  >MPICH_RANK_ORDER 2>mpiplace.info)
./mpiplace 18 8 2 3 2 2 >MPICH_RANK_ORDER 2>mpiplace.info

#Set MPI environment variables
export MPICH_RANK_REORDER_METHOD=3      #=3 tells MPI to use the generated rank placement
#export MPICH_MAX_SHORT_MSG_SIZE=512000 #messages shorted than this size (in bytes) are sent faster; increasing this doesn't seem to speed things up
#export MPICH_UNEX_BUFFER_SIZE=2M       #buffer size for storing unclaimed (by MPI_irecv()) messages; 2Mb is too short, so leave at default (~60Mb)
#export MPICH_PTL_UNEX_EVENTS=20        #no. of allowed queued events.  did not play with it.
#export MPICH_ENV_DISPLAY=1             #prints out env info

#Run the MPI program (uses MPICH_RANK_ORDER file if present)
aprun -n 288 ./grmhd 18 8 2

You can downlaod the source of the program, mpiplace.c, here:

http://www.cfa.harvard.edu/~atchekho/research/mpiplace.c

You can compile it as usual (compiler does not matter since it is not run with MPI):  

gcc -O3 -o mpiplace mpiplace.c

8. What else could be sped up 
With the new grid, pre-computations of grid Jacobian/connection are taking quite a bit of time: 1 hour for a 16^3 tile size.  (This is because the new grid has more involved expressions.  This, however, only affects the initialization when the information is precomputed and stored for future use.)  Given that code run time on Kraken is limited to 24 hours, these initial pre-computations eat up about 1/24 = 4% of total runtime.  A simple way to speed such computations up is to make use of axisymmetry of the grid and compute the connection, etc. only for one meridional slice.  This would bring 4% down to <0.4%.  This is a plan for future.

However, right now you can go ahead and make use of the grid.  Please let me know if you have any questions or if something does not work for you.

Cheers,
Sasha

