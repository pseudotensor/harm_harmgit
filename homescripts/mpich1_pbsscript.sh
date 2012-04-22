Script:

#PBS -l nodes=4:compute:ppn=2,walltime=10:00:00
#PBS -N xhpl

# prog name
PROG=xhpl.$(uname -m)
PROGARGS=""

NODES=$PBS_NODEFILE

# How many proc do I have?
NP=$(wc -l $NODES | awk '{print $1}')

# create pgfile with rank 0 node with one less
# process because it gets one by default
ME=$(hostname -s)
N=$(egrep "^$ME\$" $NODES | wc -l | awk '{print $1}')
N=$(($N - 1))
if [ "$N" = "0" ]
then
        >pgfile
else
        echo "$ME $N $PWD/$PROG" >pgfile
fi

# add other nodes to pgfile
for i in $(cat $NODES | egrep -v "^$ME\$" | sort | uniq)
do
        N=$(egrep "^$i\$" $NODES | wc -l | awk '{print $1}')
        ARCH=$(ssh $i uname -m)
        echo "$i $N $PWD/xhpl.$ARCH"
done >>pgfile

# MPICH path
# mpirun is a script, no worries
MPICH=/usr/local/mpich/1.2.6..13/gm/x86_64/smp/pgi64/ssh/bin
PATH=$MPICH/bin:$PATH

export LD_LIBRARY_PATH=/usr/local/goto/lib

set -x

# cd into the directory where I typed qsub
if [ "$PBS_ENVIRONMENT" = "PBS_INTERACTIVE" ]
then
        mpirun.ch_gm \
                -v \
                -pg pgfile \
                --gm-kill 5 \
                --gm-no-shmem \
                LD_LIBRARY_PATH=/usr/local/goto/lib \
                $PROG $PROGARGS
else
        cd $PBS_O_WORKDIR
        cat $PBS_NODEFILE >hpl.$PBS_JOBID

        mpirun.ch_gm \
                -pg pgfile \
                --gm-kill 5 \
                --gm-no-shmem \
                LD_LIBRARY_PATH=/usr/local/goto/lib \
                $PROG $PROGARGS >>hpl.$PBS_JOBID
fi

exit 0
