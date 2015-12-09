#!/bin/bash/

# submitting dependency list
# http://beige.ucs.indiana.edu/I590/node45.html

filereport="radsubmitreport.radtma0.8.8.dat"
rm -rf $filereport

#FIRST=`sbatch first_1.sh`
#echo $FIRST
#SECOND=`sbatch -W depend=afterok:$FIRST second_1.sh`
#echo $SECOND
#THIRD=`sbatch -W depend=afterok:$SECOND third_1.sh`
#echo $THIRD
#FOURTH=`sbatch -W depend=afterok:$THIRD fourth_1.sh`
#echo $FOURTH
#exit 0

# first job id can come from qstat | grep jmckinne if not doing dep list for very first run (usually want to ensure first run looks good before setting up full dep list)
#JOBID=`sbatch batch.sbatch.kraken.radtma0.8a`
# if job already completed, then skip up to that point, and remove dependency at that point
JOBID=4915427
echo "$JOBID for batch.slurm.kraken.radtma0.8a" >> $filereport

JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8b | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8b" >> $filereport
JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8c | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8c" >> $filereport
JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8d | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8d" >> $filereport
JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8e | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8e" >> $filereport
JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8f | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8f" >> $filereport
JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8g | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8g" >> $filereport
JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8h | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8h" >> $filereport
JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8i | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8i" >> $filereport
JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8j | tail -1 | awk '{print $4}'`
echo "$JOBID for batch.slurm.kraken.radtma0.8j" >> $filereport

#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8k | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8k" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8l | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8l" >> $filereport
#JOBID=`sbatch batch.slurm.kraken.radtma0.8m | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8m" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8n | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8n" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8o | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8o" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8p | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8p" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8q | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8q" >> $filereport
#JOBID=`sbatch batch.slurm.kraken.radtma0.8r | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8r" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8s | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8s" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8t | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8t" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8u | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8u" >> $filereport
#JOBID=`sbatch batch.slurm.kraken.radtma0.8v | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8v" >> $filereport
#JOBID=`sbatch batch.slurm.kraken.radtma0.8w | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8w" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8x | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8x" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8y | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8y" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8z | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8z" >> $filereport
	    
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nexta | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nexta" >> $filereport
#
#JOBID=`sbatch batch.slurm.kraken.radtma0.8nextb | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextb" >> $filereport
#
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextc | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextc" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextd | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextd" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nexte | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nexte" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextf | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextf" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextg | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextg" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nexth | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nexth" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nexti | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nexti" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextj | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextj" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextk | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextk" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextl | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextl" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextm | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextm" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextn | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextn" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nexto | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nexto" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextp | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextp" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextq | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextq" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextr | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextr" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nexts | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nexts" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextt | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextt" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextu | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextu" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextv | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextv" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextw | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextw" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextx | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextx" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nexty | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nexty" >> $filereport
#JOBID=`sbatch --dependency=afterany:$JOBID batch.slurm.kraken.radtma0.8nextz | tail -1 | awk '{print $4}'`
#echo "$JOBID for batch.slurm.kraken.radtma0.8nextz" >> $filereport
#


exit 0



