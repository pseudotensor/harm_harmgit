#!/bin/bash/

# submitting dependency list
# http://beige.ucs.indiana.edu/I590/node45.html

filereport="radsubmitreport.rada0.94.dat"

#FIRST=`qsub first_1.sh`
#echo $FIRST
#SECOND=`qsub -W depend=afterok:$FIRST second_1.sh`
#echo $SECOND
#THIRD=`qsub -W depend=afterok:$SECOND third_1.sh`
#echo $THIRD
#FOURTH=`qsub -W depend=afterok:$THIRD fourth_1.sh`
#echo $FOURTH
#exit 0

# first job id can come from qstat | grep jmckinne if not doing dep list for very first run (usually want to ensure first run looks good before setting up full dep list)
JOBID=`qsub batch.qsub.kraken.rada0.94a`
# if job already completed, then skip up to that point, and remove dependency at that point
#JOBID=3950778
echo "$JOBID for batch.qsub.kraken.rada0.94a" >> $filereport

JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94b`
# removed dependency because got done with "rada0.94a" while setting up these files.
#JOBID=`qsub batch.qsub.kraken.rada0.94b`
echo "$JOBID for batch.qsub.kraken.rada0.94b" >> $filereport
                         
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94c`
echo "$JOBID for batch.qsub.kraken.rada0.94c" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94d`
echo "$JOBID for batch.qsub.kraken.rada0.94d" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94e`
echo "$JOBID for batch.qsub.kraken.rada0.94e" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94f`
echo "$JOBID for batch.qsub.kraken.rada0.94f" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94g`
echo "$JOBID for batch.qsub.kraken.rada0.94g" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94h`
echo "$JOBID for batch.qsub.kraken.rada0.94h" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94i`
echo "$JOBID for batch.qsub.kraken.rada0.94i" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94j`
echo "$JOBID for batch.qsub.kraken.rada0.94j" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94k`
echo "$JOBID for batch.qsub.kraken.rada0.94k" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94l`
echo "$JOBID for batch.qsub.kraken.rada0.94l" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94m`
echo "$JOBID for batch.qsub.kraken.rada0.94m" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94n`
echo "$JOBID for batch.qsub.kraken.rada0.94n" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94o`
echo "$JOBID for batch.qsub.kraken.rada0.94o" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94p`
echo "$JOBID for batch.qsub.kraken.rada0.94p" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94q`
echo "$JOBID for batch.qsub.kraken.rada0.94q" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94r`
echo "$JOBID for batch.qsub.kraken.rada0.94r" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94s`
echo "$JOBID for batch.qsub.kraken.rada0.94s" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94t`
echo "$JOBID for batch.qsub.kraken.rada0.94t" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94u`
echo "$JOBID for batch.qsub.kraken.rada0.94u" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94v`
echo "$JOBID for batch.qsub.kraken.rada0.94v" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94w`
echo "$JOBID for batch.qsub.kraken.rada0.94w" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94x`
echo "$JOBID for batch.qsub.kraken.rada0.94x" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94y`
echo "$JOBID for batch.qsub.kraken.rada0.94y" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94z`
echo "$JOBID for batch.qsub.kraken.rada0.94z" >> $filereport
	    
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nexta`
echo "$JOBID for batch.qsub.kraken.rada0.94nexta" >> $filereport

JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextb`
echo "$JOBID for batch.qsub.kraken.rada0.94nextb" >> $filereport

JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextc`
echo "$JOBID for batch.qsub.kraken.rada0.94nextc" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextd`
echo "$JOBID for batch.qsub.kraken.rada0.94nextd" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nexte`
echo "$JOBID for batch.qsub.kraken.rada0.94nexte" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextf`
echo "$JOBID for batch.qsub.kraken.rada0.94nextf" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextg`
echo "$JOBID for batch.qsub.kraken.rada0.94nextg" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nexth`
echo "$JOBID for batch.qsub.kraken.rada0.94nexth" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nexti`
echo "$JOBID for batch.qsub.kraken.rada0.94nexti" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextj`
echo "$JOBID for batch.qsub.kraken.rada0.94nextj" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextk`
echo "$JOBID for batch.qsub.kraken.rada0.94nextk" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextl`
echo "$JOBID for batch.qsub.kraken.rada0.94nextl" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextm`
echo "$JOBID for batch.qsub.kraken.rada0.94nextm" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextn`
echo "$JOBID for batch.qsub.kraken.rada0.94nextn" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nexto`
echo "$JOBID for batch.qsub.kraken.rada0.94nexto" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextp`
echo "$JOBID for batch.qsub.kraken.rada0.94nextp" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextq`
echo "$JOBID for batch.qsub.kraken.rada0.94nextq" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextr`
echo "$JOBID for batch.qsub.kraken.rada0.94nextr" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nexts`
echo "$JOBID for batch.qsub.kraken.rada0.94nexts" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextt`
echo "$JOBID for batch.qsub.kraken.rada0.94nextt" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextu`
echo "$JOBID for batch.qsub.kraken.rada0.94nextu" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextv`
echo "$JOBID for batch.qsub.kraken.rada0.94nextv" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextw`
echo "$JOBID for batch.qsub.kraken.rada0.94nextw" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextx`
echo "$JOBID for batch.qsub.kraken.rada0.94nextx" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nexty`
echo "$JOBID for batch.qsub.kraken.rada0.94nexty" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0.94nextz`
echo "$JOBID for batch.qsub.kraken.rada0.94nextz" >> $filereport



exit 0



