#!/bin/bash/

# submitting dependency list
# http://beige.ucs.indiana.edu/I590/node45.html

filereport="radsubmitreport.dat"

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
#JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0a`
# if job already completed, then skip up to that point, and remove dependency at that point
#JOBID=3950778
#echo "$JOBID for batch.qsub.kraken.rada0a" >> $filereport

#JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0b`
# removed dependency because got done with "rada0a" while setting up these files.
JOBID=`qsub batch.qsub.kraken.rada0b`
echo "$JOBID for batch.qsub.kraken.rada0b" >> $filereport
                         
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0c`
echo "$JOBID for batch.qsub.kraken.rada0c" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0d`
echo "$JOBID for batch.qsub.kraken.rada0d" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0e`
echo "$JOBID for batch.qsub.kraken.rada0e" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0f`
echo "$JOBID for batch.qsub.kraken.rada0f" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0g`
echo "$JOBID for batch.qsub.kraken.rada0g" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0h`
echo "$JOBID for batch.qsub.kraken.rada0h" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0i`
echo "$JOBID for batch.qsub.kraken.rada0i" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0j`
echo "$JOBID for batch.qsub.kraken.rada0j" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0k`
echo "$JOBID for batch.qsub.kraken.rada0k" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0l`
echo "$JOBID for batch.qsub.kraken.rada0l" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0m`
echo "$JOBID for batch.qsub.kraken.rada0m" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0n`
echo "$JOBID for batch.qsub.kraken.rada0n" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0o`
echo "$JOBID for batch.qsub.kraken.rada0o" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0p`
echo "$JOBID for batch.qsub.kraken.rada0p" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0q`
echo "$JOBID for batch.qsub.kraken.rada0q" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0r`
echo "$JOBID for batch.qsub.kraken.rada0r" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0s`
echo "$JOBID for batch.qsub.kraken.rada0s" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0t`
echo "$JOBID for batch.qsub.kraken.rada0t" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0u`
echo "$JOBID for batch.qsub.kraken.rada0u" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0v`
echo "$JOBID for batch.qsub.kraken.rada0v" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0w`
echo "$JOBID for batch.qsub.kraken.rada0w" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0x`
echo "$JOBID for batch.qsub.kraken.rada0x" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0y`
echo "$JOBID for batch.qsub.kraken.rada0y" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0z`
echo "$JOBID for batch.qsub.kraken.rada0z" >> $filereport
	    
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nexta`
echo "$JOBID for batch.qsub.kraken.rada0nexta" >> $filereport

JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextb`
echo "$JOBID for batch.qsub.kraken.rada0nextb" >> $filereport

JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextc`
echo "$JOBID for batch.qsub.kraken.rada0nextc" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextd`
echo "$JOBID for batch.qsub.kraken.rada0nextd" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nexte`
echo "$JOBID for batch.qsub.kraken.rada0nexte" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextf`
echo "$JOBID for batch.qsub.kraken.rada0nextf" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextg`
echo "$JOBID for batch.qsub.kraken.rada0nextg" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nexth`
echo "$JOBID for batch.qsub.kraken.rada0nexth" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nexti`
echo "$JOBID for batch.qsub.kraken.rada0nexti" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextj`
echo "$JOBID for batch.qsub.kraken.rada0nextj" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextk`
echo "$JOBID for batch.qsub.kraken.rada0nextk" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextl`
echo "$JOBID for batch.qsub.kraken.rada0nextl" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextm`
echo "$JOBID for batch.qsub.kraken.rada0nextm" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextn`
echo "$JOBID for batch.qsub.kraken.rada0nextn" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nexto`
echo "$JOBID for batch.qsub.kraken.rada0nexto" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextp`
echo "$JOBID for batch.qsub.kraken.rada0nextp" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextq`
echo "$JOBID for batch.qsub.kraken.rada0nextq" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextr`
echo "$JOBID for batch.qsub.kraken.rada0nextr" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nexts`
echo "$JOBID for batch.qsub.kraken.rada0nexts" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextt`
echo "$JOBID for batch.qsub.kraken.rada0nextt" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextu`
echo "$JOBID for batch.qsub.kraken.rada0nextu" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextv`
echo "$JOBID for batch.qsub.kraken.rada0nextv" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextw`
echo "$JOBID for batch.qsub.kraken.rada0nextw" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextx`
echo "$JOBID for batch.qsub.kraken.rada0nextx" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nexty`
echo "$JOBID for batch.qsub.kraken.rada0nexty" >> $filereport
JOBID=`qsub -W depend=afterany:$JOBID batch.qsub.kraken.rada0nextz`
echo "$JOBID for batch.qsub.kraken.rada0nextz" >> $filereport



exit 0



