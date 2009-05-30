#!/bin/bash
# uses an old rdump(or ending-bad rdump, just need header format), an old dump used to restart from, and a new name for the rdump
#
#e.g. dump2rdump.sh rdump.old dump.old rdump.new

cat $1 | head -1 > $3
sed -e '1d' < $2 | awk '{print $7" "$8" "$9" "$10" "$11" "$12" "$13" "$14}' >> $3

# then use dump header to fix new header's:
# t, nstep, dump_cnt,image_cnt,rdump_cnt,dt,realnstep,restartsteps[0,1]
