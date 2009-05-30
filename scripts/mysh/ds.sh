#!/bin/bash
for fil in `/bin/ls *.c`
do
  echo $fil
  diff $fil ../grmhd-new-1cpu/ > $fil.diff
done
diff defs.h ../grmhd-new-1cpu/ > defs.h.diff
diff decs.h ../grmhd-new-1cpu/ > decs.h.diff

