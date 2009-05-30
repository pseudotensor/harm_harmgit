#!/bin/bash
for fil in `/bin/ls *.c`
do
  echo $fil
  diff $fil ./old/ > $fil.diff
done
diff defs.h ./old/ > defs.h.diff
diff decs.h ./old/ > decs.h.diff
diff mpidefs.h ./old/ > mpidefs.h.diff
diff mpidecs.h ./old/ > mpidecs.h.diff
diff global.h ./old/ > global.h.diff
diff mympi.h ./old/ > mympi.h.diff

