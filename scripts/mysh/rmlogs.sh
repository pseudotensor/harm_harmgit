for fil in `cat pg.lst`; do echo $fil; ssh $fil "rm `pwd`/0_*"; done
