for fil in `cat pg.lst`; do echo $fil; ssh $fil "ls -al `pwd`/*fail*"; done
