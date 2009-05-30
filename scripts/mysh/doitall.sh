mydir=`pwd`
myprog=`basename $mydir`
nx=$1
ny=$2
rm *.o *.il ; make ; mv grmhd $myprog ; sh kill.sh $myprog ; sh rmdat.sh ; ./mympirun.sh 1 pg.lst $myprog $nx $ny > grmhdoutput.txt



